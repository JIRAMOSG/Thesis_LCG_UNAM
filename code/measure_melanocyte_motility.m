% This code has been designed to try to measure motility in confocal results, however I recommed use more sophisticated approaches like IMARIS to this aim. 
%  - Loads time-lapse TIFF (single file)
%  - Detects tdTomato+ cells per frame (2D max-projection if z-stack)
%  - Links detections into tracks (simple assignment + gap closing)
%  - Computes motility metrics and saves CSVs/plots
%
% Requirements: Image Processing Toolbox, Statistics Toolbox (for pdist2/hungarian)
% Optional: install munkres() from File Exchange for Hungarian if you don't want built-in.
%
clear; clc;

%% ---------- USER PARAMETERS ----------
tif_path = 'your_timelapse.tif';     % path to single multi-page tiff
channel_idx = 1;                     % if file contains multiple channels in pages, adapt reading
pixel_size_um = 0.33;                % microns per pixel (edit)
time_interval_sec = 60;              % seconds between frames (edit)
max_displacement_um = 10;            % max allowed displacement per frame (um)
max_gap_frames = 2;                  % allow linking across up to this many missing frames
min_track_length_frames = 5;         % filter out very short tracks
output_prefix = 'melanocytes_analysis';
%% -------------------------------------

% Convert units
time_interval_min = time_interval_sec / 60;
max_disp_px = max_displacement_um / pixel_size_um;

%% ------------- Load TIFF stack -------------
% We attempt to read with tifflib via Tiff or tifread
info = imfinfo(tif_path);
nPages = numel(info);

% Try to detect if file has Z per time or just frames:
% If ImageDescription or PageName encodes time/z it could be complex.
% Simple heuristic: if image height/width same for all pages we assume pages are frames
h = info(1).Height; w = info(1).Width;

% Load all pages into a cell array first
fprintf('Reading %d pages from %s ...\n', nPages, tif_path);
stack_pages = cell(nPages,1);
t = Tiff(tif_path,'r');
for p = 1:nPages
    t.setDirectory(p);
    stack_pages{p} = t.read();
end
t.close();

% Heuristic: if images are grouped in (z x time) ordering you'll need to set z per time.
% We'll try to detect constant z: if nPages equals nTime * nZ you'll need to specify nZ.
% For now, ask the user to set nZ if needed (but we won't prompt — assume 2D or already projected)
% To be robust, allow user to set nZ here:
nZ = 1; % EDIT: set to number of z-slices per timepoint if your tif is zstack per time
% If you set nZ>1, ensure nPages is divisible by nZ
if nZ > 1 && mod(nPages,nZ) ~= 0
    error('nPages (%d) not divisible by nZ (%d) — adjust nZ or supply a different file.', nPages, nZ);
end

if nZ==1
    nFrames = nPages;
    % Convert to 3D array (time, y, x)
    sample = stack_pages{1};
    stack = zeros(nFrames, size(sample,1), size(sample,2), 'like', sample);
    for f=1:nFrames
        stack(f,:,:) = stack_pages{f};
    end
else
    nFrames = nPages / nZ;
    % Make max projection per time point
    sample = stack_pages{1};
    stack = zeros(nFrames, size(sample,1), size(sample,2), 'like', sample);
    idx = 1;
    for tt = 1:nFrames
        tmp = zeros(nZ, size(sample,1), size(sample,2), 'like', sample);
        for z = 1:nZ
            tmp(z,:,:) = stack_pages{idx};
            idx = idx + 1;
        end
        stack(tt,:,:) = squeeze(max(tmp,[],1)); % max projection per time
    end
end
fprintf('Loaded %d frames (after optional Z projection)\n', nFrames);

%% --------- Preprocess: bleach correction & background subtraction ----------
% Simple median-based bleach correction: divide by smooth curve of mean intensity
means = zeros(nFrames,1);
for f=1:nFrames
    means(f) = mean(stack(f,:,:),'all');
end
% Fit a smooth curve (lowess) - for bleaching correction we normalize to the 95th percentile
x = (1:nFrames)';
sm = smooth(x, means, 0.1, 'loess'); % adjust smoothing window if needed
baseline = sm;
% Avoid division by zero
baseline(baseline<=0) = 1;
% Apply correction: scale each frame to compensate
stack_corr = stack;
for f=1:nFrames
    stack_corr(f,:,:) = stack(f,:,:) ./ baseline(f) * median(means);
end

% Background subtraction: morphological opening to estimate background
bg_est = zeros(1,size(stack_corr,2),size(stack_corr,3), 'like', stack_corr);
se = strel('disk', 20); % choose radius bigger than cells
for f=1:nFrames
    frame = squeeze(stack_corr(f,:,:));
    bg = imopen(frame, se);
    stack_corr(f,:,:) = frame - bg;
end

%% --------- Detection per frame ----------
detections = cell(nFrames,1);
min_area_px = round((3 / pixel_size_um)^2); % example: ignore objects smaller than ~3um^2 (adjust)
for f=1:nFrames
    frame = squeeze(stack_corr(f,:,:));
    % Enhance contrast
    frame_adap = imadjust(mat2gray(frame));
    % Adaptive threshold (Sauvola/Adaptive)
    T = adaptthresh(frame_adap, 0.45); % sensitivity param: adjust 0.3-0.6
    bw = imbinarize(frame_adap, T);
    % Clean
    bw = bwareaopen(bw, min_area_px);
    bw = imclose(bw, strel('disk',3));
    bw = imfill(bw, 'holes');
    % Separate touching objects: distance/watershed if needed
    D = -bwdist(~bw);
    D(~bw) = -Inf;
    L = watershed(D);
    bw_sep = bw;
    bw_sep(L==0) = 0;

    props = regionprops(bw_sep, frame, 'Centroid','Area','MeanIntensity');
    centroids = reshape([props.Centroid],2,[])';
    areas = [props.Area]';
    meanInt = [props.MeanIntensity]';
    detections{f} = table(centroids(:,1), centroids(:,2), areas, meanInt, ...
                         'VariableNames', {'x','y','area','meanInt'});
end
fprintf('Detection complete: average detections/frame = %.2f\n', mean(cellfun(@height, detections)));

%% --------- Link detections into tracks (frame-to-frame assignment) ----------
% We'll build a track table: columns: trackID, frame, x, y
nextTrackID = 1;
tracks = table([],[],[],[], 'VariableNames', {'trackID','frame','x','y'});

% Keep active tracks: struct array with fields id, last_pos, last_frame, gap_count
active = [];

for f=1:nFrames
    cur = detections{f};
    if isempty(cur)
        % increment gap for active tracks
        for ai = 1:numel(active)
            active(ai).gap = active(ai).gap + 1;
        end
        % remove those with gap > max_gap_frames
        active = active([active.gap] <= max_gap_frames);
        continue;
    end

    cur_pos = [cur.x, cur.y];

    if isempty(active)
        % start new tracks for all detections
        for i=1:size(cur_pos,1)
            tracks = [tracks; {nextTrackID, f, cur_pos(i,1), cur_pos(i,2)}];
            active(end+1).id = nextTrackID; %#ok<SAGROW>
            active(end).pos = cur_pos(i,:);
            active(end).last_frame = f;
            active(end).gap = 0;
            nextTrackID = nextTrackID + 1;
        end
        continue;
    end

    % Build cost matrix between active last positions and current detections
    active_pos = vertcat(active.pos);
    D = pdist2(active_pos, cur_pos); % Euclidean pixel distances
    % Invalidate distances > max_disp_px
    D(D > max_disp_px) = Inf;

    % Hungarian assignment
    % Use munkres algorithm if available; otherwise use built-in assignment via matchpairs
    try
        [assignA, assignB] = matchpairs(D, Inf); % returns pairs (row indices, col indices)
    catch
        % fallback: greedy nearest assignment
        assignA = [];
        assignB = [];
        [r,c] = find(D < Inf);
        if ~isempty(r)
            % sort pairs by distance
            pairs = [r, c, D(sub2ind(size(D), r, c))];
            pairs = sortrows(pairs, 3);
            usedA = false(size(active_pos,1),1);
            usedB = false(size(cur_pos,1),1);
            for k=1:size(pairs,1)
                a = pairs(k,1); b = pairs(k,2);
                if ~usedA(a) && ~usedB(b)
                    assignA(end+1) = a; %#ok<AGROW>
                    assignB(end+1) = b; %#ok<AGROW>
                    usedA(a) = true; usedB(b) = true;
                end
            end
        end
    end

    % Mark assigned active tracks and detections
    assignedActive = false(numel(active),1);
    assignedDet = false(size(cur_pos,1),1);

    for k = 1:size(assignA,1)
        a = assignA(k); b = assignB(k);
        if isfinite(D(a,b))
            tid = active(a).id;
            tracks = [tracks; {tid, f, cur_pos(b,1), cur_pos(b,2)}];
            % update active
            active(a).pos = cur_pos(b,:);
            active(a).last_frame = f;
            active(a).gap = 0;
            assignedActive(a) = true;
            assignedDet(b) = true;
        end
    end

    % For unassigned detections: start new tracks
    for i=1:size(cur_pos,1)
        if ~assignedDet(i)
            tracks = [tracks; {nextTrackID, f, cur_pos(i,1), cur_pos(i,2)}];
            active(end+1).id = nextTrackID;
            active(end).pos = cur_pos(i,:);
            active(end).last_frame = f;
            active(end).gap = 0;
            nextTrackID = nextTrackID + 1;
        end
    end

    % For active tracks that were not assigned, increment gap
    for a = 1:numel(active)
        if ~assignedActive(a)
            active(a).gap = active(a).gap + 1;
        end
    end
    % Remove those with gap > max_gap_frames
    active = active([active.gap] <= max_gap_frames);
end

% Convert tracks to table with numeric arrays
tracks_sorted = sortrows(tracks, {'trackID','frame'});
trackIDs = unique(tracks_sorted.trackID);
fprintf('Total tracks before filtering: %d\n', numel(trackIDs));

%% --------- Compute per-track metrics ----------
% Build a table of per-track metrics
trk_metrics = table('Size',[0 8], 'VariableTypes', ...
    {'double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'trackID','nFrames','pathLength_um','netDisp_um','meanSpeed_um_per_min','directionality','meanIntensity','MSD_alpha'});

for k=1:numel(trackIDs)
    tid = trackIDs(k);
    rows = tracks_sorted(tracks_sorted.trackID==tid, :);
    if height(rows) < min_track_length_frames
        continue
    end
    xs = rows.x; ys = rows.y; frs = rows.frame;
    % compute displacements in pixels
    dx = diff(xs); dy = diff(ys);
    dists_px = sqrt(dx.^2 + dy.^2);
    % times: use frame difference * time_interval_min
    dt_min = diff(frs) * time_interval_min;
    instantaneous_speed = dists_px ./ dt_min; % px per min
    instantaneous_speed_um = instantaneous_speed * pixel_size_um;
    pathLength_um = sum(dists_px) * pixel_size_um;
    netDisp_um = sqrt((xs(end)-xs(1))^2 + (ys(end)-ys(1))^2) * pixel_size_um;
    meanSpeed = mean(instantaneous_speed_um);
    directionality = netDisp_um / pathLength_um;
    % mean intensity: we need to retrieve per-frame mean intensity if available
    % We stored per-frame mean intensity in detections but not in tracks table; approximate by sampling detections
    meanIntensity = NaN;
    % MSD: compute for lags up to half track length, estimate exponent by fitting MSD ~ tau^alpha
    pos = [xs, ys] * pixel_size_um;
    n = size(pos,1);
    maxlag = floor(n/2);
    msd = zeros(maxlag,1);
    taus = (1:maxlag)' * time_interval_min; % minutes
    for L = 1:maxlag
        diffs = pos(1+L:end,:) - pos(1:end-L,:);
        msd(L) = mean(sum(diffs.^2,2));
    end
    % fit log-log slope for alpha, robust only if enough points
    alpha = NaN;
    if maxlag >= 3
        p = polyfit(log(taus), log(msd), 1);
        alpha = p(1);
    end

    trk_metrics = [trk_metrics; {tid, height(rows), pathLength_um, netDisp_um, meanSpeed, directionality, meanIntensity, alpha}];
end

fprintf('Tracks after filtering (length >= %d): %d\n', min_track_length_frames, height(trk_metrics));

%% --------- Save outputs ----------
writetable(tracks_sorted, [output_prefix '_tracks.csv']);
writetable(trk_metrics, [output_prefix '_track_metrics.csv']);
fprintf('Saved %s_tracks.csv and %s_track_metrics.csv\n', output_prefix, output_prefix);

%% --------- Quick plots for QC ----------
figure;
histogram(trk_metrics.meanSpeed, 30);
xlabel('Mean speed (um / min)'); ylabel('count'); title('Distribution of track mean speeds');

figure;
scatter(trk_metrics.pathLength_um, trk_metrics.netDisp_um, 30, trk_metrics.meanSpeed, 'filled');
xlabel('Path length (um)'); ylabel('Net displacement (um)');
colorbar; title('Path length vs net displacement (color=mean speed)');

% show sample tracks overlay on first frame
figure; imshow(mat2gray(squeeze(stack(1,:,:)))); hold on;
colors = lines(min(20, height(trk_metrics)));
cnt=0;
for k=1:min(20, height(trk_metrics))
    tid = trk_metrics.trackID(k);
    rows = tracks_sorted(tracks_sorted.trackID==tid,:);
    plot(rows.x, rows.y, '-', 'LineWidth', 1.5, 'Color', colors(mod(k-1,size(colors,1))+1,:));
    cnt = cnt+1;
end
title(sprintf('Example tracks overlay (first %d tracks)', cnt));

fprintf('Done.\n');
