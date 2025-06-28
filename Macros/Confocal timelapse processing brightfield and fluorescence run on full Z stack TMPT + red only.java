//This macro produces maximum intensity projections of fluorescence channels, along with a guassian based stack focussed projection for the brightfield, all in timelapsed files. 
//Input: CZI files in a selected folder
//Process: MIP channels 2,3,4. Best focus projection for brighfield channel 1
//Output: AVI and Tiff to source folder

dir = getDirectory("Choose a folder");
list = getFileList(dir);
//setBatchMode(true);
for (loop=0; loop<lengthOf(list);loop++){
if(endsWith(list[loop],".czi")){

// open image #i
run("Bio-Formats", "open=["+dir+list[loop]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");

Title = getTitle();
//separate channels
run("Split Channels");
selectWindow("C1-"+Title);
rename("TPMT");
//Focussed T-PMT image
run("Gaussian-based stack focuser", "radius_of_gaussian_blur=3");
run("8-bit");
rename("TPMT2");
selectWindow("TPMT");
close();
selectWindow("C2-"+Title);
rename("Red");
run("Z Project...", "projection=[Max Intensity] all");
/*selectWindow("C3-"+Title);
rename("Green");
run("Z Project...", "projection=[Max Intensity] all");
selectWindow("C4-"+Title);
rename("Red");
run("Z Project...", "projection=[Max Intensity] all");*/
//run("Merge Channels...", "c1=MAX_Red c2=MAX_Green c3=MAX_Blue c4=TMPT2 create");
run("Merge Channels...", "c1=MAX_Red c4=TPMT2 create ignore");
saveAs("Tiff", dir+Title+".tif");
filename = dir+Title+".AVI";
run("AVI... ", "compression=None frame=5 save=[filename]");
run("Close All");
}
}