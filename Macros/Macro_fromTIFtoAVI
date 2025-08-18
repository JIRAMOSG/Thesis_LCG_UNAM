// This macro converts all .tif files in a selected folder to .avi videos.
// Assumes each TIFF is a stack (e.g. time-lapse or Z-stack).
// Author: Johana Itzel Ramos-Galguera
// Output: .avi files saved in the same directory.

dir = getDirectory("Choose the folder containing TIFF files");
list = getFileList(dir);

setBatchMode(true);

for (i = 0; i < list.length; i++) {
    filename = list[i];

    // Only process .tif 
    if (endsWith(filename, ".tif")) {
        fullPath = dir + filename;
        print("Processing: " + filename);

        open(fullPath);
        title = getTitle();

        // Check that it's a stack
        if (nSlices() < 2) {
            print("File " + filename + " is not a stack. Skipping.");
            close();
            continue;
        }

        // Convert to 8-bit if necessary
        if (bitDepth != 8) {
            run("8-bit");
        }
        run("Grays");

        // Get base name without extension
        dotIndex = lastIndexOf(filename, ".");
        if (dotIndex > 0)
            baseName = substring(filename, 0, dotIndex);
        else
            baseName = filename;

        aviPath = dir + baseName + ".avi";

        // Save as AVI (JPEG compression, 5 fps)
        run("AVI... ", "compresion=JPEG frame=5  save=[" + aviPath + "]");
        print("Saved AVI: " + aviPath);
        close();
    }
}

setBatchMode(false);
print("All conversions completed.");
