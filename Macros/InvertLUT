// Invert LUT for .tif files.
// Johana Itzel Ramos-Galguera

dir = getDirectory("Choose a folder");
list = getFileList(dir);

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".tiff")) {
        path = dir + list[i];

        // Abrir archivo .tif
        open(path);
        wait(200);

        // Gene filenames
        base = substring(list[i], 0, lastIndexOf(list[i], "."));

       
        run("Split Channels");

        // brightfield (grises)
        selectWindow("C1-" + base);
        run("Grays");

        // fluorescence (rojo)
        selectWindow("C2-" + base);
        run("Red");

        // Merge (C1 y C2)
        run("Merge Channels...", "c1=[C1-" + base + "] c2=[C2-" + base + "] c3=None create");

        // Save as .avi
        saveAs("AVI", dir + base + "_merged.avi");
        close(); // merged
        selectWindow("C1-" + base); close();
        selectWindow("C2-" + base); close();
    }
}
