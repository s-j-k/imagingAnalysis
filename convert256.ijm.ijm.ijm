run("Image Sequence...", 
"open=[D:/Research/Scott/raw data/20210124/TSeries-01242021-0952-003/TSeries-01242021-0952-003_Cycle00001_Ch2_000001.ome.tif] sort use");
run("Z Project...", "projection=[Average Intensity]");
selectWindow("TSeries-01242021-0952-003");
run("Size...", "width=256 height=256 depth=20000 constrain average interpolation=Bilinear");
saveAs("Tiff", "D:/Research/Scott/raw data/20210124/TSeries3/TSeries-01242021-256-003.tif");
selectWindow("AVG_TSeries-01242021-0952-003");
saveAs("Tiff", "D:/Research/Scott/raw data/20210124/AVG_TSeries-01242021-003.tif");
close();
close();
selectWindow("TSeries-01242021-256-002.tif");
close();