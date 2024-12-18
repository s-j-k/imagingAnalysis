 // this is a macro that will make a grid of ROIs along a selected FOV
 // ctrl + A to select the entire FOV prior to running the macro
 // macro will ask for total number of boxes
 //in the case of a 720 x 512px image where you want 4x4 squares, the total number is then
 // 512/4 = 128 squares
 
 
 outputDir = getDirectory("Choose directory for saving grids");
    	
    	roiManager("reset");

        mainWindow = getTitle();
        close("\\Others"); // close all but selected window
        h = getHeight;
        w = getWidth;
        area = w*h;

        getSelectionCoordinates(x, y);
     
        gridNum = getNumber("How many grid boxes in total?",128);
        boxArea = area/gridNum; // area of each grid box
        boxSide = sqrt(boxArea); // side length of each grid box
        numBoxY = floor(h/boxSide); // number of boxes that will fit along the width
        numBoxX = floor(w/boxSide); // "" height
        remainX = (w - (numBoxX*boxSide))/2; // remainder distance left when all boxes fit
        remainY = (h - (numBoxY*boxSide))/2;
        for (i=0; i<numBoxY; i++) { // draws rectangles in a grid, centred on the X and Y axes and adds to ROI manager
        	for (j=0; j<numBoxX; j++) {
        		makeRectangle((remainX + (j*boxSide)), (remainY+(i*boxSide)), boxSide, boxSide); 
        		roiManager("add");
        	}
        }

//        for (i=0; i<gridNum; i++) { // generate random number to select a square in the grid, makes sure the square has not already been chosen
//        	selectWindow(mainWindow); 
//        	roiManager("select", i);
//        	flag = 0;
        	
//    	    for (j=0; j<x.length; j++) {
//          for (k=0; k<y.length; k++) {
//    	    		if (selectionContains(x[j], y[k])==1) {
//   	    			flag = 1;
//    	    		}
//    	    	}
//    	    }
    	    
//    	    if (flag==1) {
//    	    	run("Duplicate...", "duplicate"); // duplicates the random square in the grid box
//    	        saveAs("Tiff", outputDir + mainWindow + "_Grid" + i);
//    	    }
        
//        }