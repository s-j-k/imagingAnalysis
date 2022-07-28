function [newratio,cmap]=overlayKK (ratio,intensity,counter,cmap,nfactor,xbitdepth)
% takes in a ratio image (or any other indexed image) and then uses the jet
% colormap to pseudocolor the image.  It then converts it into an HSV
% image.  The intensity of the pixels is then defined by another image
% (i.e. brightness of a pixel will be defined by the intensity image and
% the color of the pixel will be defined by the ratio image).
%--------------------------------------------------------------------------

minval = 0.95;
maxval = 1.1;
% 
% %convert ratio into an 8-bit image
ratiouint8=uint8(((ratio-minval)/(maxval-minval))*255);

blueColorMap = [(zeros(1, 132)), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), (zeros(1, 132))];
colorMap = flipud([redColorMap; (zeros(1, 256)); blueColorMap]');

cmap = colorMap;%polarmap(256);%jet(256);

cmap = jet(256);
% load('c:\Users\kuchik01\Dropbox\KishoreDropbox\ActivePassive\colormaps_for_figures\newmap.mat');
% cmap=newmap;
%load the appropriate colormap
% load('Z:\bacskai\Michal\MATLAB programs michal\colormaps\finalcm.mat');
% cmap=cm;


%For the first image in the sequence, assigns the colormap from above to
%the image (to save time afterwards, however, it uses the same colormap for
%all subsequent images)
if counter == 1
    imshow(ratiouint8);
    colormap(cmap);
    cmap=colormap;
    else cmap=cmap;
end 

ratiorgb=ind2rgb(ratiouint8,cmap); %convert the indexed image to RGB
ratiohsv=rgb2hsv(ratiorgb); %convert the RGB image to HSV
% if xbitdepth == 16 %olympus fix because 16 bit images only have 12 bit intensity range
%     xbitdepth=12;
% end;
%intensity=double(imadjust(intensity));
%intensity=intensity/(2^xbitdepth-1); %normalize intensity to maximum possible of yellow + blue channels
ratiohsv(:,:,3)=intensity; %map the intensity image to the "value" of HSV
newratio=hsv2rgb(ratiohsv); %new image keeps the color from the RGB image and the intensity from the intensity image
%disp('done');

    