%% generate image-based PSTH and save locally

temp=reshape(img,697,403,50,100); %reshape your matrix into trials
psth=mean(temp,4); %average across all trials
psth16=uint16(psth); %convert double back to uint16 to write images

% you can add a circshift here so that you can have the last 5 frames go
% first

%write images in TIFF to explore in ImageJ
for i=1:size(psth16,3)
imwrite(psth16(:,:,i),['p0_g_wn',num2str(i)],'Tiff','Compression','none');
end;

% in imagej, optimize your image and then save the intensity/mean image

intensity=imread('AVG_plane0_green_white.tif'); %load your intensity/mean image
stim=mean(psth16(:,:,4:6),3); %calculate your stimulus mean
base=mean(psth16(:,:,30:50),3); %calculate your baseline mean
ratio=stim./base; %create a ratio

%pixel-baesd visualization of a ratio

%convert your intensity image to a double and adjust the contrast
intensity = double(squeeze(mean(intensity,3)));
intensity = intensity/max(max(intensity));
intensity = imadjust(intensity);

%set the threshold of your df/F
minval = 0.8;
maxval = 1.3;

%convert ratio into an 8-bit image
ratiouint8=uint8(((ratio-minval)/(maxval-minval))*255);
cmap = jet(256); %placeholder colormap, can adjust as needed

ratiorgb=ind2rgb(ratiouint8,cmap); %convert the indexed image to RGB
ratiohsv=rgb2hsv(ratiorgb); %convert the RGB image to HSV
ratiohsv(:,:,3)=intensity; %map the intensity image to the "value" of HSV
newratio=hsv2rgb(ratiohsv); %new image keeps the color from the RGB image and the intensity from the intensity image

imwrite(newratio,['whitenoise_plane0_green_white.tif']);