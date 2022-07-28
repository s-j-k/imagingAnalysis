%lines to edit: graystackimage, ensure gray image is the only image in working directory
clear clc close all
datafolder = 'V:\Aaron\aw019\z-stack';
graystack='test_gray_reg.tif';
align=0; %1 for stackreg, 0 if no alignment desired
invert=0; %1 if threshold inversion needed, 0 if not needed
cd(datafolder)

%%
tic
graystackinfo=imfinfo(graystack);
slicenb=size(graystackinfo,1);


disk='E:';
javaaddpath([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\mij.jar'])
javaaddpath([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\ij-1.52i.jar'])
MIJ.start([disk '\KishoreLab\Shared\Matlab\Fiji.app'])

graystackpath=strcat('path=[',datafolder,'\',graystack,']');
insertAfter(graystackpath,"\","\")
MIJ.run('Open...', graystackpath);


MIJ.run('Enhance Contrast...', 'saturated=0.3 process_all');

if align==1
i=1;
while i < slicenb;
     MIJ.run('Next Slice [>]');
     i=i+2;
end
  
MIJ.run('StackReg ', 'transformation=[Rigid Body]');
else
end

MIJ.run('Gaussian Blur...', 'sigma=2 stack');
MIJ.run('8-bit');
MIJ.run('Auto Local Threshold...', 'method=Otsu radius=15 parameter_1=0 parameter_2=0 white stack')

if invert==1;
MIJ.run('Invert', 'stack');
else
end

MIJ.run('Analyze Particles...', 'size=80-Infinity circularity=0.80-1.00 add stack');

ij.IJ.runMacro(strcat('newImage("Untitled", "16-bit black", 796, 512, ', num2str(slicenb),');'));
ij.IJ.runMacro('roiManager("Show All");');
ij.IJ.runMacro('roiManager("Fill");');
fijipath=insertBefore(path,"\","\");
savepath=strcat('saveAs("Tiff", "', fijipath, '\\rawplaque.tif");');
ij.IJ.runMacro(savepath)
toc