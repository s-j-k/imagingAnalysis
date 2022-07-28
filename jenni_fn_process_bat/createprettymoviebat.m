function createprettymoviebat()
load('N:\Jenni\RoofBuddy2\000_h5\site1\31.25\TonotopyMeanEvokedMovie.mat')

disk = 'E:';
javaaddpath ([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\mij.jar'])
javaaddpath ([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\ij-1.52i.jar'])
MIJ.start([disk '\KishoreLab\Shared\Matlab\Fiji.app'])
normmovie = cellfun(@(x) (x)./median(median(median(x(100:400,100:400,:)))),MeanEvokedMovie,'UniformOutput',0)%-min(min(min(x))))./max(max(max(x))),MeanEvokedMovie,'UniformOutput',0);
tmpmovieevokedtone = cat(4,normmovie{:});
meanmovie = mean(tmpmovieevokedtone,4);
 
MIJ.createImage(meanmovie); %(100:400,100:400,:)

end