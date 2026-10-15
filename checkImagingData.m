% check the imaging data
close all
animalDay='D7'; % checked D1, D2, D3, D4, D5
A = fn_readBinSuite2p(['O:\rockfish\su\Data\sk318_reg\' animalDay '\suite2p\plane0\'],'data.bin',1000);
% implay(Anew)
Anew=A/prctile(A(:),99.0);
figure; imagesc(mean(Anew,3));
colormap(gray);
caxis([0 1]);
cd(['O:\rockfish\su\Data\sk318_reg\' animalDay '\suite2p\plane0']);
saveas(gcf, [animalDay '.png']);
saveas(gcf, [animalDay '.fig']);
% tells you between each frame the offset used to correct for motion
figure;subplot(1,2,1); title('Motion correction X offset'); hold on;
plot(ops.xoff);
subplot(1,2,2); plot(ops.yoff);title('Motion correction Y offset');
saveas(gcf, [animalDay '_offset.fig']);
% values of >20 are bad

% for sk318 ROI size is 656
% 775 1463