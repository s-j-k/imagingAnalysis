% check the imaging data

animalDay='D5'; % checked D1, D2, D3, D4, D5
A = fn_readBinSuite2p(['O:\rockfish\su\Data\sk318_reg\' animalDay '\suite2p\plane0\','data.bin'],1000);
% implay(Anew)
Anew=A/prctile(A(:),99.0);
figure; imagesc(mean(Anew,3));
colormap(gray);
caxis([0 1]);
cd(['O:\rockfish\su\Data\sk318_reg\' animalDay '\suite2p\plane0']);
saveas(gcf, [animalDay '.png']);
