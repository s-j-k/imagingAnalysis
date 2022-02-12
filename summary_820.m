function summary_820
% you need to change directory to the one that has the .mat files
cd '/Users/bbs/Dropbox/Kinematic_Data/Bruker/raw data/0104data/gcamp'
load('DATA820.mat','ISO_F')
load('data920.mat','F_dff')
% remove inactive cells
 cells=ones(181,1);
 remove=[24; 45; 47; 65; 115; 154; 106; 105; 120; 123; 129; 133;144;146;148;149;165;181];
 cells(remove)=0;
 realcells=find(cells);
time=(0:1800)/30;
% plot the data
figure;
subplot(2,3,2)
plot(time,ISO_F(58,1000:2800))
ylim([-.25 2.5])
xlim([0 1800/30])
title('Single Cell (820nm)')
ylabel('DF/F')
xlabel('Seconds')

subplot(2,3,1)
plot(time,F_dff(58,1000:2800))
xlim([0 1800/30])
ylim([-.25 2.5])
title('Single Cell (920nm)')
ylabel('DF/F')
xlabel('Seconds')

subplot(2,3,5)
imagesc(ISO_F(realcells,1000:2800),[0 1])
title('Population (820nm)')
ylabel('Cell')
xlabel('Seconds')

set(gca,'ytick',[1 length(realcells)])
set(gca,'xtick',[1 20*30 40*30 1800])
set(gca,'xticklabel',{'0';'20';'40'; '60'})
set(gca,'tickdir','out')

subplot(2,3,4)
imagesc(F_dff(realcells,1000:2800),[0 1])
ylabel('Cell')
xlabel('Seconds')
set(gca,'ytick',[1 length(realcells)])
set(gca,'xtick',[1 20*30 40*30 1800])
set(gca,'xticklabel',{'0';'20';'40'; '60'})
set(gca,'tickdir','out')
title('Population (920nm)')

subplot(2,3,6)
colorbar
box off
axis off


