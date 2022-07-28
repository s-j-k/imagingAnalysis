function megabat 
listid = {'K:\Jenni\Gray2\roistructure\megaroiscoordRIC.mat',...
    'N:\Jenni\RoofBuddy2\roistructure\megaroiscoordRIC.mat',...
    'N:\Jenni\RoofBuddy2\roistructure\megaroiscoordLIC.mat',...
    'N:\Jenni\RoofBuddy1\roistructure\megaroiscoordLIC.mat'}
for ii = 1:length(listid)
    data(ii)=load(listid{ii});
end
midlinecoord = [];
antepostcoord = [];
depth = [];
tuning = [];
for ii = 1:length(listid)
    midlinecoord = [midlinecoord data(ii).megaroiscoord.midlinecoord];
    antepostcoord = [antepostcoord data(ii).megaroiscoord.anteropost];
    depth = [depth data(ii).megaroiscoord.depth];
    tuning = [tuning data(ii).megaroiscoord.tuning];
end
tunedcells = find(tuning>0);
brainside = [ones(1,size(data(1).megaroiscoord.midlinecoord,2)+size(data(2).megaroiscoord.midlinecoord,2)),zeros(1,size(data(3).megaroiscoord.midlinecoord,2)+size(data(4).megaroiscoord.midlinecoord,2))];

newtonecolormap = colormaptonotopy(10); % 10 is nb of tones
figure(); hold on;
for i = 1:length(tunedcells)
    plot3(midlinecoord(tunedcells(i)),antepostcoord(tunedcells(i)),-depth(tunedcells(i)),'.','Color',newtonecolormap(tuning(tunedcells(i)),:))
end
ylabel('Relative AP coordinate [\mum]');xlabel(['ML coordinate [\mum]', ' n_c_e_l_l_s' num2str(length(cattuning)), 'n_s_i_t_e_s' num2str(st)]);zlabel('Depth [\mum]');
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_Tonotopy3dcoord.pdf'],'-dpdf')
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_Tonotopy3dcoord.svg'],'-dsvg')


figure(); hold on;
scatter3(midlinecoord(intersect(tunedcells,find(brainside==1))),-(antepostcoord(intersect(tunedcells,find(brainside==1)))),-depth(intersect(tunedcells,find(brainside==1))),5,newtonecolormap(tuning(intersect(tunedcells,find(brainside==1))),:),'filled');
ylabel('AP axis');xlabel('ML axis');zlabel('Depth')
title('All RIC n = 2')

figure(); hold on;
scatter3(midlinecoord(intersect(tunedcells,find(brainside==0)))-max(midlinecoord(intersect(tunedcells,find(brainside==0)))),antepostcoord(intersect(tunedcells,find(brainside==0))),-depth(intersect(tunedcells,find(brainside==0))),5,newtonecolormap(tuning(intersect(tunedcells,find(brainside==0))),:),'filled');
ylabel('AP axis');xlabel('ML axis');zlabel('Depth')
title('All LIC n = 2')
end