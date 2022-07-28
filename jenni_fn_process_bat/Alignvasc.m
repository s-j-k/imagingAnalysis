animal = 'RoofBuddy1';
brainside = 'LIC';%'RIC';

datatable = readtable(['N:\Jenni\' animal '\Zstackinfopersite.xlsx']);


Correspondingsitesind = find(cellfun(@(x) strcmp(x,brainside),datatable.Side));

if strcmp(animal,'RoofBuddy2')
    fixed = imread(['N:\Jenni\' animal '\z-stack\20210810_vascref.JPG']);
    Mostanteriorsite = 10;%num2str([datatable.Site(Correspondingsitesind)])
elseif strcmp(animal,'RoofBuddy1')
    fixed = imread(['N:\Jenni\' animal '\z-stack\20210130_vascref.JPG']);
    Mostanteriorsite = 9;%9;%num2str([datatable.Site(Correspondingsitesind)])
    if strcmp(brainside,'LIC')% remove site 1 that used different stims
        Correspondingsitesind(find(Correspondingsitesind)==1) = [];
    end
elseif strcmp(animal,'Gray2')
    fixed = imread(['N:\Jenni\' animal '\z-stack\20220203_vascref.JPG']);
    Mostanteriorsite = 1;
end
fixed = uint16(imadjust(fixed(:,:,1)));

Mostanteriorind = find(datatable.Site==Mostanteriorsite);
Correspondingsitesind(find(datatable.Site(Correspondingsitesind)==Mostanteriorsite)) = [];


h = figure(); hold on;
imagesc(flipud(fixed)); colormap gray; axis tight

vascim = imread(['N:\Jenni\' animal '\' datatable.VascImage{Mostanteriorind} '.PNG']);%(uint16(fliplr(slice_green{zstackimagecorr})));
vascim = squeeze(vascim(:,:,1));
if strcmp(animal,'RoofBuddy2')
    vascim2 = imread(['N:\Jenni\' animal '\081121_Vasc_LIC_Site1_Start.PNG']);% site 10(uint16(fliplr(slice_green{zstackimagecorr})));
elseif  strcmp(animal,'RoofBuddy1')
    vascim2 = imread(['N:\Jenni\' animal '\020321_Site2_Vasc_RIC_End.PNG']);%site 21 (uint16(fliplr(slice_green{zstackimagecorr})));
elseif strcmp(animal,'Gray2')
    vascim2 = imread(['N:\Jenni\' animal '\020322_LIC_Site3_Vasc_End.PNG']);%site 21 (uint16(fliplr(slice_green{zstackimagecorr})));
end
vascim2 = squeeze(vascim2(:,:,1));
Widefieldsize = [1210,900]; %x is (690px/0.57 microns) and y is 512/0.57 if microns to pix are the same in widefield (But seems like it from the size estomation of the surgery pic)
Widefieldpix = [size(vascim2,1),size(vascim2,2)];
Vascpix = [size(fixed,1),size(fixed,2)];% rough estimate
Vascsize = [1210*4,900*2.5];
Vascratio = Vascsize./fliplr(Vascpix);
meanmeanimage9 = imresize(vascim, 1/4, 'bilinear'); % was 4 for the other bats
meanmeanimage92 = imresize(vascim2, 1/4, 'bilinear');

moving = (uint16(fliplr((meanmeanimage9))));
moving2 = (uint16(fliplr((meanmeanimage92))));
if strcmp(animal,'RoofBuddy2')
    imagesc([1:size(moving,2)]+550,[1:size(moving,1)]+160,flipud(moving),'AlphaData',0.4)
    imagesc([1:size(moving2,2)]+175,[1:size(moving2,1)]+130,flipud(moving2),'AlphaData',0.4)
elseif  strcmp(animal,'RoofBuddy1')
    imagesc([1:size(moving,2)]+180,[1:size(moving,1)]+260,flipud(moving),'AlphaData',0.4)
    imagesc([1:size(moving2,2)]+620,[1:size(moving2,1)]+160,flipud(moving2),'AlphaData',0.4)
elseif  strcmp(animal,'Gray2')
    imagesc([1:size(moving,2)]+900,[1:size(moving,1)]+580,flipud(moving),'AlphaData',0.4)
    imagesc([1:size(moving2,2)]+230,[1:size(moving2,1)]+580,flipud(moving2),'AlphaData',0.4)
end


set(gca,'xtick',[])
set(gca,'ytick',[])

plot([size(fixed,2)/2,size(fixed,2)/2],[1,size(fixed,1)],'k-')
if strcmp(animal,'RoofBuddy2')
    RICrefcoord = [(550-size(fixed,2)/2)*Vascratio(1),0]; % first coordinate from midline
    LICrefcoord = [((size(fixed,2)/2)-175)*Vascratio(1),(160-130)*Vascratio(2)]; % first coordinate from midline second fom RIC most anterior site
elseif  strcmp(animal,'RoofBuddy1')
    LICrefcoord = [((size(fixed,2)/2)-180)*Vascratio(1),0]; % first coordinate from midline
    RICrefcoord = [(620-size(fixed,2)/2)*Vascratio(1),(260-160)*Vascratio(2)];% first coordinate from midline second fom RIC most anterior site
elseif  strcmp(animal,'Gray2')
    RICrefcoord = [(900-size(fixed,2)/2)*Vascratio(1),0]; % first coordinate from midline
    LICrefcoord = [((size(fixed,2)/2)-230)*Vascratio(1),(580-580)*Vascratio(2)]; % first coordinate from midline second fom RIC most anterior site

end

% save(['N:\Jenni\' animal '\z-stack\RICrefcoord.mat'],'RICrefcoord')
% save(['N:\Jenni\' animal '\z-stack\LICrefcoord.mat'],'LICrefcoord')
print(['N:\Jenni\' animal '\z-stack\Vascoverlayonsurgerypic.pdf'],'-dpdf')
print(['N:\Jenni\' animal '\z-stack\Vascoverlayonsurgerypic.svg'],'-dsvg')


% try for all sites


figure
imshowpair(imadjust(fixed), imadjust(moving))
[optimizer, metric] = imregconfig('multimodal')
optimizer.InitialRadius = 6.25e-30;
optimizer.Epsilon = 1.5e-10;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 500;
tform = imregtform(imadjust(moving),imadjust(fixed),'similarity',optimizer,metric);
warp=imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
figure; hold on
imshowpair(fixed,warp)


moving = imhistmatch(moving,fixed);

% Estimate the transformation needed to bring the two images into alignment.

[~,movingReg] = imregdemons(moving,fixed,[500 400 200],...
    'AccumulatedFieldSmoothing',1.3);
imshowpair(fixed,movingReg)
% figure
%         imshowpair(fixed, moving2)
%         [optimizer, metric] = imregconfig('multimodal')
%         optimizer.InitialRadius = 6.25e-3;
%         optimizer.Epsilon = 1.5e-10;
%         optimizer.GrowthFactor = 1.01;
%         optimizer.MaximumIterations = 300;
%         tform = imregtform(moving2,fixed,'affine',optimizer,metric);
%         warp=imwarp(moving2,tform,'OutputView',imref2d(size(fixed)));
%         figure
%         imshowpair(fixed,warp)
