close all

animal = 'RoofBuddy1';
brainside = 'LIC';
hardriveani = 'N';
datatable = readtable([hardriveani ':\Jenni\' animal '\Zstackinfopersite.xlsx']);


Correspondingsitesind = intersect(find(cellfun(@(x) strcmp(x,brainside),datatable.Side)),find(datatable.Rois==1));
Mostanteriorsite = datatable.Site(intersect(find(datatable.VascRef==1),Correspondingsitesind));

% if strcmp(animal,'RoofBuddy2') % specify the side of the brain too
%     Mostanteriorsite = datatable.Site(intersect(find(datatable.VascRef==1),Correspondingsitesind));
% elseif strcmp(animal,'RoofBuddy1') % specify the side of the brain too
%     Mostanteriorsite = datatable.Site(intersect(find(datatable.VascRef==1),Correspondingsitesind));
%     if strcmp(brainside,'LIC')% remove site 1 that used different stims
% %         Correspondingsitesind(find(Correspondingsitesind)==1) = [];
%     end
% end
% remove ref site from corresponding site 
Mostanteriorind = find(datatable.Site==Mostanteriorsite);
% Correspondingsitesind(Correspondingsitesind==Mostanteriorind) = [];
% Correspondingsitesind(find(datatable.Site(Correspondingsitesind)==10)) = [];
fixed = imread([datatable.datadrive{Mostanteriorind} ':\Jenni\' animal '\' datatable.VascImage{Mostanteriorind} '.PNG']);%(uint16(fliplr(slice_green{zstackimagecorr})));
fixed = squeeze(fixed(:,:,1));
Widefieldsize = [1210,900]; %x is (690px/0.57 microns) and y is 512/0.57 if microns to pix are the same in widefield (But seems like it from the size estomation of the surgery pic)
Widefieldpix = [size(fixed,2),size(fixed,1)];
foursize = [round(750/2.55),round(512/2.66)]; % X AXIS IS CORRECTED
twosize = [round(750/1.32),round(512/1.37)];
ratiowidefieldfour = round(foursize)./round(Widefieldsize);
ratiowidefieldtwo = round(twosize)./round(Widefieldsize);
pixratiofour =  round(Widefieldpix.*ratiowidefieldfour); % size of fov in pixel for widefield
pixratiotwo = round(Widefieldpix.*ratiowidefieldtwo);
roissites = find(datatable.Rois==1);
widefieldratio = (Widefieldpix)./(Widefieldsize); % less pixel i y but mor distance
load([hardriveani ':\Jenni\' animal '\z-stack\' brainside 'refcoord.mat'],[brainside 'refcoord']) % That should be generalized 
col = [1,1,1];
% plot overlaying site put the coord in microns from most anterior site
% recorded and from the midline (see ref coordinate) % for LIC midline
% coord from top left so make an excpetion to remove the site of the site
% from this measure
% 2X is not exactely starting from the center of widefield I think x is
% pretty faithful (see vasculature in end images) for RFB2 LIC 21 and 23
% some ref vasc overlap 
if  strcmp(brainside,'LIC')
yoffset = -37;%-69;%in -163; % in microns (anterior) computed on 2X image chack that it's the same for 4X % closer ot anterior
elseif strcmp(brainside,'RIC')
yoffset = 0;%-69;%in -163; % in microns (anterior) computed on 2X image chack that it's the same for 4X % closer ot anterior
end
if  strcmp(brainside,'LIC') && strcmp(animal,'RoofBuddy2')
xoffset =-30;
elseif strcmp(brainside,'RIC') && strcmp(animal,'RoofBuddy2')
xoffset = -30;
elseif strcmp(brainside,'RIC') && strcmp(animal,'Gray2')
xoffset = -30;
elseif strcmp(brainside,'LIC') && strcmp(animal,'RoofBuddy1')
xoffset =-50;
yoffset = 30;
end
if strcmp(brainside,'LIC')
    frommidline = LICrefcoord(1)-Widefieldsize(1);
    frommostanterior = LICrefcoord(2);
elseif strcmp(brainside,'RIC')
    frommidline = RICrefcoord(1);
    frommostanterior = RICrefcoord(2);
end
h = figure(); hold on;
if  strcmp(brainside,'LIC') % NPO it's a mirror
    imagesc(([1:size(fixed,2)]/widefieldratio(1))+frommidline,([1:size(fixed,1)]/widefieldratio(2))+frommostanterior,(fixed)); colormap gray; axis tight
elseif strcmp(brainside,'RIC')
    imagesc(([1:size(fixed,2)]/widefieldratio(1))+frommidline,([1:size(fixed,1)]/widefieldratio(2))+frommostanterior,fliplr(fixed)); colormap gray; axis tight
end
set(gca,'YDir','reverse')
if  strcmp(brainside,'LIC') % NPO it's a mirror
    set(gca,'XDir','reverse')
end
for st =1:length(Correspondingsitesind)
    clear ytemp xtemp
    moving = imread([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\' datatable.VascImage{Correspondingsitesind(st)} '.PNG']);
    moving = squeeze(moving(:,:,1));
    
    field(st) = load([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(st))) '_fieldnewcoordpx.mat'],'fieldnewcoordpx')
    fov(st) = load([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(st))) '_fovdnewcoordpx.mat'],'fovnewcoordpx')
    if strcmp(brainside,'RIC')
        imagesc((([1:size(fixed,2)]-field(st).fieldnewcoordpx(1))/widefieldratio(1))+frommidline,(([1:size(fixed,1)]+field(st).fieldnewcoordpx(2))/widefieldratio(2))+frommostanterior,fliplr((moving)),'AlphaData',0.5); colormap gray
    elseif strcmp(brainside,'LIC')
        imagesc((([1:size(fixed,2)]-field(st).fieldnewcoordpx(1))/widefieldratio(1))+frommidline,(([1:size(fixed,1)]+field(st).fieldnewcoordpx(2))/widefieldratio(2))+frommostanterior,((moving)),'AlphaData',0.5); colormap gray
    end
end


for st = 1:length(Correspondingsitesind)
    % for RFB2 LIC, vasc picture can not be trusted because it's a start
    % one (forgot to take the end vasc) not an end one I specifically computed the  diff using the 2X 2P
    % image in between site 21 and 23 x = -47 microns(lateral) y = 161
    % microns (more posterior) in 2X px it's x = -62 and y = 221;
    if strcmp(brainside,'LIC') && strcmp(animal, 'RoofBuddy2') && sum(ismember(datatable.Site(Correspondingsitesind(st)),[23,24]))
        field(st).fieldnewcoordpx(1) = -35;% in widefield
        field(st).fieldnewcoordpx(2) = 130;% in widefield
    end
    field(st).fieldnewcoordpx(1)
    clear xvalues yvalues
    xorigine = floor(((Widefieldpix(1)/2)-(fov(st).fovnewcoordpx(3)/2))); % in widefield
    yorigine = floor(((Widefieldpix(2)/2)-(fov(st).fovnewcoordpx(4)/2))); % in widefield most anterior y coord
    
    if  field(st).fieldnewcoordpx(1)>=0 % positive towards the midline and so needs to be removed for microns
        if strcmp(brainside,'LIC')
            lateralfovx(st) = (size(fixed,2)-xorigine)-field(st).fieldnewcoordpx(1)+xoffset;
            midlinefovx(st) = lateralfovx(st)-fov(st).fovnewcoordpx(3);
        elseif strcmp(brainside,'RIC')
            lateralfovx(st) = (size(fixed,2)-xorigine)-field(st).fieldnewcoordpx(1)-xoffset;
            midlinefovx(st) = lateralfovx(st)-fov(st).fovnewcoordpx(3);
        end
    elseif field(st).fieldnewcoordpx(1)<0
        if strcmp(brainside,'LIC')
            lateralfovx(st) = (size(fixed,2)-xorigine)-field(st).fieldnewcoordpx(1)+xoffset;
            midlinefovx(st) = lateralfovx(st)-fov(st).fovnewcoordpx(3);
        elseif strcmp(brainside,'RIC')
            lateralfovx(st) = (size(fixed,2)-xorigine)-field(st).fieldnewcoordpx(1)-xoffset; % off set per zoom
            midlinefovx(st) = lateralfovx(st)-fov(st).fovnewcoordpx(3);
        end
    end
    xvalues = (([lateralfovx(st),lateralfovx(st),midlinefovx(st),midlinefovx(st),lateralfovx(st)]/widefieldratio(1))+frommidline) % size(fixed,2)-xorigine lateral start of fov
    %     elseif field(st).fieldnewcoordpx(1)<0
    %         lateralfovx(st) = (size(fixed,2)-xorigine)-field(st).fieldnewcoordpx(1);
    %         midlinefovx(st) = lateralfovx(st)-fov(st).fovnewcoordpx(3);
    %         xvalues = (([lateralfovx(st),lateralfovx(st),midlinefovx(st),midlinefovx(st),lateralfovx(st)]/widefieldratio(1))+frommidline) % size(fixed,2)-xorigine lateral start of fov
    %     end
    %     yvalues = (([fov(st).fovnewcoordpx(2),(fov(st).fovnewcoordpx(2)+fov(st).fovnewcoordpx(4)),(fov(st).fovnewcoordpx(2)+fov(st).fovnewcoordpx(4)),fov(st).fovnewcoordpx(2),fov(st).fovnewcoordpx(2)])/widefieldratio(2))+frommostanterior;
    anterofovy(st) = (yorigine)+field(st).fieldnewcoordpx(2)+yoffset;
    posteriorfovy(st) = anterofovy(st)+fov(st).fovnewcoordpx(4);
    yvalues = (([anterofovy(st),posteriorfovy(st),posteriorfovy(st),anterofovy(st),anterofovy(st)]/widefieldratio(2))+frommostanterior) % size(fixed,2)-xorigine lateral start of fov
    
    plot(xvalues,yvalues,'--','Color',[1,1,1]);
end
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_SiteoverVac'],'-dpdf')
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_SiteoverVac.svg'],'-dsvg')



h = figure(); hold on;
if  strcmp(brainside,'LIC') % NPO it's a mirror
    imagesc(([1:size(fixed,2)]/widefieldratio(1))+frommidline,([1:size(fixed,1)]/widefieldratio(2))+frommostanterior,(fixed)); colormap gray; axis tight
elseif strcmp(brainside,'RIC')
    imagesc(([1:size(fixed,2)]/widefieldratio(1))+frommidline,([1:size(fixed,1)]/widefieldratio(2))+frommostanterior,fliplr(fixed)); colormap gray; axis tight
end
set(gca,'YDir','reverse')
if  strcmp(brainside,'LIC') % NPO it's a mirror
    set(gca,'XDir','reverse')
end

for st =1:length(Correspondingsitesind)
    clear xvalues yvalues
    moving = imread([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\' datatable.VascImage{Correspondingsitesind(st)} '.PNG']);
    moving = squeeze(moving(:,:,1));
    
    field(st) = load([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(st))) '_fieldnewcoordpx.mat'],'fieldnewcoordpx')
    fov(st) = load([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(st))) '_fovdnewcoordpx.mat'],'fovnewcoordpx')
    
  if strcmp(brainside,'RIC')
        imagesc((([1:size(fixed,2)]-field(st).fieldnewcoordpx(1))/widefieldratio(1))+frommidline,(([1:size(fixed,1)]+field(st).fieldnewcoordpx(2))/widefieldratio(2))+frommostanterior,fliplr((moving)),'AlphaData',0.5); colormap gray
    elseif strcmp(brainside,'LIC')
        imagesc((([1:size(fixed,2)]-field(st).fieldnewcoordpx(1))/widefieldratio(1))+frommidline,(([1:size(fixed,1)]+field(st).fieldnewcoordpx(2))/widefieldratio(2))+frommostanterior,((moving)),'AlphaData',0.5); colormap gray
  end
end
for st =1:length(Correspondingsitesind)
    clear xvalues yvalues
    xvalues = (([lateralfovx(st),lateralfovx(st),midlinefovx(st),midlinefovx(st),lateralfovx(st)]/widefieldratio(1))+frommidline) % size(fixed,2)-xorigine lateral start of fov
    yvalues = (([anterofovy(st),posteriorfovy(st),posteriorfovy(st),anterofovy(st),anterofovy(st)]/widefieldratio(2))+frommostanterior) % size(fixed,2)-xorigine lateral start of fov
    plot(xvalues,yvalues,'--','Color',[1,1,1]);
end

newtonecolormap = colormaptonotopy(10); % 10 is nb of tones
    % plot Rois onto map
    %load roiscoord % transform everything in microns
anterofovymicrons = (anterofovy/widefieldratio(2))+frommostanterior;
midlinefovxmicrons = (midlinefovx/widefieldratio(1))+frommidline;
for st = 1:length(Correspondingsitesind)
    clear newcoord
    load([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\' datatable.subfolder{Correspondingsitesind(st)} '\site' num2str(datatable.Site(Correspondingsitesind(st))) '\31.25\31.25_rois_coord_plane0.mat']);
    
    if strcmp(datatable.zoomRecording{Correspondingsitesind(st)},'2X')
        if  strcmp(brainside,'LIC') % NPO it's a mirror
            newcoord = cellfun(@(x) [(((x(:,1)))/1.32)+midlinefovxmicrons(st),((512-x(:,2))/1.37)+anterofovymicrons(st)],roisCoord{1},'UniformOutput',0);
        elseif   strcmp(brainside,'RIC') % NPO it's a mirror
            newcoord = cellfun(@(x) [(((750-x(:,1)))/1.32)+midlinefovxmicrons(st),((512-x(:,2))/1.37)+anterofovymicrons(st)],roisCoord{1},'UniformOutput',0);
        end
    else strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
         if  strcmp(brainside,'LIC') % NPO it's a mirror
            newcoord = cellfun(@(x) [(((x(:,1)))/1.32)+midlinefovxmicrons(st),((512-x(:,2))/1.37)+anterofovymicrons(st)],roisCoord{1},'UniformOutput',0);
        elseif   strcmp(brainside,'RIC') % NPO it's a mirror
            newcoord = cellfun(@(x) [(((750-x(:,1)))/2.55)+midlinefovxmicrons(st),((512-x(:,2))/2.66)+anterofovymicrons(st)],roisCoord{1},'UniformOutput',0);
        end
%         newcoord = cellfun(@(x) [((750-(x(:,1)))/1.32)+midlinefovxmicrons(st),(x(:,2)/1.37)+anterofovymicrons(st)],roisCoord{1},'UniformOutput',0);
    end
    
    cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',col),newcoord)
    load([datatable.datadrive{1} ':\Jenni\' animal '\z-stack\Tonotopystruct.mat']);
    
    tunedcells = find(Tonotopystruct(datatable.Site(Correspondingsitesind(st))).cellstuning(2,:));
    bftunedcellsind{st} =  Tonotopystruct(datatable.Site(Correspondingsitesind(st))).cellstuning(3,tunedcells);
    
    cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',newtonecolormap(y,:)),newcoord(tunedcells),mat2cell(bftunedcellsind{st},1,ones(1,length(bftunedcellsind{st}))))
  
    Relativecoord = newcoord;
    midlinecoord{st} = cellfun(@(x) median(x(:,1)),Relativecoord(tunedcells));
    anteropostlinecoord{st} = cellfun(@(x) median(x(:,2)),Relativecoord(tunedcells));
    if strcmp(animal,'RoofBuddy1') || strcmp(animal,'RoofBuddy2')
    depthtuning{st} = Tonotopystruct(datatable.Site(Correspondingsitesind((st)))).depth(tunedcells);
    depthtuningall{st} = Tonotopystruct(datatable.Site(Correspondingsitesind((st)))).depth(:);
    elseif strcmp(animal,'Gray2') % add depth manually per notes because animal has a window
    depthtuning{st} =  repmat(datatable.Depth((datatable.Site(Correspondingsitesind((st))))),1,size(tunedcells,2));  
    depthtuningall{st} =  repmat(datatable.Depth((datatable.Site(Correspondingsitesind((st))))),1,size(Tonotopystruct(datatable.Site(Correspondingsitesind(st))).cellstuning(2,:),2));  
    end
% coordinates for every rois
 midlinecoordall{st} = cellfun(@(x) median(x(:,1)),Relativecoord(:));
 anteropostlinecoordall{st} = cellfun(@(x) median(x(:,2)),Relativecoord(:));
end
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_SiteoverVascWithRois.pdf'],'-dpdf')
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_SiteoverVascWithRois.svg'],'-dsvg')



% save every coordinate + tuning + social + sweep maybe ? 
load([hardriveani ':\Jenni\' animal '\roistructure\socialroisdata.mat'])

midlinecoordalltmp = [];
anteropostlinecoordalltmp = [];
tontopytmp = [];
depthtuningtmp = [];
for st = 1:length(Correspondingsitesind)
    midlinecoordalltmp = [midlinecoordalltmp midlinecoordall{st}'];
    anteropostlinecoordalltmp = [anteropostlinecoordalltmp anteropostlinecoordall{st}'];
    tontopytmp = [tontopytmp, Tonotopystruct(datatable.Site(Correspondingsitesind(st))).cellstuning(3,:)];
    depthtuningtmp = [depthtuningtmp depthtuningall{st}'];
%     socialtmp = [socialtmp, cellfun(@median,socialroisdata(st).selectivityindexperstimmed)'];
%     socialtmpshuffle = [socialtmp, cellfun(@median,socialroisdata(st).selectivityindexperstimmed)'];

end
megaroiscoord.midlinecoord = midlinecoordalltmp;
megaroiscoord. anteropost = anteropostlinecoordalltmp;
megaroiscoord.depth = depthtuningtmp;
megaroiscoord.tuning = tontopytmp;
save([hardriveani ':/Jenni/' animal '/roistructure/megaroiscoord' brainside '.mat'],'megaroiscoord')
% mega tuning scatter
cattuning = [bftunedcellsind{:}];
catmidline = [midlinecoord{:}];
catap = [anteropostlinecoord{:}];
[catdepth] = [depthtuning{:}];
uniquedepth = unique(catdepth);
for ii = 1:length(uniquedepth)
    mediantuningdepth(ii) = median(cattuning(find(catdepth==uniquedepth(ii))));
end

figure(); hold on; subplot(1,3,1); plot(cattuning,catmidline,'.'); title('Midline distance')
subplot(1,3,2); plot(cattuning,catap,'.'); title('AP distance')
subplot(1,3,3); plot(cattuning,-catdepth,'.'); title('Depth')
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_Distancetuningsummary'],'-dpdf')

figure(); hold on; subplot(1,3,1); plot(catmidline,cattuning,'.'); title('Midline distance')
subplot(1,3,2); plot(catap,cattuning,'.'); title('AP distance')
subplot(1,3,3); plot(-catdepth,cattuning,'.'); title('Depth')
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_Distancetuningsummary2'],'-dpdf')
figure(); hold on;
for i = 1:length(catap)
    plot3(catmidline(i),catap(i),-catdepth(i),'.','Color',newtonecolormap(cattuning(i),:))
end
ylabel('Relative AP coordinate [\mum]');xlabel(['ML coordinate [\mum]', ' n_c_e_l_l_s' num2str(length(cattuning)), 'n_s_i_t_e_s' num2str(st)]);zlabel('Depth [\mum]');
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_Tonotopy3dcoord.pdf'],'-dpdf')
print([hardriveani ':\Jenni\' animal '\figure\tonotopy\' animal brainside '_Tonotopy3dcoord.svg'],'-dsvg')


figure(); hold on;
scatter3(catmidline,catap,-catdepth,5,newtonecolormap(cattuning,:),'filled');
ylabel('AP axis');xlabel('ML axis');zlabel('Depth')

predictors = [catmidline;catap;catdepth];
tabletest = table(cattuning',catap',catmidline',catdepth','VariableNames',{'Tuning','AP','MD','D'});
modelspec = 'Tuning ~ AP*MD*D - AP:MD:D';
[mdl] = fitglm(tabletest,modelspec);
[a,b,c] = glmfit(predictors',cattuning');



Y=cattuning';
X=predictors';

% get a solution,
X = [X ones(length(X),1)];
B    = inv(X'*X)*X'*Y; % just as above
Yhat = X*B;
Res  = Y - Yhat;

SStotal = norm(Y - mean(Y)).^2;
SSeffect = norm(Yhat - mean(Yhat)).^2;
SSerror  = norm(Res-mean(Res)).^2;

df      = rank(X)-1;
dferror = length(Y) - df - 1;
R2      = SSeffect / SStotal;
F       = (SSeffect / df) / (SSerror / dferror);
p       = 1 - fcdf(F,df,dferror);

% make a figure
figure('Name','Multiple Regression');
subplot(1,2,1); plot(Y,'x','LineWidth',5);  hold on
plot(Yhat,'r','LineWidth',2); grid on;
title(['R^2 ' num2str(R2) ' Data and model F=' num2str(F) ' p=' num2str(p)],'FontSize',14)
subplot(1,2,2); normplot(Res);


% ------------------------------------
% semi-partial correlation coefficient
% ------------------------------------
% let's think of the model and what it means it terms of geometry.
% the data Y can be described as a vector and a point in R_20
% a space with 20 dimensions. We then establish a model X with 5
% regressors; that is we look for a combination of these 5 vectors which
% will get as close as possible to Y. To find the actual contribution of x1
% to the data for this model one needs to look at how much x1 explains
% to the total variance, ie we want to compare the R2 between the full and
% a reduced model without x1 - the difference will be how much x1 explains in Y.

Xreduced    = X(:,2:end); % reduced model all minus 1st regressor 5c hange this to estimate other things
Breduced    = inv(Xreduced'*Xreduced)*Xreduced'*Y;
Yhatreduced = Xreduced*Breduced;
Resreduced  = Y - Yhatreduced;

dfreduced       = rank(Xreduced) -1 ;
dferrorreduced = length(Y) - dfreduced - 1;
SSeffectreduced = norm(Yhatreduced-mean(Yhatreduced)).^2;
SSerrorreduced  = norm(Resreduced-mean(Resreduced)).^2;
R2reduced       = SSeffectreduced / SStotal;
Freduced       = (SSeffectreduced / dfreduced) / (SSerrorreduced / dferrorreduced);
preduced       = 1 - fcdf(Freduced,dfreduced,dferrorreduced);

Semi_Partial_corr_coef = R2 - R2reduced;
dfe_semi_partial_coef  = df - dfreduced;
F_semi_partail_coef    = (Semi_Partial_corr_coef*dferror) / ...  % variance explained by x1
    ((1-R2)*dfe_semi_partial_coef); % unexplained variance overall
p_semi_partial_coef    = 1 - fcdf(Semi_Partial_corr_coef, df, dfe_semi_partial_coef); % note df is from the full model

% make a figure
figure('Name','Multiple Regression - Full versus reduced model');
subplot(2,2,1); plot(Y,'x','LineWidth',5);  hold on
plot(Yhat,'r','LineWidth',2); grid on;
title(['Data and model F=' num2str(F) ' p=' num2str(p)])
subplot(2,2,2); plot(Y,'x','LineWidth',5);  hold on
plot(Yhatreduced,'r','LineWidth',2); grid on;
title(['Data and reduced model F=' num2str(Freduced) ' p=' num2str(preduced)])
subplot(2,2,3); plot(Yhat,'b','LineWidth',2); grid on; hold on
plot(Yhatreduced,'r','LineWidth',2);
title('Modelled data for both models')
subplot(2,2,4); plot((Res-Resreduced).^2,'k','LineWidth',2); grid on; hold on
title('Difference of residuals squared')


% --------------------------------
% partial correlation coefficient
% --------------------------------
% As we shall see below, this is easily obtained using projections - but
% for now let just think about what we want to measure. We are interested
% in knowing the correlation between y and x1 controlling for the effect of
% the other xs, that is removing the effect of other xs. Compred to
% semi-parital coef we also want to remove the effect of xs on y or if you
% prefer we want to compute how much of x1 we need to get to the point
% defined by the xs which is the closest to Y

% above we removed X(:,2:end) from Y and got Resreduced
% we need to do the same for x1 so that we can correlate x1 and y witout xs
x = X(:,1);
B = inv(Xreduced'*Xreduced)*Xreduced'*x;
xhat = Xreduced*B;
Resx = x - xhat;

% the correlation between Resreduced and Resx is the partial coef
Partial_coef = corr(Resx,Resreduced);






























% update this with color and proper boundaries 
load(['N:\Jenni\' animal '\roistructure\socialroisdata.mat'])
% save(['N:\Jenni\' animal '\z-stack\coordpixsite.mat'],'coordpixsite')

minselect = -0.3;%boundaryhist(1);
maxselect = 0.3;%boundaryhist(2);
colormapselect = colormap((redblue(100)));
selectcorrespondance = linspace(-1,1,100);
selectcorrespondance2 = linspace(minselect,maxselect,100);
colorredwhite = flipud([ones(100,1),(log10([1.01:0.1:11])/log10(11))',(log10([1.01:0.1:11])/log10(11))']);%([0:0.01:1-0.01].^1.2)',([0:0.01:1-0.01].^1.05)'];
posselect = find(socialselectivity>0);
negselect = find(socialselectivity<=0);
restselectneg = find(socialselectivity>minselect & socialselectivity<0);
restselectpos = find(socialselectivity<maxselect & socialselectivity>=0);
vectcorrespondanceselect(posselect(find(socialselectivity(posselect)>=maxselect))) = maxselect;
vectcorrespondanceselect(negselect(find(socialselectivity(negselect)<=minselect))) = minselect;
for i = 1:length(restselectneg)
    vectcorrespondanceselect(restselectneg(i)) = selectcorrespondance2(max(find(socialselectivity(restselectneg(i))>=selectcorrespondance2)));
end
for i = 1:length(restselectpos)
    vectcorrespondanceselect(restselectpos(i)) = selectcorrespondance2(max(find(socialselectivity(restselectpos(i))>=selectcorrespondance2)));
end
for ii = 1:length(vectcorrespondanceselect)
    selectcorrespondance(ii) = find(selectcorrespondance2==vectcorrespondanceselect(ii));
end

figure(); hold on;
for i = 1:length(catap)
    plot3(catmidline(i),catap(i),-catdepth(i),'.','Color',colormapselect(selectcorrespondance(i),:))
end
ylabel('Relative AP coordinate [\mum]');xlabel(['ML coordinate [\mum]']);zlabel('Depth [\mum]');
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_Social3dcoord.pdf'],'-dpdf')
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_Social3dcoord.svg'],'-dsvg')

