close all

animal = 'RoofBuddy2';
brainside = 'LIC';

datatable = readtable(['N:\Jenni\' animal '\Zstackinfopersite.xlsx']);


Correspondingsitesind = find(cellfun(@(x) strcmp(x,brainside),datatable.Side));
if strcmp(animal,'RoofBuddy2') % specify the side of the brain too
    Mostanteriorsite = 10;%num2str([datatable.Site(Correspondingsitesind)])
elseif strcmp(animal,'RoofBuddy1') % specify the side of the brain too
    Mostanteriorsite = 9;%num2str([datatable.Site(Correspondingsitesind)])
    if strcmp(brainside,'LIC')% remove site 1 that used different stims
        Correspondingsitesind(find(Correspondingsitesind)==1) = [];
    end
end
Mostanteriorind = find(datatable.Site==Mostanteriorsite);
% Correspondingsitesind(find(datatable.Site(Correspondingsitesind)==10)) = [];
fixed = imread(['N:\Jenni\' animal '\' datatable.VascImage{Mostanteriorind} '.PNG']);%(uint16(fliplr(slice_green{zstackimagecorr})));
fixed = squeeze(fixed(:,:,1));
Widefieldsize = [1210,900]; %x is (690px/0.57 microns) and y is 512/0.57 if microns to pix are the same in widefield (But seems like it from the size estomation of the surgery pic)
Widefieldpix = [size(fixed,2),size(fixed,1)];
foursize = [round(750/2.55),round(512/2.66)]; % X AXIS IS CORRECTED
twosize = [round(750/1.32),round(512/1.37)];
ratiowidefieldfour = round(foursize)./round(Widefieldsize);
ratiowidefieldtwo = round(twosize)./round(Widefieldsize);
pixratiofour =  round(Widefieldpix.*ratiowidefieldfour); % size of fov in pixel for widefield
pixratiotwo = round(Widefieldpix.*ratiowidefieldtwo);
roissites = find(datatable.Rois(Correspondingsitesind)==1);
widefieldratio = fliplr(Widefieldpix)./(Widefieldsize);
topleftcornerrefsite = load(['N:\Jenni\' animal '\z-stack\RICrefcoord.mat'],'RICrefcoord')
col = [1,1,1];
% function megamapplot()
h = figure(); hold on;
imagesc([1:size(fixed,2)]/widefieldratio(1),[1:size(fixed,1)]/widefieldratio(2),flipud(fliplr(fixed))); colormap gray; axis tight


for st = 1:length(roissites)
    clear ytemp xtemp
    moving = imread([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\' datatable.VascImage{Correspondingsitesind(roissites(st))} '.PNG']);
    moving = squeeze(moving(:,:,1));
    
    field(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '_fieldnewcoordpx.mat'],'fieldnewcoordpx')
    fov(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '_fovdnewcoordpx.mat'],'fovnewcoordpx')
    
    if strcmp(brainside,'LIC')
        imagesc([[1:size(fixed,2)]-field(st).fieldnewcoordpx(1)]/widefieldratio(1),[[1:size(fixed,1)]+field(st).fieldnewcoordpx(2)]/widefieldratio(2),flipud(fliplr(moving)),'AlphaData',0.7); colormap gray
    else
        imagesc([[1:size(fixed,2)]-field(st).fieldnewcoordpx(1)]/widefieldratio(1),[[1:size(fixed,1)]-field(st).fieldnewcoordpx(2)]/widefieldratio(2),flipud(fliplr(moving)),'AlphaData',0.7); colormap gray
    end
end


for st = 1:length(roissites)%length(Correspondingsitesind)
    clear xvalues yvalues
    if st == 5 && strcmp(animal,'RoofBuddy2')
        ytemp = (fov(st).fovnewcoordpx(2));% y 2X seems to be offcenter
        xtemp = fov(st).fovnewcoordpx(1);
    else
        if strcmp(brainside,'LIC')
            if  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
                if field(st).fieldnewcoordpx(1)==0
                    xtemp = fov(st).fovnewcoordpx(1);
                    ytemp = fov(st).fovnewcoordpx(2);
                else
                    xtemp = fov(1).fovnewcoordpx(1)-field(st).fieldnewcoordpx(1);%(((Widefieldpix(1)/2)-round(pixratiotwo(1)/2))-field(st).fieldnewcoordpx(1));%-field(st).fieldnewcoordpx(1))/2%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
                    ytemp = fov(1).fovnewcoordpx(2)+field(st).fieldnewcoordpx(2);%((((Widefieldpix(2)/2)-round(pixratiotwo(2)/2)))+field(st).fieldnewcoordpx(2));%-round(pixratiotwo(2)/2)
                end
            elseif  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
                xtemp = ((Widefieldpix(1)/2))+field(st).fieldnewcoordpx(1);%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
                ytemp = ((((Widefieldpix(2)/2)-round(pixratiofour(2)/2)))+field(st).fieldnewcoordpx(2));
            end
        else
            % weird that they are not the same apart from the ratio being
            % different is it because 2X seemed off center ?
            if  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
                xtemp = fov(1).fovnewcoordpx(1)-field(st).fieldnewcoordpx(1);%(((Widefieldpix(1)/2)-round(pixratiotwo(1)/2))-field(st).fieldnewcoordpx(1));%-field(st).fieldnewcoordpx(1))/2%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
                ytemp = fov(1).fovnewcoordpx(2)+field(st).fieldnewcoordpx(2);%((((Widefieldpix(2)/2)-round(pixratiotwo(2)/2)))-field(st).fieldnewcoordpx(2))+70;%-round(pixratiotwo(2)/2)
            elseif  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
                xtemp = ((Widefieldpix(1)/2))-field(st).fieldnewcoordpx(1);%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
                ytemp = ((((Widefieldpix(2)/2)-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(2));
            end
        end
    end
    
    
    %
    %         if st == 5
    %             ytemp = (fov(st).fovnewcoordpx(2));
    %             xtemp = abs(fov(st).fovnewcoordpx(1)-687);
    %         end
    if st == 5 && strcmp(animal,'RoofBuddy2')
        xvalues = [xtemp,xtemp,(xtemp+fov(st).fovnewcoordpx(3)),(xtemp+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
        yvalues = [ytemp,(ytemp+fov(st).fovnewcoordpx(4)),(ytemp+fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2);
    elseif strcmp(brainside,'RIC')
        xvalues = [xtemp,xtemp,(xtemp+fov(st).fovnewcoordpx(3)),(xtemp+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
        yvalues = [ytemp,(ytemp-fov(st).fovnewcoordpx(4)),(ytemp-fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2); % used to be +
    elseif strcmp(brainside,'LIC')
        xvalues = [(xtemp),(xtemp),((xtemp)+fov(st).fovnewcoordpx(3)),((xtemp)+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
        yvalues = [ytemp,(ytemp+fov(st).fovnewcoordpx(4)),(ytemp+fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2);
    end
    plot(xvalues,yvalues,'--','Color',[1,1,1]);
    
    
    %compute coordinate from reference image for all fov
    coordpixsite(datatable.Site(Correspondingsitesind(roissites(st)))).xvaluesmicron(1) = xvalues(1) + topleftcornerrefsite.RICrefcoord(1);
    coordpixsite(datatable.Site(Correspondingsitesind(roissites(st)))).yvaluesmicron(1) = Widefieldsize(2)-yvalues(1) + topleftcornerrefsite.RICrefcoord(2);
    if strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
        coordpixsite(datatable.Site(Correspondingsitesind(roissites(st)))).xfovsize = linspace(1,foursize(1),750);
        coordpixsite(datatable.Site(Correspondingsitesind(roissites(st)))).yfovsize = linspace(1,foursize(2),512);
    elseif strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
        coordpixsite(datatable.Site(Correspondingsitesind(roissites(st)))).xfovsize = linspace(1,twosize(1),750);
        coordpixsite(datatable.Site(Correspondingsitesind(roissites(st)))).yfovsize = linspace(1,twosize(2),512);
    end
end
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_SiteoverVac'],'-dpdf')
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_SiteoverVac.svg'],'-dsvg')


% plot rois over map
% 
% h = figure(); hold on;
% 
% imagesc([1:size(fixed,2)]/widefieldratio(1),[1:size(fixed,1)]/widefieldratio(2),flipud(fliplr(fixed))); colormap gray; axis tight
% 
% for st = 1:length(roissites)%length(Correspondingsitesind)
%     clear ytemp xtemp
%     if strcmp(animal,'RoofBuddy2')
%         moving = imread([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\' datatable.VascImage{Correspondingsitesind(roissites(st))} '.PNG']);
%         moving = squeeze(moving(:,:,1));
%         
%         field(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '_fieldnewcoordpx.mat'],'fieldnewcoordpx')
%         fov(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '_fovdnewcoordpx.mat'],'fovnewcoordpx')
%         
%         %     warp(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) 'warp.mat'],'warp')
%         imagesc([[1:size(fixed,2)]-field(st).fieldnewcoordpx(1)]/widefieldratio(1),[[1:size(fixed,1)]-field(st).fieldnewcoordpx(2)]/widefieldratio(2),flipud(fliplr(moving)),'AlphaData',0.7); colormap gray
%         % end
%         % add field site
%     elseif  strcmp(animal,'RoofBuddy1')
%         field(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '_fieldnewcoordpx.mat'],'fieldnewcoordpx')
%         fov(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '_fovdnewcoordpx.mat'],'fovnewcoordpx')
%         
%         if strcmp(brainside,'LIC')
%             imagesc([[1:size(fixed,2)]-field(st).fieldnewcoordpx(1)]/widefieldratio(1),[[1:size(fixed,1)]+field(st).fieldnewcoordpx(2)]/widefieldratio(2),flipud(fliplr(moving)),'AlphaData',0.7); colormap gray
%         else
%             imagesc([[1:size(fixed,2)]-field(st).fieldnewcoordpx(1)]/widefieldratio(1),[[1:size(fixed,1)]-field(st).fieldnewcoordpx(2)]/widefieldratio(2),flipud(fliplr(moving)),'AlphaData',0.7); colormap gray
%         end
%     end
% end

h = figure(); hold on;
imagesc([1:size(fixed,2)]/widefieldratio(1),[1:size(fixed,1)]/widefieldratio(2),flipud(fliplr(fixed))); colormap gray; axis tight


for st = 1:length(roissites)
    clear ytemp xtemp
    moving = imread([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\' datatable.VascImage{Correspondingsitesind(roissites(st))} '.PNG']);
    moving = squeeze(moving(:,:,1));
    
    field(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '_fieldnewcoordpx.mat'],'fieldnewcoordpx')
    fov(st) = load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '_fovdnewcoordpx.mat'],'fovnewcoordpx')
    
    if strcmp(brainside,'LIC')
        imagesc([[1:size(fixed,2)]-field(st).fieldnewcoordpx(1)]/widefieldratio(1),[[1:size(fixed,1)]+field(st).fieldnewcoordpx(2)]/widefieldratio(2),flipud(fliplr(moving)),'AlphaData',0.7); colormap gray
    else
        imagesc([[1:size(fixed,2)]-field(st).fieldnewcoordpx(1)]/widefieldratio(1),[[1:size(fixed,1)]-field(st).fieldnewcoordpx(2)]/widefieldratio(2),flipud(fliplr(moving)),'AlphaData',0.7); colormap gray
    end
end

colormaptone = jet(10);
for tid = 1:10
    
    ToneColorMap = colormaptone(tid,:);%[linspace(1, 0, 124), (zeros(1, 132))];
    if ToneColorMap(:,1)==1
        logvalues1 = ones(1,257)';
    elseif  ToneColorMap(:,1)==0
        logvalues1 = (log([(ToneColorMap(:,1)+0.1)*10:((10-1)/256):10])/log(10))';
    else
        logvalues1 = (log([linspace(ToneColorMap(:,1)*10,10,257)])/log(10))';
    end
    if ToneColorMap(:,2)==1
        logvalues2 = ones(1,257)';
    elseif  ToneColorMap(:,2)==0
        logvalues2 = (log([(ToneColorMap(:,2)+0.1)*10:((10-1)/256):10])/log(10))';
    else
        logvalues2 = (log([linspace(ToneColorMap(:,2)*10,10,257)])/log(10))';
    end
    if ToneColorMap(:,3)==1
        logvalues3 = ones(1,257)';
    elseif  ToneColorMap(:,3)==0
        logvalues3 = (log([(ToneColorMap(:,3)+0.1)*10:((10-1)/256):10])/log(10))';
    else
        logvalues3 = (log([linspace(ToneColorMap(:,3)*10,10,257)])/log(10))';
    end
    colorMap = flipud([logvalues1,logvalues2,logvalues3]);
    cmap{tid} = colorMap;
end
newtonecolormap = cell2mat(cellfun(@(x) x(end,:),cmap,'UniformOutput',0)');

for st = 1:length(roissites)%length(Correspondingsitesind)
    clear xvalues yvalues
     if st == 5 && strcmp(animal,'RoofBuddy2')
        ytemp = (fov(st).fovnewcoordpx(2));% y 2X seems to be offcenter
        xtemp = fov(st).fovnewcoordpx(1);
    else
        if strcmp(brainside,'LIC')
            if  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
                if field(st).fieldnewcoordpx(1)==0
                    xtemp = fov(st).fovnewcoordpx(1);
                    ytemp = fov(st).fovnewcoordpx(2);
                else
                    xtemp = fov(1).fovnewcoordpx(1)-field(st).fieldnewcoordpx(1);%(((Widefieldpix(1)/2)-round(pixratiotwo(1)/2))-field(st).fieldnewcoordpx(1));%-field(st).fieldnewcoordpx(1))/2%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
                    ytemp = fov(1).fovnewcoordpx(2)+field(st).fieldnewcoordpx(2);%((((Widefieldpix(2)/2)-round(pixratiotwo(2)/2)))+field(st).fieldnewcoordpx(2));%-round(pixratiotwo(2)/2)
                end
            elseif  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
                xtemp = ((Widefieldpix(1)/2))+field(st).fieldnewcoordpx(1);%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
                ytemp = ((((Widefieldpix(2)/2)-round(pixratiofour(2)/2)))+field(st).fieldnewcoordpx(2));
            end
        else
            if  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
                xtemp = fov(1).fovnewcoordpx(1)-field(st).fieldnewcoordpx(1);%(((Widefieldpix(1)/2)-round(pixratiotwo(1)/2))-field(st).fieldnewcoordpx(1));%-field(st).fieldnewcoordpx(1))/2%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
                ytemp = fov(1).fovnewcoordpx(2)+field(st).fieldnewcoordpx(2);%((((Widefieldpix(2)/2)-round(pixratiotwo(2)/2)))-field(st).fieldnewcoordpx(2))+70;%-round(pixratiotwo(2)/2)
            elseif  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
                xtemp = ((Widefieldpix(1)/2))-field(st).fieldnewcoordpx(1);%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
                ytemp = ((((Widefieldpix(2)/2)-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(2));
            end
        end
    end
    
    
    %
    %         if st == 5
    %             ytemp = (fov(st).fovnewcoordpx(2));
    %             xtemp = abs(fov(st).fovnewcoordpx(1)-687);
    %         end
    if st == 5 && strcmp(animal,'RoofBuddy2')
        xvalues = [xtemp,xtemp,(xtemp+fov(st).fovnewcoordpx(3)),(xtemp+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
        yvalues = [ytemp,(ytemp+fov(st).fovnewcoordpx(4)),(ytemp+fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2);
    elseif strcmp(brainside,'RIC')
        xvalues = [xtemp,xtemp,(xtemp+fov(st).fovnewcoordpx(3)),(xtemp+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
        yvalues = [ytemp,(ytemp-fov(st).fovnewcoordpx(4)),(ytemp-fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2); % used to be +
    elseif strcmp(brainside,'LIC')
        xvalues = [(xtemp),(xtemp),((xtemp)+fov(st).fovnewcoordpx(3)),((xtemp)+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
        yvalues = [ytemp,(ytemp+fov(st).fovnewcoordpx(4)),(ytemp+fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2);
    end
    plot(xvalues,yvalues,'--','Color',[1,1,1]);
    
    
%     if strcmp(animal,'RoofBuddy2')
%         if st == 5
%             ytemp = (fov(st).fovnewcoordpx(2))+30;% y 2X seems to be offcenter
%             xtemp = fov(st).fovnewcoordpx(1);
%         else
%             if  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
%                 xtemp = (((Widefieldpix(1)/2)-round(pixratiotwo(1)/2))-field(st).fieldnewcoordpx(1));%-field(st).fieldnewcoordpx(1))/2%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
%                 ytemp = ((((Widefieldpix(2)/2)-round(pixratiotwo(2)/2)))-field(st).fieldnewcoordpx(2))+30;%-round(pixratiotwo(2)/2)
%             elseif  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
%                 xtemp = ((Widefieldpix(1)/2))-field(st).fieldnewcoordpx(1);%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
%                 ytemp = ((((Widefieldpix(2)/2)-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(2))+15;
%             end
%         end
%         if st == 5
%             xvalues = [xtemp,xtemp,(xtemp+fov(st).fovnewcoordpx(3)),(xtemp+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
%             yvalues = [ytemp,(ytemp+fov(st).fovnewcoordpx(4)),(ytemp+fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2);
%         else
%             xvalues = [xtemp,xtemp,(xtemp+fov(st).fovnewcoordpx(3)),(xtemp+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
%             yvalues = [ytemp,(ytemp+fov(st).fovnewcoordpx(4)),(ytemp+fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2);
%         end
%     elseif strcmp(animal,'RoofBuddy1')
%         if strcmp(brainside,'LIC')
%             if  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
%                 if field(st).fieldnewcoordpx(1)==0
%                     xtemp = fov(st).fovnewcoordpx(1);
%                     ytemp = fov(st).fovnewcoordpx(2);
%                 else
%                     xtemp = fov(1).fovnewcoordpx(1)-field(st).fieldnewcoordpx(1);%(((Widefieldpix(1)/2)-round(pixratiotwo(1)/2))-field(st).fieldnewcoordpx(1));%-field(st).fieldnewcoordpx(1))/2%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
%                     ytemp = fov(1).fovnewcoordpx(2)+field(st).fieldnewcoordpx(2);%((((Widefieldpix(2)/2)-round(pixratiotwo(2)/2)))+field(st).fieldnewcoordpx(2));%-round(pixratiotwo(2)/2)
%                 end
%             elseif  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
%                 xtemp = ((Widefieldpix(1)/2))+field(st).fieldnewcoordpx(1);%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
%                 ytemp = ((((Widefieldpix(2)/2)-round(pixratiofour(2)/2)))+field(st).fieldnewcoordpx(2));
%             end
%         else
%             if  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
%                 xtemp = fov(1).fovnewcoordpx(1)-field(st).fieldnewcoordpx(1);%(((Widefieldpix(1)/2)-round(pixratiotwo(1)/2))-field(st).fieldnewcoordpx(1));%-field(st).fieldnewcoordpx(1))/2%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
%                 ytemp = fov(1).fovnewcoordpx(2)+field(st).fieldnewcoordpx(2);%((((Widefieldpix(2)/2)-round(pixratiotwo(2)/2)))-field(st).fieldnewcoordpx(2))+70;%-round(pixratiotwo(2)/2)
%             elseif  strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
%                 xtemp = ((Widefieldpix(1)/2))-field(st).fieldnewcoordpx(1);%(round(Widefieldpix(2)/2))-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(1);
%                 ytemp = ((((Widefieldpix(2)/2)-round(pixratiofour(2)/2)))-field(st).fieldnewcoordpx(2));
%             end
%         end
%         xvalues = [(xtemp),(xtemp),((xtemp)+fov(st).fovnewcoordpx(3)),((xtemp)+fov(st).fovnewcoordpx(3)),(xtemp)]./widefieldratio(1);
%         yvalues = [ytemp,(ytemp+fov(st).fovnewcoordpx(4)),(ytemp+fov(st).fovnewcoordpx(4)),ytemp,ytemp]./widefieldratio(2);
%     end
%     plot(xvalues,yvalues,'--','Color',[1,1,1]);

    % plot Rois onto map
    %load roiscoord
    load([datatable.datadrive{Correspondingsitesind(roissites(st))} ':\Jenni\' animal '\' datatable.subfolder{Correspondingsitesind(roissites(st))} '\site' num2str(datatable.Site(Correspondingsitesind(roissites(st)))) '\31.25\31.25_rois_coord_plane0.mat']);
    
    % %     cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',[((750*widefieldratio(1))-x(:,1)/4)+(xvalues(1)),yvalues(2)+(x(:,2)/4)],'EdgeColor','none','FaceColor',col),((roisCoord{1})))
   if strcmp(animal,'RoofBuddy2') 
    if strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'2X')
        newcoord = cellfun(@(x) [((750-(x(:,1)))/1.32)+xvalues(1),(x(:,2)/1.37)+yvalues(1)],roisCoord{1},'UniformOutput',0);
    else strcmp(datatable.zoomRecording{Correspondingsitesind(roissites(st))},'4X')
        newcoord = cellfun(@(x) [((750-(x(:,1)))/2.55)+xvalues(1),(x(:,2)/2.66)+yvalues(1)],roisCoord{1},'UniformOutput',0);
    end
   elseif strcmp(animal,'RoofBuddy1') 
               newcoord = cellfun(@(x) [((890-(x(:,1)))/1.32)+xvalues(1),(x(:,2)/1.37)+yvalues(1)],roisCoord{1},'UniformOutput',0);
   end
    if st == 5 && strcmp(animal,'RoofBuddy2')
        newcoord = cellfun(@(x) [((750-(x(:,1)))/1.32)+xvalues(1),(x(:,2)/1.37)+yvalues(1)],roisCoord{1},'UniformOutput',0);
    end
    % cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',[(((750)*widefieldratio(1))-x(:,1)/4),((x(:,2)/4)+yvalues(2))],'EdgeColor','none','FaceColor',col),((roisCoord{1})))
    cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',col),newcoord)
    load([datatable.datadrive{Correspondingsitesind(roissites(1))} ':\Jenni\' animal '\z-stack\Tonotopystruct.mat']);
    
    tunedcells = find(Tonotopystruct(datatable.Site(Correspondingsitesind(roissites(st)))).cellstuning(2,:));
    bftunedcellsind{st} =  Tonotopystruct(datatable.Site(Correspondingsitesind(roissites(st)))).cellstuning(3,tunedcells);
    
    cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',newtonecolormap(y,:)),newcoord(tunedcells),mat2cell(bftunedcellsind{st},1,ones(1,length(bftunedcellsind{st}))))
    if strcmp(brainside,'RIC')
        load('N:\Jenni\RoofBuddy2\z-stack\RICrefcoord.mat')
        Relativecoord = cellfun(@(x) (x)+RICrefcoord,newcoord,'UniformOutput',0);
    else
        load('N:\Jenni\RoofBuddy2\z-stack\LICrefcoord.mat')
        Relativecoord = cellfun(@(x) [x(:,1)+LICrefcoord(1),x(:,2)-LICrefcoord(2)],newcoord,'UniformOutput',0);
    end
    midlinecoord{st} = cellfun(@(x) median(x(:,1)),Relativecoord(tunedcells));
    anteropostlinecoord{st} = cellfun(@(x) median(x(:,2)),Relativecoord(tunedcells));
    depthtuning{st} = Tonotopystruct(datatable.Site(Correspondingsitesind(roissites(st)))).depth(tunedcells);
end
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_SiteoverVacWithRois.pdf'],'-dpdf')
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_SiteoverVacWithRois.svg'],'-dsvg')

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
subplot(1,3,3); plot(cattuning,catdepth,'.'); title('Depth')
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_Distancetuningsummary'],'-dpdf')

figure(); hold on; subplot(1,3,1); plot(catmidline,cattuning,'.'); title('Midline distance')
subplot(1,3,2); plot(catap,cattuning,'.'); title('AP distance')
subplot(1,3,3); plot(catdepth,cattuning,'.'); title('Depth')
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_Distancetuningsummary2'],'-dpdf')
figure(); hold on;
for i = 1:length(catap)
    plot3(catmidline(i),-catap(i),-catdepth(i),'.','Color',newtonecolormap(cattuning(i),:))
end
ylabel('Relative AP coordinate [\mum]');xlabel(['ML coordinate [\mum]', ' n_c_e_l_l_s' num2str(length(cattuning)), 'n_s_i_t_e_s' num2str(st)]);zlabel('Depth [\mum]');
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_Tonotopy3dcoord.pdf'],'-dpdf')
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_Tonotopy3dcoord.svg'],'-dsvg')


figure(); hold on;
scatter3(catmidline,-catap,-catdepth,5,newtonecolormap(cattuning,:),'filled');
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































load('N:\Jenni\RoofBuddy2\roistructure\socialselectivity.mat')
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
    plot3(catmidline(i),-catap(i),-catdepth(i),'.','Color',colormapselect(selectcorrespondance(i),:))
end
ylabel('Relative AP coordinate [\mum]');xlabel(['ML coordinate [\mum]']);zlabel('Depth [\mum]');
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_Social3dcoord.pdf'],'-dpdf')
print(['N:\Jenni\' animal '\figure\tonotopy\' animal brainside '_Social3dcoord.svg'],'-dsvg')

