function aligmenttoref(type)
% do ratio and alignment
% animal = 'RoofBuddy2';
% brainside = 'RIC';
animal = 'Gray2';
brainside = 'RIC';

datatable = readtable(['K:\Jenni\' animal '\Zstackinfopersite.xlsx']);


Correspondingsitesind = intersect(find(cellfun(@(x) strcmp(x,brainside),datatable.Side)),find(datatable.Rois==1));
% make sure those correspond to the alignement to surgery picture see
% Alignvasc.m
if strcmp(animal,'RoofBuddy2') 
    if strcmp(brainside,'LIC')
        Mostanteriorsite = 21;%num2str([datatable.Site(Correspondingsitesind)])
    elseif strcmp(brainside,'RIC')
        Mostanteriorsite = 10;%num2str([datatable.Site(Correspondingsitesind)])
    end
elseif strcmp(animal,'RoofBuddy1')
    Mostanteriorsite = 9;%num2str([datatable.Site(Correspondingsitesind)])
    if strcmp(brainside,'LIC')% remove site 1 that used different stims
        Correspondingsitesind(find(Correspondingsitesind)==1) = [];
    end
    elseif strcmp(animal,'Gray2')
    Mostanteriorsite = 1;%num2str([datatable.Site(Correspondingsitesind)])
end
Mostanteriorind = find(datatable.Site==Mostanteriorsite);
Correspondingsitesind(find(datatable.Site(Correspondingsitesind)==Mostanteriorsite)) = [];
fixed = imread([datatable.datadrive{Mostanteriorind} ':\Jenni\' animal '\' datatable.VascImage{Mostanteriorind} '.PNG']);%(uint16(fliplr(slice_green{zstackimagecorr})));
fixed = squeeze(fixed(:,:,1));
if strcmp(brainside,'LIC')
    fixed = fliplr(fixed);
end
Widefieldsize = [1210,900]; %x is (690px/0.57 microns) and y is 512/0.57 if microns to pix are the same in widefield (But seems like it from the size estomation of the surgery pic)
Widefieldpix = [size(fixed,1),size(fixed,2)];
foursize = [round(750/2.55),round(512/2.66)]; % X AXIS IS CORRECTED
twosize = [round(750/1.32),round(512/1.37)];
ratiowidefieldfour = round(foursize)./round(Widefieldsize);
ratiowidefieldtwo = round(twosize)./round(Widefieldsize);
pixratiofour =  round(Widefieldpix.*ratiowidefieldfour);
pixratiotwo = round(Widefieldpix.*ratiowidefieldtwo);
roissites = find(datatable.Rois(Correspondingsitesind)==1);
widefieldratio = fliplr(Widefieldpix)./(Widefieldsize);

% do it for the ref site
fieldnewcoordpx = [0,0];
y1 = round((Widefieldpix(1)/2)-round(pixratiotwo(1)/2));
x1 = round((Widefieldpix(2)/2)-round(pixratiotwo(2)/2));
fovnewcoordpx = [x1,y1, pixratiotwo(2),pixratiotwo(1)];
warp = fixed;

save([datatable.datadrive{(Mostanteriorsite)} ':\Jenni\' animal '\z-stack\' 'site' num2str(Mostanteriorsite) '_fieldnewcoordpx.mat'],'fieldnewcoordpx')
save([datatable.datadrive{(Mostanteriorsite)} ':\Jenni\' animal '\z-stack\' 'site' num2str(Mostanteriorsite) '_fovdnewcoordpx.mat'],'fovnewcoordpx')
save([datatable.datadrive{(Mostanteriorsite)} ':\Jenni\' animal '\z-stack\' 'site' num2str(Mostanteriorsite) 'warp.mat'],'warp')


for st = 1:length(Correspondingsitesind)%3:max([datatable.Site]') for site 21 %[1:3,5]% intersect that with Rois drawn
    clear warp vecstartpy vecstartpx fieldnewcoordpx fovdnewcoordpx cropmoving
    recordingX = str2num(datatable.zoomRecording{Correspondingsitesind(st)}(1));
    switch type
        case 'automatic'
            if strcmp(animal,'RoofBuddy1')
                figure(); hold on;  subplot(1,3,1);imagesc(fixed); colormap gray
                moving = imread([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\' datatable.VascImage{Correspondingsitesind(st)} '.PNG']);
                moving = fliplr(squeeze(moving(:,:,1)));
                subplot(1,3,2);imagesc(moving); colormap gray
                refx = datatable.Xrefpoint(Mostanteriorsite);
                refy = datatable.Yrefpoint(Mostanteriorsite);
                diffx = datatable.Xrefpoint(Correspondingsitesind(st))-datatable.Xrefpoint(Mostanteriorsite);
                diffy = datatable.Yrefpoint(Correspondingsitesind(st))-datatable.Yrefpoint(Mostanteriorsite);
                cropmoving = moving(1:end+diffy,1:end-abs(diffx));
                warp = zeros(size(fixed,1),size(fixed,2));
                if size(cropmoving,1)==240
                    warp(abs(diffy):end,(abs(diffx)):end) = cropmoving;
                elseif size(cropmoving,1)==481
                    warp(abs(diffy)+2:end,(abs(diffx)-1):end) = cropmoving;
                elseif size(cropmoving,1)==341
                    warp(abs(diffy)+1:end,(abs(diffx)):end) = cropmoving;
                else
                    warp(abs(diffy)+1:end,(abs(diffx)+1):end) = cropmoving;
                end
                subplot(1,3,3); imagesc(warp)
                figure(); hold on; imshowpair(uint16(fixed),uint16(warp))
                fieldnewcoordpx(1) = diffx; % in widefield pixel from the ref point either minus or plus depending if closest to the midline or not (minus if closes to the midline)
                fieldnewcoordpx(2) = diffy; % always minus because the ref site is the most anterior
                if recordingX == 4
                    x1warp = diffx;
                    y1warp = diffy;
                    fovnewcoordpx = round([x1warp, y1warp, pixratiotwo(2), pixratiotwo(1)]);
                    
                elseif recordingX ==2
                    x1warp = diffx;
                    y1warp = diffy;
                    fovnewcoordpx = round([x1warp, y1warp, pixratiotwo(2), pixratiotwo(1)]);
                end
            else
                % make an exception for sites that have the same ref image (RFB2
                % site 21 and 22 for example or 14,15,16) % LIC not aligning
                % properly to the fixed image
                moving = imread([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\' datatable.VascImage{Correspondingsitesind(st)} '.PNG']);
                moving = squeeze(moving(:,:,1));
                if strcmp(brainside,'LIC')
                    moving = fliplr(moving);
                end
                figure
                imshowpair(imadjust(fixed), imadjust(moving))
                [optimizer, metric] = imregconfig('multimodal')
                optimizer.InitialRadius = 6.25e-10;
                optimizer.Epsilon = 1.5e-6;
                optimizer.GrowthFactor = 1.01;
                optimizer.MaximumIterations = 500;
                tform = imregtform(imadjust(moving),imadjust(fixed),'translation',optimizer,metric);
                warp=imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
                figure; hold on
                imshowpair(fixed,warp)
                tx=tform.T(3,1);
                ty=tform.T(3,2);
                
                [fieldnewcoordpx(1),fieldnewcoordpx(2)] = transformPointsForward(tform,1,1);
                
                %         for ii = 1:size(warp,1)
                %             if isempty(find(warp(ii,:)>0))
                %                 vecstartpx(ii) = size(warp,2);
                %             else
                %                 vecstartpx(ii) = min(find(warp(ii,:)>0));
                %             end
                %         end
                %         for ii = 1:size(warp,2)
                %             if isempty(find(warp(:,ii)>0))
                %                 vecstartpy(ii) = size(warp,1);
                %             else
                %                 vecstartpy(ii) = max(find(warp(:,ii)==0));
                %             end
                %         end
                
                if recordingX == 4
                    vecstartpxtmp = vecstartpx;
                    vecstartpxtmp(vecstartpxtmp==Widefieldpix(2)) = 0;
                    vecstartpytmp = (vecstartpy);
                    vecstartpytmp(vecstartpytmp==Widefieldpix(1)) = 0;
                    
                    %              fielddisplacementx = vecstartpxtmp(1,(Widefieldpix(1)/2)-round(pixratiofour(1)/2):((Widefieldpix(1)/2)-round(pixratiofour(1)/2))+pixratiofour(1));
                    %              fielddisplacementy = vecstartpytmp(1,(Widefieldpix(2)/2)-round(pixratiofour(2)/2):((Widefieldpix(2)/2)-round(pixratiofour(2)/2))+pixratiofour(2));
                    
                    
                    figure(); hold on;
                    imagesc(flipud(moving)); colormap gray
                    x1 = round((Widefieldpix(1)/2)-round(pixratiofour(1)/2));
                    x2 = round(x1+pixratiofour(1));
                    y1 = round((Widefieldpix(2)/2)-round(pixratiofour(2)/2));
                    y2 = round(y1+pixratiofour(2));
                    
                    x = [x1, x2, x2, x1, x1];
                    y = [y1, y1, y2, y2, y1];
                    plot(y, x, 'k--');
                    axis tight
                    
                    figure(); hold on;
                    imagesc(flipud(warp)); colormap gray; axis tight
                    
                    [x1warp,y1warp] = transformPointsForward(tform,y1,(x1));%x1/(size(warp,1)/ty); first argument is x second is y
                    %              [x2warp,y2warp] = transformPointsForward(tform,y2,(x2));%x1/(size(warp,1)/ty);
                    
                    
                    x = [x1warp, x1warp, x1warp+pixratiofour(2), x1warp+pixratiofour(2),x1warp];
                    y = size(warp,1)- [y1warp, y1warp+pixratiofour(1), y1warp+pixratiofour(1), y1warp, y1warp];
                    plot(x, y,'k--');
                    fovnewcoordpx = [x1warp,size(warp,1)-y1warp, pixratiofour(2), pixratiofour(1)];
                    
                elseif recordingX ==2
                    %              vecstartpxtmp = vecstartpx;
                    %              vecstartpxtmp(vecstartpxtmp==Widefieldpix(2)) = 0;
                    %              vecstartpytmp = vecstartpy;
                    %              vecstartpytmp(vecstartpytmp==Widefieldpix(1)) = 0;
                    %
                    %              fielddisplacementx = vecstartpxtmp(1,(Widefieldpix(1)/2)-round(pixratiotwo(1)/2):((Widefieldpix(1)/2)-round(pixratiotwo(1)/2))+pixratiotwo(1));
                    %              fielddisplacementy = vecstartpytmp(1,(Widefieldpix(2)/2)-round(pixratiotwo(2)/2):((Widefieldpix(2)/2)-round(pixratiotwo(2)/2))+pixratiotwo(2));
                    
                    
                    figure(); hold on;
                    imagesc(flipud(moving)); colormap gray
                    x1 = round((Widefieldpix(1)/2)-round(pixratiotwo(1)/2));
                    x2 = round(x1+pixratiotwo(1));
                    y1 = round((Widefieldpix(2)/2)-round(pixratiotwo(2)/2));
                    y2 = round(y1+pixratiotwo(2));
                    
                    x = [x1, x2, x2, x1, x1];
                    y = [y1, y1, y2, y2, y1];
                    plot(y, x, 'k--');
                    axis tight
                    
                    figure(); hold on;
                    imagesc(flipud(warp)); colormap gray; axis tight
                    
                    [x1warp,y1warp] = transformPointsForward(tform,y1,(x1));%x1/(size(warp,1)/ty); first argument is x second is y
                    %              [x2warp,y2warp] = transformPointsForward(tform,y2,(x2));%x1/(size(warp,1)/ty);
                    
                    x = [x1warp, x1warp, x1warp+pixratiotwo(2), x1warp+pixratiotwo(2),x1warp];
                    y = size(warp,1)- [y1warp, y1warp+pixratiotwo(1), y1warp+pixratiotwo(1), y1warp, y1warp];
                    plot(x, y,'k--');
                    fovnewcoordpx = [x1warp,size(warp,1)-y1warp, pixratiotwo(2), pixratiotwo(1)];
                end
            end
        case 'manual'
        % x displacement towards midline from ref image in pixel
        fieldnewcoordpx(1) = datatable.Manualx(Correspondingsitesind(st)) % positive towards midline negative to lateral % for LIC
        fieldnewcoordpx(2) = datatable.Manualy(Correspondingsitesind(st)) % positive towards posterior
        if recordingX == 4
            x1 = round((Widefieldpix(1)/2)-round(pixratiofour(1)/2));
            x2 = round(x1+pixratiofour(1));
            y1 = round((Widefieldpix(2)/2)-round(pixratiofour(2)/2));
            y2 = round(y1+pixratiofour(2));
            fovnewcoordpx = [y1+fieldnewcoordpx(2),x1+fieldnewcoordpx(1), pixratiofour(2), pixratiofour(1)];

        elseif recordingX == 2
            x1 = round((Widefieldpix(1)/2)-round(pixratiotwo(1)/2));
            x2 = round(x1+pixratiotwo(1));
            y1 = round((Widefieldpix(2)/2)-round(pixratiotwo(2)/2));
            y2 = round(y1+pixratiotwo(2));
            fovnewcoordpx = [y1+fieldnewcoordpx(2),x1+fieldnewcoordpx(1), pixratiotwo(2), pixratiotwo(1)];
        end
        
    end
    
    
    
    %         startpixelx(Correspondingsitesind(st)) = find(sum(warp,2),1);
    %         startpixely(Correspondingsitesind(st)) = find(sum(warp,1),1);
    %         movedimage(Correspondingsitesind(st),:,:) = warp;
    %(find(warp>0,1));
    %         nbpixcorr(Correspondingsitesind(st)) = [size(moving,1),size(moving,2)]-1;
    %         megamapx = size(warp,1)+max(vecstartpx(vecstartpx<size(warp,2)));%warp;
    %         megamapy = size(warp,2)+max(vecstartpy(vecstartpy<size(warp,1)));%warp;
    %         newmap = zeros(megamapx,megamapy);
    % add the recording site FOV in widefield image
    
    save([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(st))) '_fieldnewcoordpx.mat'],'fieldnewcoordpx')
    save([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(st))) '_fovdnewcoordpx.mat'],'fovnewcoordpx')
%     save([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(st))) 'warp.mat'],'warp')
    %         save([datatable.datadrive{Correspondingsitesind(st)} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(Correspondingsitesind(st))) 'tform.mat'],'tform')
end
% bla = vecstartpy(Correspondingsitesind(st),:);
% bla(find(vecstartpy(Correspondingsitesind(st),:)==size(fixed,1)))=nan;
% bli =  vecstartpx(Correspondingsitesind(st),:);
% bli(find(vecstartpx(Correspondingsitesind(st),:)==size(fixed,2)))=nan;
% figure(); hold on;
% imagesc(1:687,1:514,flipud(fixed)); colormap gray;
% imagesc((1:687)+bla,(1:514)+bli,flipud(warp));
%
end