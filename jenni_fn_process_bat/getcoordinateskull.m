function getcoordinateskull()
animal = 'RoofBuddy2';
datatable = readtable(['N:\Jenni\' animal '\Zstackinfopersite.xlsx']);

for st = 23%3:max([datatable.Site]') for site 21
    
    
    if strcmp(datatable.CorrespondingZstack{st},'NA')
    else
        cd([datatable.datadrive{st} ':\Jenni\' animal '\z-stack'])%load('RoofBuddy2_013_011_green.mat')
        load([datatable.CorrespondingZstack{st} '_green.mat'])
        meanimage = double(read(Tiff([datatable.datadrive{st} ':\Jenni\' animal '\' datatable.subfolder{st} '\site' num2str([datatable.Site(st)]) '\MeanImg31.25.tiff'])));
        recordingX = str2num(datatable.zoomRecording{st}(1));
        zstackX = str2num(datatable.zoomZstack{st}(1));

        meanmeanimage = imresize(meanimage, zstackX/recordingX, 'bilinear');
        
        zstackimagecorr = (datatable.SliceSiteNb(st));%29; % for site 21, 19 for site 22 (about 50 microns under)
        slice_green_mat = cat(3,slice_green{:});
        micronstep = datatable.step(st);
        endofskullslice = (datatable.SliceMaxDepth(st));
        
        %get image correspondance 1X to 2X
        figure(); hold on
        imagesc(fliplr(imadjust(uint16(slice_green{zstackimagecorr})))); axis tight; colormap gray
        figure(); hold on
        imagesc(fliplr(imadjust(uint16((meanmeanimage))))); axis tight; colormap gray
        if strcmp(animal,'RoofBuddy2') && st ==8
            meanimage9 = double(read(Tiff([datatable.datadrive{st} ':\Jenni\' animal '\' datatable.subfolder{st} '\site9\MeanImg31.25.tiff'])));
            meanmeanimage9 = imresize(meanimage, zstackX/recordingX, 'bilinear');
            moving = (uint16(fliplr((meanmeanimage9))));
            fixed = (uint16(fliplr(slice_green{(datatable.SliceSiteNb(st+1))})));
        else
            if strcmp(datatable.Side{st},'RIC')
            fixed = (uint16(fliplr((slice_green{zstackimagecorr}))));
            moving = (uint16(fliplr(((meanmeanimage)))));
            else
            fixed = (uint16(fliplr(flipud(slice_green{zstackimagecorr}))));
            moving = (uint16(fliplr(flipud((meanmeanimage)))));
            end
        end
        figure
        imshowpair(imadjust(fixed), moving)
        [optimizer, metric] = imregconfig('multimodal')
        optimizer.InitialRadius = 6.25e-3;
        optimizer.Epsilon = 1.5e-10;
        optimizer.GrowthFactor = 1.01;
        optimizer.MaximumIterations = 300;
        tform = imregtform(moving,imadjust(fixed),'affine',optimizer,metric);
        warp=imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
        figure
        imshowpair(fixed,warp)
        startpixelx = find(sum(warp,2),1); % site 15 1
        startpixely= find(sum(warp,1),1);%(find(warp>0,1)); site 15 60;%
        nbpixcorr = [size(moving,1),size(moving,2)]-1;
        % medianvalue = median(reshape(slice_green_mat,1,size(slice_green_mat,1)*size(slice_green_mat,2)*size(slice_green_mat,3)));
        % stdvalue = mean(reshape(slice_green_mat,1,size(slice_green_mat,1)*size(slice_green_mat,2)*size(slice_green_mat,3)));
        % maxvalue =  prctile(reshape(slice_green_mat,1,size(slice_green_mat,1)*size(slice_green_mat,2)*size(slice_green_mat,3)),[95]);
        % minvalue =  min(reshape(slice_green_mat,1,size(slice_green_mat,1)*size(slice_green_mat,2)*size(slice_green_mat,3)));
        maxvaluedepth = max(squeeze(max(permute(slice_green_mat(startpixelx:startpixelx+nbpixcorr(1),startpixely:startpixely+nbpixcorr(2),zstackimagecorr:end),[3,1,2]))));
        minvaluedepth = min(squeeze(min(permute(slice_green_mat(startpixelx:startpixelx+nbpixcorr(1),startpixely:startpixely+nbpixcorr(2),1:zstackimagecorr),[3,1,2]))));
        % % meanvaluedepth = squeeze(mean(permute(slice_green_mat,[3,1,2])));
        % % stdvaluedepth = squeeze(std(permute(slice_green_mat,[3,1,2])));
        %
        norm_slice_green = cellfun(@(x) (x(startpixelx:startpixelx+nbpixcorr(1),startpixely:startpixely+nbpixcorr(2),:)-minvaluedepth)./(maxvaluedepth-minvaluedepth),slice_green,'UniformOutput',0);
        rawslice_green = cellfun(@(x) (x(startpixelx:startpixelx+nbpixcorr(1),startpixely:startpixely+nbpixcorr(2),:)),slice_green,'UniformOutput',0);
        
        norm_slice_green_mat = cat(3,norm_slice_green{:});
        rawslice_green_mat = cat(3,rawslice_green{:});
        % plot the corresponding image from the zstack, compare this
        figure();
        hold on; subplot(1,2,1);imagesc(fliplr(flipud(imadjust(uint16((meanmeanimage)))))); axis tight; colormap gray
        hold on; subplot(1,2,2); imagesc(imadjust(uint16(fliplr(flipud(slice_green_mat(startpixelx:startpixelx+nbpixcorr(1),startpixely:startpixely+nbpixcorr(2),zstackimagecorr))))));
        colormap gray
        
        % bla = (squeeze(norm_slice_green_mat(500,500,:))');
        %
        %
        % [b,S,mu]  = polyfit(15:length(bla), bla(15:length(bla)), 3);
        % fy = polyval(b,15:length(bla),S,mu);
        % y = fy;
        %
        % % plot diff per pixel as a function of depth, rapid increase innorm
        % % fluorescence indeicates skull :)
        % imagesc(diff(squeeze(norm_slice_green_mat(500,1:720,:)))); caxis([0,1])
        %
        % skullpoint = cellfun(@(x) find(x>0.5),norm_slice_green,'UniformOutput',0);
        % skullpoint_corrected =  cellfun(@(x) x(find(diff(x)<2)),skullpoint,'UniformOutput',0);
        % test = zeros(size(norm_slice_green{1},1),size(norm_slice_green{1},2));
        % test_corrected = zeros(size(norm_slice_green{1},1),size(norm_slice_green{1},2));
        % test(skullpoint{33}) = 1;
        % test_corrected(skullpoint_corrected{33}) = 1;
        %
        % figure(); hold on;
        % subplot(1,2,1); hold on; imagesc(test); axis tight
        % subplot(1,2,2); hold on; imagesc(test_corrected); axis tight
        distrinorm = reshape(norm_slice_green_mat,1,size(norm_slice_green_mat,1)*size(norm_slice_green_mat,2)*size(norm_slice_green_mat,3));
        distriraw = reshape(rawslice_green_mat,1,size(rawslice_green_mat,1)*size(rawslice_green_mat,2)*size(rawslice_green_mat,3));
        limdistriraw = prctile(distriraw,95);
        limdistrinorm = prctile(distrinorm,95);
        topimage = imagesc(flipud(fliplr(rawslice_green_mat(:,:,end)))); colormap gray
        
        % try to color in blue all the pixel that have a value above the threshold
        % plot distribtuion for 10 top planes see emerging differences
        Thrseholdtopplanenorm = prctile(reshape(norm_slice_green_mat(:,:,end),1,size(norm_slice_green_mat(:,:,end),1)*size(norm_slice_green_mat(:,:,end),2)),50)
        Thrseholdtopplane = prctile(reshape(rawslice_green_mat(:,:,end),1,size(norm_slice_green_mat(:,:,end),1)*size(norm_slice_green_mat(:,:,end),2)),75)
        thresholdallraw = prctile(distriraw,90); % 95 for RoofBuddy2
        figure(); hold on;
        for tp = 1:10
            [a,b]= hist(reshape(norm_slice_green_mat(:,:,end-tp),1,size(norm_slice_green_mat(:,:,end-tp),1)*size(norm_slice_green_mat(:,:,end-tp),2)),50);
            plot(b,a)
        end
        % plot image with skull in nan
        for tp = 1:size(slice_green_mat,3)-1
            clear skullvaltmp
            norm_slice_green_mattmp{tp} = rawslice_green_mat(:,:,end-tp);
            skullvaltmp = find(norm_slice_green_mattmp{tp}>=thresholdallraw);
            %      [skullvaltmpy,skullvaltmpx] = find(norm_slice_green_mattmp>=thresholdallraw);
            norm_slice_green_mattmp{tp}(skullvaltmp) = nan;
            %      norm_slice_green_mattmp() = nan;
%             figure(); hold on; subplot(1,2,1)
%             imagesc(flipud(fliplr(norm_slice_green_mattmp{tp}))); axis tight; colormap gray
%             hold on; subplot(1,2,2); imagesc(imadjust(uint16(fliplr(flipud(slice_green_mat(startpixelx:startpixelx+nbpixcorr(1),startpixely:startpixely+nbpixcorr(2),end-tp))))));
        end
        bla = cellfun(@(x) isnan(x),norm_slice_green_mattmp,'UniformOutput',0);
        blamat = cat(3,bla{:});
        sumbla = sum(blamat,3);
        h = figure(); hold on; ax(1)=subplot(1,4,1);imagesc(flipud(fliplr(sumbla))); colorbar; title('skull points');
        sumbla(find(blamat(:,:,1)==0)) = 0; % removed for RoofBuddy 1 but
%         is that reasonable ?
        ax(2)=subplot(1,4,2);imagesc(flipud(fliplr(sumbla))); colorbar; title('skull points corrected');
        [newdepthmatx,newdepthmaty] = find(sumbla>0);
        % find the depth for ones that are not zeros
        blx = mat2cell([newdepthmatx],ones(1,size(newdepthmatx,1)),1);
        bly = mat2cell([newdepthmaty],ones(1,size(newdepthmatx,1)),1);
        skullthicknessestimate = cellfun(@(x,y) max(find(blamat(x,y,:)==1)*micronstep),blx,bly);
        skull = zeros(size(norm_slice_green_mat(:,:,end-tp),1),size(norm_slice_green_mat(:,:,end-tp),2));
        skull(find(sumbla>0)) = skullthicknessestimate;
        skull(find(skull>(size(norm_slice_green_mat,3)-endofskullslice)*micronstep)) = (size(norm_slice_green_mat,3)-endofskullslice)*micronstep;
        ax(3)=subplot(1,4,3);imagesc(flipud(fliplr(skull))); colorbar; title('skull thickness');
        depth = ones(size(norm_slice_green_mat(:,:,end-tp),1),size(norm_slice_green_mat(:,:,end-tp),2))*(size(slice_green_mat,3)-zstackimagecorr)*micronstep;
        depth = depth-skull; % correct for negative value
        depth(find(depth<0)) = 0; % correct for negative value and replace by max skull depth (estimated by me from z-stack)
        ax(4)=subplot(1,4,4);imagesc(flipud(fliplr(depth))); colorbar;  title('estimated depth');
        colormap(ax(1),gray)
        colormap(ax(2),gray)
        colormap(ax(3),jet)
        colormap(ax(4),jet)
        save([datatable.datadrive{st} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(st)) '_depth.mat'],'depth')
        print([datatable.datadrive{st} ':\Jenni\' animal '\z-stack\' 'site' num2str(datatable.Site(st)) '_DepthEstimate'],'-dpdf')

        % coordpixx = repmat(1:size(rawslice_green_mat,1),1,size(rawslice_green_mat,2));
        % coordpixy = repmat(1:size(rawslice_green_mat,2),1,size(rawslice_green_mat,1));
        
        
        % try threshold approach:
        % % bli = cellfun(@(x) find(x>prctile(distri,95)),norm_slice_green,'UniformOutput',0);
        % for xx = 1:size(norm_slice_green_mat,1)
        %     for yy = 1:size(norm_slice_green_mat,2)
        %         for dd = 1:size(norm_slice_green_mat,3)
        %             if norm_slice_green_mat(xx,yy,dd)<prctile(distrinorm,90)
        %                 skullpointmat(xx,yy,dd) = 0;
        %             else
        %                 skullpointmat(xx,yy,dd) = 1;
        %             end
        %         end
        %         if ~isempty(find(squeeze(skullpointmat(xx,yy,:))==1,1))
        %             firstskullpoint(xx,yy) = length(norm_slice_green);
        %         else
        %             firstskullpoint(xx,yy) = find(squeeze(skullpointmat(xx,yy,:))==1,1);
        %         end
        %         depthvalue(xx,yy) = (firstskullpoint(xx,yy) -(zstackimagecorr))*micronstep;
        %         %         if  depthvalue(xx,yy) <0
        % %            depthvalue(xx,yy) = 50;
        % %         end
        %     end
        % end
        % % newdepthvalue = imresize(depthvalue,2);
        % save(['depthvalue_RoofBuddy2_013_011.mat'],'depthvalue')
    end
end
end