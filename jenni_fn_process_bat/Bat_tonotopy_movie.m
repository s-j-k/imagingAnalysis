function Bat_tonotopy_movie()
animal = 'RoofBuddy2';
datatable = readtable(['N:\Jenni\' animal '\ProcessingProgress.xlsx']);

manualroisdone = find(table2array(datatable(:,13))==1);%find(~cellfun(@isempty,cellfun(@(x) find(x=='1'),table2array(datatable(:,13)),'UniformOutput',0)));%find(table2array(datatable(:,13))==1);
daytraining = table2array(datatable(:,1));
preprocessedrois = unique(daytraining(manualroisdone));
site = table2array(datatable(:,5));
framerate = table2array(datatable(:,8));
protocol = table2array(datatable(:,7));
path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,10));
uniquedaytraining = unique(daytraining);
uniquesite = unique(site);
for ii = 1:length(uniquesite)
    analysefilelistind(ii) = find(site==uniquesite(ii),1);
end

for fn = 8%1:length(analysefilelistind)
    close all
    
    clear TC roiMatevoked tunedcells signicelltuning meantuningdf roiMatevokedstim
    % read frames from h5
    h5list = dir(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' '*.h5']);
    nbframestmp(1) = 0;
    for filenb = 1:length(h5list)
        clear hinfo
        hinfo = hdf5info(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' h5list(filenb).name]);
        nbframestmp(filenb+1) = hinfo.GroupHierarchy.Datasets.Dims(3);
    end
    
    % load behavior file
    load(['N:\Jenni\' animal '\behavior\' behaviorfile{analysefilelistind(fn)} '.mat']); %RoofBuddy1tonotopy02-03-2021-12-50
    startstim = batimagingdata.frametone;%1:400:44000;
    
    seqid = batimagingdata.durationseq;%batimagingdata.;%repmat(1:11,1,10);
    toneid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
    uniquetone = unique(toneid);
    uniqueduration = unique(seqid); % also order of presentation
    % nb of repetition
    nbrep = (length(seqid)/length(uniquetone))/length(uniqueduration);
    xvalues = [-10:40]/31.25;
    for tid = 1:length(uniquetone)
        for dd = 1:length(uniqueduration)
            ordertone(tid,dd,:) = intersect(find(toneid==uniquetone(tid)),find(seqid==uniqueduration(dd))); % tone id duration and trial rep
        end
    end
    % load ROIs
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_TC_plane0.mat']);
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\Fall.mat'])
    load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\TonotopyMeanEvokedMovie.mat'])
    minval = 1;% prctile(reshape(ratio,1,size(ratio,1)*size(ratio,2)),[5]);
    maxval = 1.3;%prctile(reshape(ratio,1,size(ratio,1)*size(ratio,2)),[95]); across all trials
    psth=MeanEvokedMovie; %average across all trials
    %     intersity2 = mean(MeanEvokedMovie{1},3);
    xvalues2 = [-10:44]*31.25;
    if exist((['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2pmeanimg.tif']))
        intensity=imread(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2pmeanimg.tif']); %load your intensity/mean image
    else
        load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\Fall.mat'],'ops')
        intensity=uint16(ops.meanImg); %load your intensity/mean image
        %         meanimg = imadjust(uint16(ops.meanImg));
        %         h=figure(2); hold on; imshow(meanimg); axis tight;
        %         saveas(h,['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2pmeanimg.tif'])
    end
    colormaptone = jet(length(psth));
    
    % create a tonotopic map by replacing pixels with a small dot that has a bf
    % create the tone colormap
    for tid = 1:length(psth)
        clear ToneColorMap
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
        TonescaledcolorMap(tid,:,:) = flipud([logvalues1,logvalues2,logvalues3]);
    end
    
    % define bf per pixel and get the amplitude of the response for that
    % pixel: nb pixel, x*y
    xpixelcoordinates = mat2cell(1:size(psth{1},1),1,ones(1,size(psth{1},1)));
    ypixelcoordinates = mat2cell(1:size(psth{1},2),1,ones(1,size(psth{1},2)));
    matpasth = cat(4,psth{:});
    matpasth = matpasth./median(matpasth(:,:,1:10,:),3);
    
    
    matmaxpeakamp = squeeze(max(matpasth(:,:,11:15,:),[],3));
    [bfpixelamp,bfpixel] = max(matmaxpeakamp,[],3);
    %     for x = 1: size(psth{1},1)
    %         for y = 1:size(psth{1},2)
    %             [bfpixelamp(x,y),bfpixel(x,y)] = max(cell2mat(cellfun(@(z) max(z(x,y,11:15)./mean(z(x,y,5:10),3),[],3),psth,'UniformOutput',0)));
    %
    %         end
    %     end
    
    %     minamp = min(reshape(bfpixelamp,1,size(bfpixelamp,1)*size(bfpixelamp,2)));
    %     maxamp = max(reshape(bfpixelamp,1,size(bfpixelamp,1)*size(bfpixelamp,2)));
    %     meanamp = mean(reshape(bfpixelamp,1,size(bfpixelamp,1)*size(bfpixelamp,2)));
    %     tempbfpixelamp = bfpixelamp;
    %     tempbfpixelamp(bfpixelamp>2)=2; % fix the max so the color is not too white
    %     bfpixelampnorm = ((tempbfpixelamp-minamp)/(2-minamp));
    % %     bfpixelampnorm(bfpixelampnorm<0)=0;
    % %     bfpixelampnorm(bfpixelampnorm>1.2)=1.2;
    %     bfpixelampnormcol=bfpixelampnorm*257;
    %     bfpixelampnormcol(bfpixelampnormcol<=0) = 1;
    %
    intensitydouble = double(squeeze(mean(intensity,3)));
    intensitynorm = (intensitydouble-min(min(intensitydouble)))/(max(max(intensitydouble))--min(min(intensitydouble)));
    intensitynorm(intensitynorm<0.1)=0;
    intensitynorm(intensitynorm>=0.1)=1;
    
    bfpixel = bfpixel.*intensitynorm';
    %         bfpixel(bfpixel==0) = nan;
    h25=figure(25); hold on;
    %     imagesc(1:size(bfpixel,2),1:size(bfpixel,1),(bfpixel')); colormap([[1,1,1];colormaptone;]);
    imagesc(1:size(bfpixel,2),1:size(bfpixel,1),modefilt(bfpixel',[7,7])); colormap([[1,1,1];colormaptone;]);caxis([0,10]);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    % imagesc(1:size(bfpixel,2),1:size(bfpixel,1),imgaussfilt(bfpixel',1.2,'FilterDomain','spatial')); colormap([[1,1,1];colormaptone]);     set(gca,'xtick',[]); set(gca,'ytick',[]);
    %     set(gca,'YDir','reverse'); axis tight
    saveas(h25,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(uniquesite(fn)) '\Movie\site' num2str(uniquesite(fn)) '_tonotopymoviemap.pdf'])
    saveas(h25,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(uniquesite(fn)) '\Movie\site' num2str(uniquesite(fn)) '_tonotopymoviemap.svg'])
    
    
    if numel(intersect(manualroisdone,analysefilelistind(fn)))>0
        col = [0.5,0.5,0.5];
        
        load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_TC_plane0.mat']);
        load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
        
        for cc = 1:size(TC,1)
            for ii = 1:length(startstim)
                roiMatevoked(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40)./median(TC(cc,startstim(ii)-10:startstim(ii)-1),2);
            end
        end
        % cc is roi, ii is tone id, 3rd dimension is trial nb, 4th dimension is
        % time
        responsewindow = 11:15;
        for cc = 1:size(TC,1)
            for ii = 1:length(uniquetone)
                roiMatevokedstim(cc,ii,:,:) = squeeze(roiMatevoked(cc,find(toneid == uniquetone(ii)),:));%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
                roiMatevokedtest{cc}(ii,:) = mean(squeeze(roiMatevoked(cc,reshape(ordertone(ii,:,:),1,length(uniqueduration)*nbrep),responsewindow)),2)'; % % test repmat([1:5]',1,10)' according to this first 10
            end
        end
        for cc =  1:size(TC,1)
            for ii = 1:length(uniquetone)
                meantuningdf{cc}(ii,:) = squeeze(mean(roiMatevokedstim(cc,ii,:,responsewindow),4));
                meanbaselinedf{cc}(ii,:) = squeeze(mean(roiMatevokedstim(cc,ii,:,5:10),4));
                mediantuningdf{cc}(ii,:) = squeeze(median(roiMatevokedstim(cc,ii,:,responsewindow),4));
                stdtuningdf{cc}(ii,:) = std(squeeze(roiMatevokedstim(cc,ii,:,responsewindow))');
            end
        end
        
        
        [~,signibaseline] = cellfun(@(x,y) ttest(reshape(x,size(x,1)*size(x,2),1),reshape(y,size(y,1)*size(y,2),1)),meanbaselinedf,meantuningdf);
        for cc =  1:size(TC,1)
            [Pval(cc,:),~,stats(cc)] = anova2(roiMatevokedtest{cc}',nbrep,'off'); % tones are columns and duration is repeated in 10 trials
        end
        %     signicelltuning = ;%cellfun(@(x) friedman(x',nbrep,'off'),meantuningdf);
        tunedcells = find(signibaseline<0.05);
        freqmaxresp = (cellfun(@(x) find(mean(x,2)==max(mean(x,2))),meantuningdf)); % only for single psth neuron purpose
        bftunedcellsind = (cellfun(@(x) find(mean(x,2)==max(mean(x,2))),meantuningdf(tunedcells)));
        bftunedcells = uniquetone(cellfun(@(x) find(mean(x,2)==max(mean(x,2))),meantuningdf(tunedcells)));
        h20=figure(20); hold on;
        imagesc(1:size(bfpixel,1),1:size(bfpixel,2),modefilt(bfpixel',[7,7])); colormap([[1,1,1];colormaptone;]);caxis([0,10]);
        cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor',[1,1,1],'FaceColor',col),roisCoord{1})
        cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor',[1,1,1],'FaceColor',colormaptone(y,:)),roisCoord{1}(tunedcells),mat2cell(bftunedcellsind,1,ones(1,length(bftunedcellsind))))
        axis tight;     set(gca,'xtick',[]); set(gca,'ytick',[]);
        
        saveas(h20,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(uniquesite(fn)) '\Movie\site' num2str(uniquesite(fn)) '_tonotopymoviemapwithrois.pdf'])
        
        h21=figure(21); hold on;
        cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',col),roisCoord{1})
        cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',colormaptone(y,:)),roisCoord{1}(tunedcells),mat2cell(bftunedcellsind,1,ones(1,length(bftunedcellsind))))
        axis tight; set(gca,'xtick',[]); set(gca,'ytick',[]); axis off
        
        saveas(h21,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(uniquesite(fn)) '\Population\site' num2str(uniquesite(fn)) '_tonotopyonlyrois.pdf'])
        saveas(h21,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(uniquesite(fn)) '\Population\site' num2str(uniquesite(fn)) '_tonotopyonlyrois.svg'])
        
    else
    end
    %caxis([1,10])
    
    % attempt scatter
    % %     figure(11); hold on;
    %     xcoord = [];%1:size(bfpixel,1);
    %     ycoord =  [];
    %     colorid = [];
    %     for tid = 1:length(psth)
    %         clear xcoordtmp ycoordtmp
    %         [xcoordtmp,ycoordtmp] = find(bfpixel==tid);
    %         colorid = [colorid; repmat(colormaptone(tid,:),numel(find(bfpixel==tid)),1)];
    %         xcoord = [xcoord, xcoordtmp'];
    %         ycoord = [ycoord, ycoordtmp'];
    %     end
    %     scatter((xcoord),(ycoord),1,colorid); axis tight
    
    %     figure(1); hold on;
    %     for x = 200:500%size(psth{1},1)
    %         for y = 200:500%size(psth{1},2)
    %             colorpix = squeeze(TonescaledcolorMap(bfpixel(x,y),round(bfpixelampnormcol(x,y)),:))';
    %             plot(x,y,'.','MarkerSize',3,'Color',colorpix);
    %         end
    %     end
    %     bfpixel = cellfun(@(x,y,z) x(y,z),psth(1),xpixelcoordinates,ypixelcoordinates);
    
    for tid = 1:length(psth)
        psth16=uint16(psth{tid});
        stim=mean(psth16(:,:,11:16),3); %calculate your stimulus mean
        base=mean(psth16(:,:,5:10),3); %calculate your baseline mean
        ratio=stim./base; %create a ratio
        %pixel-baesd visualization of a ratio
        %convert your intensity image to a double and adjust the contrast
        intensity = double(squeeze(mean(intensity,3)));
        intensity = intensity/max(max(intensity));
        intensity = imadjust(intensity);
        %set the threshold of your df/F
        minval = 0.9;
        maxval = 1.1;
        %convert ratio into an 8-bit image
        ratiouint8=uint8(((ratio-minval)/(maxval-minval))*255);
        cmap = jet(256); %placeholder colormap, can adjust as needed
        ratiorgb=ind2rgb(ratiouint8,cmap); %convert the indexed image to RGB
        ratiohsv=rgb2hsv(ratiorgb); %convert the RGB image to HSV
        ratiohsv(:,:,3)=intensity'; %map the intensity image to the "value" of HSV
        newratio=hsv2rgb(ratiohsv); %new image keeps the color from the RGB image and the intensity from the intensity image
        h3 = figure(1); hold on; subplot(3,length(psth),[tid,tid+length(psth)]); imagesc(newratio); axis square; set(gca,'xtick',[]); set(gca,'ytick',[]);
        title([num2str(uniquetone(tid))])
        subplot(3,length(psth),tid+length(psth)*2); hold on; plot(xvalues2,squeeze(mean(mean(double(psth16)./double(base))))); plot([0,0],[0.9,1.1],'--k');ylim ([0.9,1.1]);
        title([num2str(uniquetone(tid))])
        
        %axis square;
        
        
        
        minval = 1;
        maxval = 1.2;
        %convert ratio into an 8-bit image
        ratiouint8=uint8(((ratio-minval)/(maxval-minval))*255);
        h = figure(tid+8);% hold on; subplot(2,length(psth)/2,[tid]);
        if tid == 1
            imagesc(intensity);
            %     colormap gray;
        end
        % [a,b]=find(ratio>maxval);
        WhiteColorMap = ones(1,3)*1;%[(zeros(0, 132)), linspace(0, 1, 124)];
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
        %colorMap = flipud([linspace(ToneColorMap(:,1),WhiteColorMap(:,1),256);linspace(ToneColorMap(:,2),WhiteColorMap(:,2),256);linspace(ToneColorMap(:,3),WhiteColorMap(:,3),256)]');%flipud([redColorMap; (zeros(1, 256));% blueColorMap]')'; % from gray to tone
        % blueColorMap = [(zeros(1, 132)), linspace(0, 1, 124)];
        % redColorMap = [linspace(1, 0, 124), (zeros(1, 132))];
        % colorMap = flipud([redColorMap; (zeros(1, 256)); blueColorMap]');
        
        cmap = colorMap;%polarmap(256);%jet(256);
        
        % cmap = jet(256);
        % load('c:\Users\kuchik01\Dropbox\KishoreDropbox\ActivePassive\colormaps_for_figures\newmap.mat');
        % cmap=newmap;
        %load the appropriate colormap
        % load('Z:\bacskai\Michal\MATLAB programs michal\colormaps\finalcm.mat');
        % cmap=cm;
        
        
        %For the first image in the sequence, assigns the colormap from above to
        % the image (to save time afterwards, however, it uses the same colormap for
        % all subsequent images)
        % if tid == 1
        ratiorgb=ind2rgb(ratiouint8,cmap); %convert the indexed image to RGB
        ratiohsv=rgb2hsv(ratiorgb); %convert the RGB image to HSV
        ratiohsv(:,:,3)=intensity'; %map the intensity image to the "value" of HSV
        newratio=hsv2rgb(ratiohsv);
        imshow(permute(newratio,[2 1 3])); axis tight;
        colormap(cmap);
        % cmap=colormap;
        %     else cmap=cmap;
        % end
        %
        % ratiorgb=ind2rgb(ratiouint8,cmap); %convert the indexed image to RGB
        % ratiohsv=rgb2hsv(ratiorgb); %convert the RGB image to HSV
        % if xbitdepth == 16 %olympus fix because 16 bit images only have 12 bit intensity range
        %     xbitdepth=12;
        % end;
        % intensity=double(imadjust(intensity));
        % intensity=intensity/(2^xbitdepth-1); %normalize intensity to maximum possible of yellow + blue channels
        % ratiohsv(:,:,3)=intensity'; %map the intensity image to the "value" of HSV
        % newratio=hsv2rgb(ratiohsv);
        % if tid ==1
        % disk = 'E:';
        % javaaddpath ([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\mij.jar'])
        % javaaddpath ([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\ij-1.52i.jar'])
        % MIJ.start([disk '\KishoreLab\Shared\Matlab\Fiji.app'])
        % end
        % MIJ.createImage(newratio(:,:,3))
        % end
        saveas(h,['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\' 'intensityprojection_' num2str(tid) '.tif'])
        
        h2= figure(); hold on;
        imshow(permute(ratiorgb,[2 1 3])); caxis([1,1.1]); axis tight%colormap(cmap);
        
        
        
        saveas(h2,['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\' 'ratiorgb_' num2str(tid) '.tif'])
        
        if tid == length(psth)
            h= figure(21); imshow(permute(newratio(:,:,1),[2,1,3]));
            saveas(h,['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\' 'allintensityprojection_' num2str(tid) '.tif'])
        end
        
        % aggregate image
        % aggimgrgb(:,:,tid) = squeeze(newratio(:,:,1));
        % save(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\' 'allrgb_' num2str(tid) '.tif'],'aggimgrgb')
        
        % figure(25); hold on;
        % for ii = 1:10
        % clear im
        %     im = imshow(aggimgrgb(:,:,ii));
        %     imAlphaData = repmat(0:1/size(im,2):1-1/size(im,2),size(im,1),1);
        %     set(im,'AlphaData',imAlphaData);
        % end
    end
    saveas(h3,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(uniquesite(fn)) '\Movie\site' num2str(uniquesite(fn)) '_tonotopymovieaverage.pdf'])
end

close all
end