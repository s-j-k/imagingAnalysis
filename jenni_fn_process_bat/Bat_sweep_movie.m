function Bat_sweep_movie()
animal = 'RoofBuddy1';
datatable = readtable(['N:\Jenni\' animal '\ProcessingProgress.xlsx']);

manualroisdone = find(table2array(datatable(:,13))==1);%find(~cellfun(@isempty,cellfun(@(x) find(x=='1'),table2array(datatable(:,13)),'UniformOutput',0)));%
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

for fn = 11:length(analysefilelistind)%1:length(analysefilelistind)
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
    load(['N:\Jenni\' animal '\behavior\' behaviorfile{intersect(intersect(find(site == uniquesite(fn)),find(~cellfun(@isempty,cellfun(@(x) strfind(x,'LinearSweep')>0,behaviorfile,'UniformOutput',0)))'),find(framerate==31.25))} '.mat']); %RoofBuddy1tonotopy02-03-2021-12-50
    startstim = batimagingdata.frametone;%1:400:44000;
    
    seqid = batimagingdata.durationseq;%batimagingdata.;%repmat(1:11,1,10);
    toneid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
    uniquetone = unique(toneid);
    uniqueduration = unique(seqid); % also order of presentation
    % nb of repetition
    nbrep = (length(seqid)/length(uniquetone))/length(uniqueduration);
    xvalues = [-10:40]/31.25;
    % load ROIs
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_TC_plane0.mat']);
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\Fall.mat'])
    load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\LinearMeanEvokedMovie.mat'])
    matduration = load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\LinearDurationMeanEvokedMovie.mat']);

    minval = 1;% prctile(reshape(ratio,1,size(ratio,1)*size(ratio,2)),[5]);
    maxval = 1.3;%prctile(reshape(ratio,1,size(ratio,1)*size(ratio,2)),[95]); across all trials
    psth=MeanEvokedMovie; %average across all trials
    %     intersity2 = mean(MeanEvokedMovie{1},3);
    xvalues2 = ([-10:44]*31.25)/1000;
    if exist((['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2pmeanimg.tif']))
        intensity=imread(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2pmeanimg.tif']); %load your intensity/mean image
    else
        load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\Fall.mat'],'ops')
        intensity=uint16(ops.meanImg); %load your intensity/mean image
        %         meanimg = imadjust(uint16(ops.meanImg));
        %         h=figure(2); hold on; imshow(meanimg); axis tight;
        %         saveas(h,['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2pmeanimg.tif'])
    end
    colromaptone = jet(length(psth));
    stimlist = {'White Noise', 'Linear Upsweep', 'Linear Downsweep'};
    for tid = 1:length(psth)
        psth16=uint16(psth{tid});
        stim=mean(psth16(:,:,11:16),3); %calculate your stimulus mean
        base=mean(psth16(:,:,1:10),3); %calculate your baseline mean
        ratio=stim./base; %create a ratio
        %pixel-baesd visualization of a ratio
        %convert your intensity image to a double and adjust the contrast
        intensity = double(squeeze(mean(intensity,3)));
        intensity = intensity/max(max(intensity));
        intensity = imadjust(intensity);
        %set the threshold of your df/F
        minval = 0.95;
        maxval = 1.12;
        %convert ratio into an 8-bit image
        ratiouint8=uint8(((ratio-minval)/(maxval-minval))*255);
        cmap = jet(256); %placeholder colormap, can adjust as needed
        ratiorgb=ind2rgb(ratiouint8,cmap); %convert the indexed image to RGB
        ratiohsv=rgb2hsv(ratiorgb); %convert the RGB image to HSV
        ratiohsv(:,:,3)=intensity'; %map the intensity image to the "value" of HSV
        newratio=hsv2rgb(ratiohsv); %new image keeps the color from the RGB image and the intensity from the intensity imagexlabel()
        h3 = figure(1); hold on; subplot(3,length(psth),[tid,tid+length(psth)]); imagesc(newratio); axis square; 
        title(stimlist{tid})
        %subplot(3,length(psth),tid+length(psth)*2); hold on; plot(xvalues2,squeeze(mean(mean(psth{tid})))./double(base)))),'k','LineWidth',3); plot([0,0],[0.9,1.1],'--k');ylim ([0.9,1.1]);
        %axis square;
        % add psth duration
        for dt = 1:length(uniqueduration)
            psthduration = matduration.MeanEvokedMovieDuration{tid}(:,:,:,dt);
            basedur=mean(psthduration(:,:,1:10),3); %calculate your baseline mean
            %subplot(3,length(psth),tid+length(psth)*2); hold on; plot(xvalues2,squeeze(mean(mean((psthduration./mean(psthduration(:,:,1:10),3))))),'Color',[1,1,1]-((1/(length(uniqueduration)+1))*dt)); %plot([0,0],[0.9,1.1],'--k');ylim ([0.9,1.1]);
            subplot(3,length(psth),tid+length(psth)*2); hold on; plot(xvalues2,squeeze(mean(mean(double(psthduration)./double(basedur)))),'Color',[1,1,1]-((1/(length(uniqueduration)+1))*dt)); ylim ([0.95,1.2]); xlim([xvalues2(1),xvalues2(end)])
            %plot([0,0],[0.9,1.2],'--k');
        end
        xlabel('Time from sound onset [s.]'); ylabel('df/f')
        legend(cellfun(@(x) num2str(x), mat2cell(uniqueduration,1,ones(1,length(uniqueduration))),'UniformOutput',0))
        legend boxoff
    end
    print(['N:\Jenni\' animal '\figure\tonotopy\site' num2str(uniquesite(fn)) '\Movie\site' num2str(uniquesite(fn)) '_SweepWNmovieaverage.pdf'],'-dpdf','-fillpage')
%   saveas(h3,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(uniquesite(fn)) '_SweepWNmovieaverage.pdf'])
end
close all
end