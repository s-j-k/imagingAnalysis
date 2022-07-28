function BatFrameRateRoisConversion()

animal = 'RoofBuddy2';
datatable = readtable(['N:\Jenni\' animal '\ProcessingProgress.xlsx']);

manualroisdone = find(table2array(datatable(:,13))==1);%find(~cellfun(@isempty,cellfun(@(x) find(x=='1'),table2array(datatable(:,13)),'UniformOutput',0)));%find(table2array(datatable(:,13))==1);
daytraining = table2array(datatable(:,1));
siteind = table2array(datatable(:,5));
preprocessedrois = unique(siteind(manualroisdone));
framerate = table2array(datatable(:,8));
protocol = table2array(datatable(:,7));
path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,10));
tonotopyfilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'tonotopy')>0,behaviorfile,'UniformOutput',0))));
sweepfilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'LinearSweep')>0,behaviorfile,'UniformOutput',0))));
angiefilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Angie')>0,behaviorfile,'UniformOutput',0))));
echofilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Echo')>0,behaviorfile,'UniformOutput',0))));
melfilesind =  find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Mel')>0,behaviorfile,'UniformOutput',0))));

lowframeind = find(framerate==31.25);

for ii = 1:length(preprocessedrois)
    behaviorfileind{ii} = intersect(find(siteind==preprocessedrois(ii)),lowframeind);
    analysefilelistind(ii) = find(siteind==preprocessedrois(ii),1);
    tonotopyind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),tonotopyfilesind);%find(daytraining==preprocessedrois(ii),1),find(behaviorfile)); % change this to compare the name of the file and 'Mel' strcmp()
    sweepind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),sweepfilesind);
    if numel(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind))==0
        angieind(ii) = nan;
    else
        angieind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),angiefilesind);
    end
    echoind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),echofilesind);
    if numel(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind))>1 % case of when there was 2 mel files
        melind(ii) = max(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind));
    elseif numel(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind))==0
        melind(ii) = nan;
    else
        melind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind);
    end
end



for fn = 1%2%1:length(analysefilelistind)
    close all
    
    clear TC roiMatevoked tcnorm tunedcells signicelltuning meantuningdf roiMatevokedstim
    % read frames from h5
    h5list = dir(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' '*.h5']);
    nbframestmp(1) = 0;
    for filenb = 1:length(h5list)
        clear hinfo
        hinfo = hdf5info(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' h5list(filenb).name]);
        nbframestmp(filenb+1) = hinfo.GroupHierarchy.Datasets.Dims(3);
    end
    
    % load behavior file
    load(['N:\Jenni\' animal '\behavior\' behaviorfile{tonotopyind(fn)} '.mat']); %RoofBuddy1tonotopy02-03-2021-12-50
    startstim = batimagingdata.frametone;%1:400:44000;
    
    seqid = batimagingdata.durationseq;%batimagingdata.;%repmat(1:11,1,10);
    toneid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
    uniquetone = unique(toneid);
    uniqueduration = unique(seqid); % also order of presentation
    for tid = 1:length(uniquetone)
        for dd = 1:length(uniqueduration)
            ordertone(tid,dd,:) = intersect(find(toneid==uniquetone(tid)),find(seqid==uniqueduration(dd))); % tone id duration and trial rep
        end
    end
    % nb of repetition
    nbrep = (length(seqid)/length(uniquetone))/length(uniqueduration);
    
    xvalues = [-10:40]/31.25;
    % load ROIs
    load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_TC_plane0.mat']);
    load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
    load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\Fall.mat'])
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\TonotopyMeanEvokedMovie.mat'])
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\TonotopyDurationMeanEvokedMovie.mat'])
    responsewindow = 11:15;
    
    
end


end