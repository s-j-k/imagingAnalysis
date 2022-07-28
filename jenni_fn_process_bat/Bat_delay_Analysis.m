function Bat_delay_Analysis()
animal = 'RoofBuddy2';%'RoofBuddy1';
datatable = readtable(['N:\Jenni\' animal '\ProcessingProgress.xlsx']);

manualroisdone = find(table2array(datatable(:,13))==1);%find(~cellfun(@isempty,cellfun(@(x) find(x=='1'),table2array(datatable(:,13)),'UniformOutput',0)));%find(table2array(datatable(:,13))==1);
daytraining = table2array(datatable(:,1));
siteind = table2array(datatable(:,5));
preprocessedrois = unique(siteind(manualroisdone));
framerate = table2array(datatable(:,8));
protocol = table2array(datatable(:,7));
path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,10));
stimselection = find(table2array(datatable(:,14))==1);
datadrivecolumn = (table2array(datatable(:,15)));
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



for fn = 6%1:3%1:length(analysefilelistind)%6:7 for RoofBuddy2
    close all
    
    clear TC ampcells peaklatencycells roiMatevoked ordertone depthvalue Pval stats tcnorm tunedcells tuningrespcells signicelltuning meantuningdf roiMatevokedstim meanbaselinedf dfmaxvalue dfmaxind mediandeppthrois depthroi mediantuningdf stdtuningdf
    % read frames from h5
    h5list = dir([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' '*.h5']);
    nbframestmp(1) = 0;
    for filenb = 1:length(h5list)
        clear hinfo
        hinfo = hdf5info([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' h5list(filenb).name]);
        nbframestmp(filenb+1) = hinfo.GroupHierarchy.Datasets.Dims(3);
    end
    cumsumframe = cumsum(nbframestmp);
    % load behavior file
    load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\behavior\' behaviorfile{echoind(fn)} '.mat']); %RoofBuddy1tonotopy02-03-2021-12-50
    startstim = batimagingdata.frametone;
    startstim = startstim+cumsumframe((echoind(fn)-analysefilelistind(fn))+1);
    
    seqid = batimagingdata.durationseq;%batimagingdata.;%repmat(1:11,1,10);
    toneid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
    uniquetone = unique(toneid); % check wav names for delay 
    uniqueduration = unique(seqid); % also order of presentation
    if length(batimagingdata.wavnames)==22
    stimlist = [num2cell(cell2mat(cellfun(@(x) sscanf(x,'PE%d'),batimagingdata.wavnames,'UniformOutput',0))),{'P'},{'E'},num2cell(cell2mat(cellfun(@(x) sscanf(x,'PE%d'),batimagingdata.wavnames,'UniformOutput',0)))];
    else
    end
    
    if ~isempty(intersect(stimselection,analysefilelistind(fn))) % remove tone above 10ms if stim needs to be corrected
        uniqueduration(uniqueduration>=10) = [];
        startstim(seqid>6)=[];
        toneid(seqid>6)=[];
        seqid(seqid>6)=[];
    else
    end
    % get stim id PE 1:11 echoduration, then single then reverse
    stimid = [num2cell(cell2mat(cellfun(@(x) sscanf(x,'PE%d'),batimagingdata.wavnames,'UniformOutput',0))),'simple call','simple echo',num2cell(cell2mat(cellfun(@(x) sscanf(x,'PE%d'),batimagingdata.wavnames,'UniformOutput',0)))];
    for tid = 1:length(uniquetone)
            ordertone(tid,:) = (find(toneid==uniquetone(tid))); % tone id duration and trial rep
    end
    % nb of repetition
    nbrep = (length(toneid)/length(uniquetone));
    
    xvalues = [-10:40]/31.25;
    % load ROIs
    load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_TC_plane0.mat']);
    load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
    load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\Fall.mat'])
    responsewindow = 11:25; % from 0 to 0.5s (roughly the tail end of the mean psths)
    
    for cc = 1:size(TC,1)
        for ii = 1:length(startstim)
            roiMatevoked(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40)./median(TC(cc,:),2); % startstim(ii)-10:startstim(ii)-1) that was for baseline
        end
    end
    
    % cc is roi, ii is tone id, 3rd dimension is trial nb, 4th dimension is
    % time
    
    for cc = 1:size(TC,1)
        roiMatevokedtest{cc} = [];
        for ii = 1:length(uniquetone)
            roiMatevokedstim(cc,ii,:,:) = squeeze(roiMatevoked(cc,find(toneid == uniquetone(ii)),:));%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
        end
    end
    
    
    
    
    % try testing duration and tuning with anova2
     % identify sound evoked cells per cell paired ttest or anova 2
    for cc = 1:size(TC,1)
        baseline = mean(roiMatevoked(cc,:,1:10),3)-1;
        evoked = mean(roiMatevoked(cc,:,11:21),3)-1;
        pttest(cc) = ttest(baseline,evoked);
        baselineanov = size(squeeze(mean(median(roiMatevokedstim(:,:,:,1:10),3),4)));
        evokedanov = size(squeeze(mean(median(roiMatevokedstim(:,:,:,11:21),3),4)));
    end
    evokedcells = find(pttest==1); % sound responsive
    
    % plot mean evoked response per roi for all trials
    [sortedroisamp,sortedroisind] = sort(mean(squeeze(median(roiMatevoked(:,:,11:15),2))'));
    h=figure(2);
    hold on; subplot(3,1,[1,2]);hold on; imagesc(xvalues,1:size(TC,1),squeeze((median(roiMatevoked(sortedroisind,:,:),2)))); colormap gray;
    ylabel('ROI'); caxis([0.95,1.1]); hold on;
    plot([xvalues(11),xvalues(11)],[1,size(TC,1)],'Color',[1,1,1],'Linewidth',1.5); ylim([1,size(TC,1)])
    xlim([xvalues(1),xvalues(36)]);
    set(gca,'xtick',[xvalues(11),xvalues(27)]);
    set(gca,'xticklabels',{'0','0.5'});
    set(gca,'YDir','normal');
    set(gca,'ytick',[100,300,500]);
    set(gca,'yticklabels',{'100','300','500'});
    set(gca,'fontsize',12); set(gca,'fontname','arial')
    
    subplot(3,1,[3]);hold on ;
    hold on; plot(xvalues,(squeeze((mean(roiMatevoked,2)))),'Color',[0.9,0.9,0.9]); xlim([xvalues(1),xvalues(end)]);
    
    hold on; plot(xvalues,median(squeeze((mean(roiMatevoked,2)))),'k','Linewidth',3); xlim([xvalues(1),xvalues(end)]);
    plot([xvalues(11),xvalues(11)],[0.95,1.2],'Color',[0,0,0],'Linewidth',1.5);
    ylim([0.95,1.2]);
%      hold on; plot(xvalues,median(squeeze((mean(roiMatevoked(evokedcells,:,:),2)))),'r','Linewidth',3); xlim([xvalues(1),xvalues(end)]);
%     plot([xvalues(11),xvalues(11)],[0.95,1.2],'Color',[0,0,0],'Linewidth',1.5);
%     ylim([0.95,1.2]);
    ylabel('Median df/f'); xlabel('Time from tone onset [s]');     xlim([xvalues(1),xvalues(36)]);
    set(gca,'xtick',[xvalues(11),xvalues(27)]);
    set(gca,'xticklabels',{'0','0.5'});
    
    
    set(gca,'fontsize',12); set(gca,'fontname','arial')
    
    [delayval,sortdelayind]=sort([stimlist{1:10}]);
   
    
    figure(); hold on;
    for iii = 1:10
        subplot(5,2,iii); hold on; plot(xvalues,mean(squeeze(median(roiMatevoked(:,find(toneid == uniquetone(sortdelayind(iii))),:),2))),'k','LineWidth',2);
        subplot(5,2,iii); hold on; plot(xvalues,mean(squeeze(median(roiMatevoked(:,find(toneid == uniquetone(sortdelayind(iii)+12)),:),2))),'Color',[0.5,0.2,0.5],'LineWidth',2);
        xlim([xvalues(1),xvalues(end)]); 
        ylim([0.95,1.1]);
    end
    % plot "delay tuning"
    for iii = 1:10
        delaytuningmaxdf(:,iii) = squeeze(max(squeeze(median(roiMatevokedstim(:,(sortdelayind(iii)),:,responsewindow),3))'));
        for cc= 1:size(TC,1)
            for rep= 1:10
            maxdftuningtrials(cc,iii,rep) = max(squeeze(squeeze(roiMatevokedstim(cc,(sortdelayind(iii)),rep,responsewindow))));
%             figure(); hold on;
%             if cc= 346
%                 plot(squeeze(median(roiMatevokedstim(346,(sortdelayind(iii)),:,responsewindow),3)),'Color',[0,0,0]+((0.09)*iii));
%             end
            end
        end
    end
    [delaymaxval,delaymaxvalind] = max(delaytuningmaxdf');
    [sortmaxdelaycellsdelayval,sortmaxdelaycells] = sort(delaymaxvalind(evokedcells));
    delaytuningmaxdfnorm = (delaytuningmaxdf-(min(delaytuningmaxdf')'))./((max(delaytuningmaxdf')')-(min(delaytuningmaxdf')'));
    % only plot sound evoked cells
    figure(); hold on; subplot(1,2,1); hold on; imagesc(delayval,1:size(roiMatevoked,1),delaytuningmaxdf(evokedcells(sortmaxdelaycells),:)); colormap gray
    axis tight
    subplot(1,2,2); hold on; imagesc(delayval,1:size(roiMatevoked,1),delaytuningmaxdfnorm(evokedcells(sortmaxdelaycells),:)); colormap gray
    axis tight
    % plot the delay tuned cells (although need to do an anova to check)
    figure(); hold on; 
    for ii = 1:length(delayval)
    subplot(2,5,ii); plot(delayval,delaytuningmaxdf(evokedcells(sortmaxdelaycells(find(sortmaxdelaycellsdelayval==ii))),:)','Color',[0.5,0.5,0.5],'LineWidth',1);
    end
    % bandwidth
    bandwidth = (sum(delaytuningmaxdfnorm(evokedcells,:)')/10);
    [sortbandval,sortbandind] = sort(bandwidth);
    % plot bandwidth vs delaytuning 
    plot(delaymaxvalind(evokedcells),bandwidth,'.')
    figure();
    hold on; imagesc(xvalues,1:size(evokedcells,1),squeeze((median(roiMatevoked(sortedroisind(ismember(sortedroisind,evokedcells)),:,:),2)))); colormap gray;
    ylabel('ROI'); caxis([0.95,1.1]); hold on; axis tight
    % test further to see if cells have tuning kruskalwallis on
    % delaymaxnorm
    for cc = 1:size(TC,1)
    vecval = reshape(squeeze(maxdftuningtrials(cc,:,:)),1,10*10); % 10 rep of  each delay
    grp = reshape(repmat((1:10),10,1),1,10*10);
   [pfrid(cc)] = anova1(vecval,grp,'off'); % columns needs to be delay
    end

     subplot(1,3,3); hold on; plot(mean(squeeze(median(roiMatevoked(:,find(toneid == uniquetone(11)),:),2))),'r');      
     subplot(1,3,3); hold on; plot(mean(squeeze(median(roiMatevoked(:,find(toneid == uniquetone(12)),:),2))),'b');      

    
  
    Delaystruct(siteind(analysefilelistind(fn))).ampcells =  ampcells;
    Delaystruct(siteind(analysefilelistind(fn))).peaklatency =  peaklatencycells;
    Delaystruct(siteind(analysefilelistind(fn))).uniquetone =  uniquetone;
    Delaystruct(siteind(analysefilelistind(fn))).uniqueduration =  uniqueduration;
    Delaystruct(siteind(analysefilelistind(fn))).uniquerate =  uniquerate;
    Delaystruct(siteind(analysefilelistind(fn))).uniqueratekhz =  uniqueratekhz;
    Delaystruct(siteind(analysefilelistind(fn))).nbrep = nbrep;

end
save(['N:\Jenni\' animal '\z-stack\Sweepstruct.mat'],'Sweepstruct')

disp('')
