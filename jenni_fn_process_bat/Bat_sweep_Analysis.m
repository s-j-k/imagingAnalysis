function Bat_tonotopy_Analysis()
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



for fn = 6:8%1:length(analysefilelistind)%6:7 for RoofBuddy2
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
    
    % load behavior file
    load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\behavior\' behaviorfile{sweepind(fn)} '.mat']); %RoofBuddy1tonotopy02-03-2021-12-50
    startstim = batimagingdata.frametone;%1:400:44000;
    
    seqid = batimagingdata.durationseq;%batimagingdata.;%repmat(1:11,1,10);
    toneid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
    uniquetone = unique(toneid); % 1 is upsweep, 2 is downsweep, 3 is white noise
    uniqueduration = unique(seqid); % also order of presentation
    stimlist = {'white noise','upsweep','downsweep'};
    if ~isempty(intersect(stimselection,analysefilelistind(fn))) % remove tone above 10ms if stim needs to be corrected
        uniqueduration(uniqueduration>=10) = [];
        startstim(seqid>6)=[];
        toneid(seqid>6)=[];
        seqid(seqid>6)=[];
    else
    end
    
    for tid = 1:length(uniquetone)
        for dd = 1:length(uniqueduration)
            ordertone(tid,dd,:) = intersect(find(toneid==uniquetone(tid)),find(seqid==uniqueduration(dd))); % tone id duration and trial rep
        end
    end
    % nb of repetition
    nbrep = (length(toneid)/length(uniquetone))/length(uniqueduration);
    
    xvalues = [-10:40]/31.25;
    % load ROIs
    load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_TC_plane0.mat']);
    load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
    load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\Fall.mat'])
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\TonotopyMeanEvokedMovie.mat'])
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\TonotopyDurationMeanEvokedMovie.mat'])
    responsewindow = 11:25; % from 0 to 0.5s (roughly the tail end of the mean psths)
    
    startstim = startstim+nbframestmp((sweepind(fn)-analysefilelistind(fn))+1);
    %TC = TC;%./median(TC')';
    for cc = 1:size(TC,1)
        for ii = 1:length(startstim)
            roiMatevoked(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40)./median(TC(cc,startstim(ii)-10:startstim(ii)-1),2);
            % roiMatevokednorm(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40)./median(TC(cc,startstim(ii)-10:startstim(ii)),2);
        end
    end
    % remove duration corrupted > than 6ms (10 and 20)
    
    % cc is roi, ii is tone id, 3rd dimension is trial nb, 4th dimension is
    % time
    for cc = 1:size(TC,1)
        roiMatevokedtest{cc} = [];
        for ii = 1:length(uniquetone)
            roiMatevokedstim(cc,ii,:,:) = squeeze(roiMatevoked(cc,find(toneid == uniquetone(ii)),:));%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
            %             for dd = 1:length(uniqueduration)
            roiMatevokedtest{cc}(ii,:) = mean(squeeze(roiMatevoked(cc,reshape(ordertone(ii,:,:),1,length(uniqueduration)*nbrep),responsewindow)),2)'; % % test repmat([1:5]',1,10)' according to this first 10
            %             trials are duration1 then 10 trials for duration 2 ect...
            %             end
        end
    end
    
    
    
    
    % try testing duration and tuning with anova2
    
    
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
    
    ylabel('Median df/f'); xlabel('Time from tone onset [s]');     xlim([xvalues(1),xvalues(36)]);
    set(gca,'xtick',[xvalues(11),xvalues(27)]);
    set(gca,'xticklabels',{'0','0.5'});
    set(gca,'fontsize',12); set(gca,'fontname','arial')
    
    %
%     print(['K:\Jenni\' animal '\figure\sweep\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_MeanEvokedResponse'],'-dpdf')
%     %     print(['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_MeanEvokedResponse'],'-dsvg')
%     saveas(h,['K:\Jenni\' animal '\figure\sweep\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_MeanEvokedResponse.svg'])  
  
    % Plot population PSTH for wn,us,ds
    rangeoct = 4.5; % from 4000 to 90000 in 0.5 oct space % 4000*2^4.5 = 9.0510e+04
    uniquerate = rangeoct./(uniqueduration/1000);% octave range is 4.5 always but in diff length unique duration is in ms so divide by 1000
    uniqueratekhz = round(((90.510-4)./uniqueduration)*10)./10; % kHZ inmn ms to comapre to Morrisson 2018 paper
    h = figure(3); hold on;
    for stid = 1:length(uniquetone)
        subplot(3,3,stid);
        imagesc(xvalues,1:size(roiMatevokedstim,1),squeeze(median(roiMatevokedstim(:,stid,:,:),3))); colormap gray; %caxis([1,1.1])
        ylabel('ROIs'); caxis([0.95,1.15])
        title(stimlist{stid})
        %   plot diff durations
        subplot(3,3,stid+(length(uniquetone))); hold on;
        for dt = 1:length(uniqueduration)
            PSTHdur = (squeeze(median(roiMatevoked(:,intersect(find(toneid==stid),find(seqid == uniqueduration(dt))),:),2)));
            a = plot(xvalues,mean(PSTHdur),'Color',[0,0,0]+0.1*dt);
            ylabel('df/f'); ylim([1,1.2]);
            amprate(stid,dt) = mean(max(squeeze(median(roiMatevoked(:,intersect(find(toneid==stid),find(seqid == uniqueduration(dt))),responsewindow),2))'));
            stdrate(stid,dt) = std(max(squeeze(median(roiMatevoked(:,intersect(find(toneid==stid),find(seqid == uniqueduration(dt))),responsewindow),2))'));
            %     Responsedurationnorm(stid,dt) = sum((PSTHdur(:,11:40))./(max(PSTHdur(:,11:40)')'));%mean(sum((PSTHdur)./(max(PSTHdur)')')/length(xvalues));
            %     stdResponsedurationnorm(stid,dt) = std(sum((squeeze(median(roiMatevoked(:,intersect(find(toneid==stid),find(seqid == uniqueduration(dt))),11:30),2))'))./(max(squeeze(median(roiMatevoked(:,intersect(find(toneid==stid),find(seqid == uniqueduration(dt))),:),2))')));
        end
        if stid == 1 % duration for white noise (1) and oct/s in 
            legend(cellfun(@(x) num2str(x), mat2cell(uniqueduration,1,ones(1,length(uniqueduration))),'UniformOutput',0)); legend box off
        else
            legend(cellfun(@(x) num2str(x), mat2cell(uniqueratekhz,1,ones(1,length(uniqueduration))),'UniformOutput',0)); legend box off
        end

        subplot(3,3,(stid+(length(uniquetone))*2)); hold on;
        if stid == 1
        shadedErrorBar(uniqueduration,amprate(stid,:),stdrate(stid,:),{'Color',[0,0,0],'LineWidth',2},1);     
        xticks(uniqueduration);
        xlim([uniqueduration(1),uniqueduration(end)])
        xticklabels(mat2cell(uniqueduration,1,ones(1,length(uniqueduration))))
        xlabel('Duration [ms]')
        else
            shadedErrorBar(uniqueduration,amprate(stid,:),stdrate(stid,:),{'Color',[0,0,0],'LineWidth',2},1);
            xticks(uniqueduration);
            xlim([uniqueduration(1),uniqueduration(end)])
            xticklabels(mat2cell((uniqueratekhz),1,ones(1,length(uniqueduration))))
            xlabel('Rate [kHz/ms]')
        end
        ylabel('Mean df/f')
        %     subplot(4,3,(stid+(length(uniquetone))*3)); hold on;
        %     shadedErrorBar(uniquerate,Responsedurationnorm(stid,:),stdResponsedurationnorm(stid,:),{'Color',[0.5,0.2,0.5],'LineWidth',2},1);
        %     xticks(uniquerate);
        %     xlim([uniquerate(1),uniquerate(end)])
        %     xticklabels(mat2cell(uniquerate,1,ones(1,length(uniqueduration))))
    end
    %     set(gca,'fontsize',12); set(gca,'fontname','arial')
    saveas(h,['N:\Jenni\' animal '\figure\sweep\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_SweepRateSummary.pdf'])
    print(['N:\Jenni\' animal '\figure\sweep\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_SweepRateSummary'],'-dpdf')
    % look at the shape of the 'direction' selectivity for all cells
    
    figure(4); hold on;
    for stid = 1:length(uniquetone)
        for dt = 1:length(uniqueduration)
            PSTHdur = (squeeze(median(roiMatevoked(:,intersect(find(toneid==stid),find(seqid == uniqueduration(dt))),:),2)));
            ampcells(stid,dt,:) = max(PSTHdur(:,responsewindow)');
            peaklatencycells(stid,dt,:) = cellfun(@(x,y) find(x==y),mat2cell((PSTHdur(:,responsewindow)'),length(responsewindow),ones(1,size(ampcells,3))),num2cell(max(PSTHdur(:,responsewindow)')));
            %         ampcells = max(PSTHdur(:,responsewindow)');            
        end
        subplot(2,3,stid); hold on;
        imagesc(1:length(uniqueduration),1:size(ampcells,3),squeeze(ampcells(stid,:,:))'); colormap gray; caxis([1,1.3])
        ylim([1,size(ampcells,3)]);
        xticks(1:length(uniqueduration));
        if stid ==1
            xticklabels(mat2cell(uniqueduration,1,ones(1,length(uniqueduration))))
            xlabel('Duration [ms]')
        else
            xticklabels(mat2cell((uniqueratekhz),1,ones(1,length(uniqueduration))))
            xlabel('Rate [kHz/ms]')
        end
%         xlim([1,length(uniqueduration)])
    end
    % compute direction selectivity
    % first upsweep downsweep
%     redgreenmapdur = copper(length(uniqueduration));
%     figure(5); hold on; subplot(2,2,1); hold on; 
%     for dt = 1:length(uniqueduration)
%         sweepselectivity(dt,:) = (ampcells(3,dt,:)-ampcells(2,dt,:))./(ampcells(3,dt,:)+ampcells(2,dt,:));
%         [sortedselectivityvalue(dt,:),sortedselectivity(dt,:)] = sort(sweepselectivity(dt,:));
%         plot(sortedselectivityvalue(dt,:),1:size(sweepselectivity,2),'Color',redgreenmapdur(dt,:));% ,'Color',[0,0,0]+0.1*dt
%         % intersection point 
%         intersectionpoint(dt,:) = find(sortedselectivityvalue(dt,:)>0,1);
%     end
%    legend(cellfun(@(x) num2str(x), mat2cell(uniqueratekhz,1,ones(1,length(uniqueduration))),'UniformOutput',0)); legend box off
%     xlim([min(min(sortedselectivityvalue)),max(max(sortedselectivityvalue))]); 
%     plot([0,0],[1,size(sortedselectivityvalue,2)],'k--')
%     ylim([1,size(sortedselectivityvalue,2)]);
%     subplot(2,2,2); hold on;
%         figure(6); hold on; 
%         for dt = 1:length(uniqueduration)
%         subplot(1,8,dt); hold on;
%         PSTHwn(dt,:,:) = (squeeze(median(roiMatevoked(:,intersect(find(toneid==1),find(seqid == uniqueduration(dt))),:),2)));
%         imagesc(squeeze(PSTHwn(dt,:,:)));colormap; caxis([1,1.2]);
%         end
    % plot gross selectivity
%     hist(reshape(sweepselectivity,1,size(sweepselectivity,1)*size(sweepselectivity,2)))
%     
%     
    Sweepstruct(siteind(analysefilelistind(fn))).ampcells =  ampcells;
    Sweepstruct(siteind(analysefilelistind(fn))).peaklatency =  peaklatencycells;
    Sweepstruct(siteind(analysefilelistind(fn))).uniquetone =  uniquetone;
    Sweepstruct(siteind(analysefilelistind(fn))).uniqueduration =  uniqueduration;
    Sweepstruct(siteind(analysefilelistind(fn))).uniquerate =  uniquerate;
    Sweepstruct(siteind(analysefilelistind(fn))).uniqueratekhz =  uniqueratekhz;
    Sweepstruct(siteind(analysefilelistind(fn))).nbrep = nbrep;
%     if strcmp(animal, 'RoofBuddy2') || strcmp(animal, 'RoofBuddy1')
%         Sweepstruct(siteind(analysefilelistind(fn))).depth = mediandeppthrois;
%     end
%     %     end
% end
end
save(['N:\Jenni\' animal '\z-stack\Sweepstruct.mat'],'Sweepstruct')

disp('')
