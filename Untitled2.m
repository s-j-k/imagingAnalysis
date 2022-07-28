datapath='P:\su\DATA\imagingData\passive\sk023\day4\suite2p\plane0';
behpath='P:\su\DATA\imagingData\passive\sk023\day4\sk023tonotopy04-05-2022-15-48.mat';
cd(datapath)
load('Fall.mat')
load(behpath)

%%
cells=F((iscell==1),:);
startstim = batimagingdata.frametone+23088;%1:400:44000;
    
    seqid = batimagingdata.durationseq;%batimagingdata.;%repmat(1:11,1,10);
    toneid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
    uniquetone = unique(toneid);
    uniqueduration = unique(seqid); % also order of presentation
    
%     if ~isempty(intersect(stimselection,analysefilelistind(fn))) % remove tone above 10ms if stim needs to be corrected
%         uniqueduration(uniqueduration>=10) = [];
%         startstim(seqid>6)=[];
%         toneid(seqid>6)=[];
%         seqid(seqid>6)=[];
%     else
%     end

% % calculate ITIs
% ITIid=[];
% for kk=1:length(toneid)
%     ITIid=toneid(kk)-seqid(kk);
% end
% 
% Attid=batimagingdata.att
% uniqueITI=unique(ITIid);
% uniqueAtt=unique(Attid);

    for tid = 1:length(uniquetone)
        for dd = 1:length(uniqueduration)
            ordertone(tid,dd,:) = intersect(find(toneid==uniquetone(tid)),find(seqid==uniqueduration(dd))); % tone id duration and trial rep
        end
    end
    % nb of repetition
%     nbrep = (length(toneid)/length(uniquetone))/length(uniqueduration);
    nbrep=3;
    xvalues = [-10:40]/31.25;
    % load ROIs
%     load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_TC_plane0.mat']);
%     load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
%     load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\Fall.mat'])
   
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\TonotopyMeanEvokedMovie.mat'])
    %     load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\TonotopyDurationMeanEvokedMovie.mat'])
    responsewindow = 11:15;
    
   %%
    TC = cells;%./median(TC')';
    for cc = 1:size(TC,1)
        for ii = 1:length(startstim)
            roiMatevoked(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40)./median(TC(cc,startstim(ii)-10:startstim(ii)-1),2);
            roiMatevokednorm(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40)./median(TC(cc,startstim(ii)-10:startstim(ii)),2);
        end
    end
    % remove duration corrupted > than 6ms (10 and 20)
    
    %%
    
    % plot individual traces
    cellnum=ll;    ii=3;
    figure;
    ans=roiMatevoked(ll,ii,:);
    ans=squeeze(ans);
    plot(1:51,ans);
    % cc is roi, ii is tone id, 3rd dimension is trial nb, 4th dimension is
    % time
    %%
    
    for cc = 1:size(TC,1)
        for ii = 1:length(uniquetone)
            roiMatevokedstim(cc,ii,:,:) = squeeze(roiMatevoked(cc,find(toneid == uniquetone(ii)),:));%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
            %             for dd = 1:length(uniqueduration)
            roiMatevokedtest{cc}(ii,:) = mean(squeeze(roiMatevoked(cc,reshape(ordertone(ii,:,:),1,[]),responsewindow)),2)'; % % test repmat([1:5]',1,10)' according to this first 10
            roiMatevokeddurt(cc,ii,:,:) = squeeze(roiMatevoked(cc,find(seqid == uniqueduration(ii)),:));
            %             trials are duration1 then 10 trials for duration 2 ect...
            %             end
        end
    end
    
    
    %%
    figure;
    imagesc(squeeze(mean(mean(roiMatevokedstim,2),3)));
    
    % try testing duration and tuning with anova2
    %%
    
    % plot mean evoked response per roi for all trials
    [sortedroisamp,sortedroisind] = sort(squeeze(median(roiMatevoked(:,:,11:15),2))');
    
    h=figure(3);
    
    hold on; subplot(3,1,[1,2]);hold on; imagesc(xvalues,1:size(TC,1),...
        squeeze((median(roiMatevoked(sortedroisind,:,:),2)))); colormap gray;
    ylabel('ROI'); caxis([0.995,1.005]); hold on; 
%     colormap(flipud(gray)); 
    title('Tone evoked response for frequency sweeps and white noise combined') % what the actual fuck
    plot([xvalues(11),xvalues(11)],[1,size(TC,1)],'Color',[1,1,1],'Linewidth',1.5); ylim([1,size(TC,1)])
    xlim([xvalues(1),xvalues(50)]);
    set(gca,'xtick',[xvalues(11),xvalues(27)]);
    set(gca,'xticklabels',{'0','0.5'});
    set(gca,'YDir','normal');
    set(gca,'ytick',[100,300,500]);
    set(gca,'yticklabels',{'100','300','500'});
    set(gca,'fontsize',12); set(gca,'fontname','arial')
% 
    subplot(3,1,[3]);hold on ;
    hold on; plot(xvalues,(squeeze((mean(roiMatevoked,2)))),'Color',[0.9,0.9,0.9]); xlim([xvalues(1),xvalues(end)]);
    
    hold on; plot(xvalues,median(squeeze((mean(roiMatevoked,2)))),'k','Linewidth',3); xlim([xvalues(1),xvalues(end)]);
    plot([xvalues(11),xvalues(11)],[0.95,1.2],'Color',[0,0,0],'Linewidth',1.5);
    ylim([0.95,1.05]);
    
    ylabel('Median df/f'); xlabel('Time from tone onset [s]');     xlim([xvalues(1),xvalues(36)]);
    set(gca,'xtick',[xvalues(11),xvalues(27)]);
    set(gca,'xticklabels',{'0','0.5'});
    set(gca,'fontsize',12); set(gca,'fontname','arial')
%     colormap cool(10)
    %%
    %
    print(['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_MeanEvokedResponse'],'-dpdf')
    %     print(['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_MeanEvokedResponse'],'-dsvg')
    saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_MeanEvokedResponse.svg'])
    %%
  
    figure; hold on;
    cells=F(iscell==1,:);
    cells=cells(tunedcells,:);
nTraces=size(cells,1);
offset=5000;
% for i = 1:nTraces
for i = 1:nTraces
    plot(1:size(cells,2),cells(i,:)+i*offset - offset)
    
end
    
% 3,5,6,10,
%%

% plot average evoked response for each ROI by frequency


figure;
for ii=1:length(uniquetone)
    gg=find(roiMatevokedstim(:,2)==uniquetone(ii));
    freqtempcells=roiMatevokedstim(:,gg,:,:);
    figure(ii)
    hold on;
     imagesc(squeeze(freqtempcells))
end

    %%

    for cc =  1:size(TC,1)
        for ii = 1:length(uniquetone)
            meantuningdf{cc}(ii,:) = squeeze(mean(roiMatevokedstim(cc,ii,:,responsewindow),4));
            meandurdf{cc}(ii,:)=squeeze(mean(roiMatevokeddurt(cc,ii,:,responsewindow),4));
            meanbaselinedf{cc}(ii,:) = squeeze(mean(roiMatevokedstim(cc,ii,:,5:10),4));
            mediantuningdf{cc}(ii,:) = squeeze(median(roiMatevokedstim(cc,ii,:,:),3));
            stdtuningdf{cc}(ii,:) = std(squeeze(roiMatevokedstim(cc,ii,:,responsewindow))');
        end
    end
    
    %%
    
    % Get 1) if cells is sound responsive 2) tuning tone 3) tuning duration
    % first do a Anova to test if the cell is responsive to any tonesbn
    % need 18 groups, each group with 8 data points for 8 trials
    % the baseline group is averaged from 8 trials
    
    %
    [~,signibaseline] = cellfun(@(x,y) ttest(reshape(x,size(x,1)*size(x,2),1),reshape(y,size(y,1)*size(y,2),1)),meanbaselinedf,meantuningdf);
    for cc =  1:size(TC,1)
        [Pval(cc,:),~,stats(cc)] = anova2(roiMatevokedtest{cc}',nbrep,'off'); % tones are columns and duration is repeated in 10 trials
    end
    %     signicelltuning = ;%cellfun(@(x) friedman(x',nbrep,'off'),meantuningdf);
    tunedcells = find(signibaseline<0.05);
    freqmaxresp = (cellfun(@(x) find(mean(x,2)==max(mean(x,2))),meantuningdf)); % only for single psth neuron purpose
    bftunedcellsind = (cellfun(@(x) find(mean(x,2)==max(mean(x,2))),meantuningdf(tunedcells)));
    bftunedcells = uniquetone(cellfun(@(x) find(mean(x,2)==max(mean(x,2))),meantuningdf(tunedcells)));
    
    %     colortuning = [(ones(1,10).*0.1:0.1:1)',(ones(1,10).*0.5)',flipud((ones(1,10).*0.1:0.1:1)')];%(cool(length(uniquetone)));
    durrespcells = find(Pval(:,2)<0.05);
    
    tuningrespcells = find(Pval(:,1)<0.05);%find(cellfun(@(x) x(1,2)<0.05,testtuningdur));
%     interactioncells = find(Pval(:,3)<0.05);%find(cellfun(@(x) x(1,3)<0.05,testtuningdur));
    
    %     tcnorm = cell2mat(cellfun(@(x) (mean(x,2)./max(mean(x,2))),meantuningdf(tuningrespcells),'UniformOutput',0));
    
    colormaptone = jet(length(uniquetone));
    for tid = 1:length(uniquetone)
        
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
        
        cmap{tid} = colorMap;%polarmap(256);%jet(256);
    end
    newtonecolormap = cell2mat(cellfun(@(x) x(end,:),cmap,'UniformOutput',0)');
    
    
    %%  normalize the data ?
    

    for cc =  1:size(TC,1)
        for ii = 1:length(uniquetone)
            [normroiMatevokedstim(cc,ii,:,responsewindow)] = math_scale_values( roiMatevokedstim(cc,ii,:,responsewindow), ...
                mean(squeeze(min(roiMatevokedstim(cc,ii,:,responsewindow)))), ...
                mean(squeeze(max(roiMatevokedstim(cc,ii,:,responsewindow)))), 0, 1.00);
            normmeantuningdf{cc}(ii,:) = squeeze(mean(roiMatevokedstim(cc,ii,:,responsewindow),4));
%             normmeanbaselinedf{cc}(ii,:) = squeeze(mean(roiMatevokedstim(cc,ii,:,5:10),4));
%             normmediantuningdf{cc}(ii,:) = squeeze(median(roiMatevokedstim(cc,ii,:,:),3));
%             normstdtuningdf{cc}(ii,:) = std(squeeze(roiMatevokedstim(cc,ii,:,responsewindow))');
        end
    end

        %% plot the tuned cells
        
        
figure(1);hold on; 
figure(2); hold on; 
            
for ccs = 1:length(tunedcells)
    if length(tunedcells)<40
            figP=figure(1);hold on;title('Mean tuning Cell #', num2str(tunedcells(ccs)));
            subplot(2,10,ccs);
            shadedErrorBar(1:length(uniquetone),mean(meantuningdf{ccs}'),std(meantuningdf{ccs}'),{'Color',[ 0.4843 0.7157 0.8882]});
            xticks([1:9])
            xticklabels(string(uniquetone)); 
    else
        if ccs <= length(tunedcells)/2
            figP=figure(1);hold on;title('Mean tuning Cell #', num2str(tunedcells(ccs)));
            subplot(4,10,ccs);
            shadedErrorBar(1:length(uniquetone),mean(normmeantuningdf{ccs}'),std(normmeantuningdf{ccs}'),{'Color',[ 0.4843 0.7157 0.8882]});
            xticks([1:9])
            xticklabels(string(uniquetone));
%         else
%             figure(2); hold on;title('Mean tuning', num2str(tunedcells(ccs)));
%             subplot(4,10,(ccs+1-ceil(length(tunedcells)/2)))
%             shadedErrorBar(1:length(uniquetone),normmean(meantuningdf{ccs}'),std(normmeantuningdf{ccs}'),{'Color',[ 0.4843 0.7157 0.8882]});
%                     xticks([1:9])
%             xticklabels(string(uniquetone));   
        end
    end
end
    %%
    % plot all cells signi tuning
    % first sort by bf
    colorredwhite = flipud([ones(100,1),(log10([1.01:0.1:11])/log10(11))',(log10([1.01:0.1:11])/log10(11))']);
    [dfmaxvalue,dfmaxind]=max(cell2mat(cellfun(@(x) mean(x'),meantuningdf(tuningrespcells),'UniformOutput',false)')');
    [sortbf,sortbfind] = sort(dfmaxind);
    h = figure(3); hold on; subplot(3,1,[1,2]); imagesc(1:length(uniquetone),1:length(tuningrespcells),cell2mat(cellfun(@(x) mean(x'),...
        meantuningdf(tuningrespcells(sortbfind)),'UniformOutput',0)'));...
%     xticklabels(mat2cell(uniquetone/1000,1,ones(1,10)));
    ylabel('significantly tuned ROIs'); title('frequency tuning');    colormap(colorredwhite);colorbar;
            xticklabels(string(uniquetone));  xlabel('frequency (Hz)')
    % tuning curve population for responsive cells
%     subplot(3,1,[3]);plot(1:length(uniquetone),cell2mat(cellfun(@(x) mean(x'),meantuningdf(tuningrespcells),'UniformOutput',0)')','k');
    % xticks([1:10]);
%     xticklabels(mat2cell(uniquetone/1000,1,ones(1,10)))
%     ylabel('df/f');
%     xlabel('frequency (Hz)')

            xticklabels(string(uniquetone));   
%%
            % same for duration tuning
    for tt = 1:length(uniqueduration)
        durationind(tt,:) = find(seqid==uniqueduration(tt));
        durationevoked(:,tt) = max(squeeze(median(roiMatevoked(:,durationind(tt,:),responsewindow),2))');
    end
    h = figure(6); hold on; subplot(3,1,[1,2]); imagesc(1:length(uniqueduration),1:length(durrespcells),...
        durationevoked); colormap cool(10); colorbar %caxis([1,1.1])
    xticks([1:length(uniqueduration)]);
    xticklabels(mat2cell(uniqueduration,1,ones(1,length(uniqueduration))))
    ylabel('Significantly tuned ROIs');
    subplot(3,1,[3]);plot(1:length(uniqueduration),durationevoked(durrespcells,:),'k');
    xticks([1:length(uniqueduration)]);
    xticklabels(mat2cell(uniqueduration,1,ones(1,length(uniqueduration))))
    ylabel('df/f');
    xlabel('Duration [ms]')
%     saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_DurationTuning.pdf'])
    %% now order the cells by max amplitude duration evoked
    [sortbd,maxdurind] = max(durationevoked');
    [rankedsortbd,rankedmaxdur]=sort(maxdurind);
    
    h = figure(7); hold on; subplot(3,1,[1,2]);
    imagesc(1:length(uniqueduration),1:length(durrespcells),...
        durationevoked(rankedmaxdur,:)); 
    colormap(colorredwhite);colorbar;
    xticks([1:length(uniqueduration)]);
    xticklabels(mat2cell(uniqueduration,1,ones(1,length(uniqueduration))))
    ylabel('significantly tuned ROIs');    title(' duration tuning');    xlabel('duration [ms]')
%     subplot(3,1,[3]);plot(1:length(uniqueduration),durationevoked((rankedmaxdur),:),'k');
%     xticks([1:length(uniqueduration)]);
%     xticklabels(mat2cell(uniqueduration,1,ones(1,length(uniqueduration))))
%     ylabel('df/f');
%     xlabel('duration [ms]')

    % seems to be some duration tuning but not a lot
    
    
    
    
%%

% relationship between freq, duration, loudness cells


% **!!! normalize tuning curves first by max and min evoked activity 
% (0 to 1) ** kindfo 


% 1) pick best 10 cells for frequency, look at duration tuning, do some
% preliminary analysis also on the loudness
% do a scatter plot of best freq by best duration
bdtunedcells=uniqueduration(cellfun(@(x) find(mean(x,2)==max(mean(x,2))),meandurdf(durrespcells)));

% find cells that are duration responsive and freq tuned
[val,pos]=intersect(durrespcells,tunedcells);
[val,pos2]=intersect(durrespcells,tunedcells);
figure;scatter(durrespcells(pos),bftunedcells(pos2))
xlabel('duration');
ylabel('frequency');
title('scatter plot of best frequency by best duration');

    %% mean tuning of cells

for ccs = 1:length(tunedcells)
    if length(tunedcells)<40
            figP=figure(1);hold on;title('Mean tuning Cell #', num2str(tunedcells(ccs)));
            subplot(2,10,ccs);
            shadedErrorBar(1:length(uniquetone),mean(meantuningdf{ccs}'),std(meantuningdf{ccs}'),{'Color',[ 0.4843 0.7157 0.8882]});
            xticks([1:9])
            xticklabels(string(uniquetone)); 
    else
        if ccs <= length(tunedcells)/2
            figP=figure(1);hold on;title('Mean tuning Cell #', num2str(tunedcells(ccs)));
            subplot(4,10,ccs);
            shadedErrorBar(1:length(uniquetone),mean(normmeantuningdf{ccs}'),std(normmeantuningdf{ccs}'),{'Color',[ 0.4843 0.7157 0.8882]});
            xticks([1:9])
            xticklabels(string(uniquetone));
%         else
%             figure(2); hold on;title('Mean tuning', num2str(tunedcells(ccs)));
%             subplot(4,10,(ccs+1-ceil(length(tunedcells)/2)))
%             shadedErrorBar(1:length(uniquetone),normmean(meantuningdf{ccs}'),std(normmeantuningdf{ccs}'),{'Color',[ 0.4843 0.7157 0.8882]});
%                     xticks([1:9])
%             xticklabels(string(uniquetone));   
        end
    end
end
%%


% how many neurons that have a best frequency have a duration tuning?

% k means using the mean tuning of freq, dur, and loudness

% if the k means doesn't work? 

% characteristic frequency = tuning cuyrve at the lowest level played
     % what frequency at the lowest dB is the cell responsive to?
     % the frequency even though when quiet has a reply from the cell
% intensity by frequency plots with spike rasters 
% https://www.jneurosci.org/content/39/35/6905 % heat maps of the response
% this plots intensity (loudness) by frequency with amplitude of repsonse


% https://www.jneurosci.org/content/32/18/6373.short
% duration tuning across vertebrates


% https://www.sciencedirect.com/science/article/pii/S0306452211008487?via%3Dihub#fig1

% isofrequency paper? 
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2440588/

%Neurons in CNIC exhibit sharp frequency tuning, have low thresholds and 
%are arranged tonotopically in isofrequency laminae with high sound
%frequencies represented ventromedially and low frequencies dorsolaterally
% (Willott and Urban, 1978; Stiebler and Ehret, 1985; Malmierca et al., 2008). 
% we see the opposite of this 

    %%
    % color site with BF
    nPlanes=1;tempRoi = stat(logical(iscell(:,1)));
    roisCoord = cell(1,nPlanes); % only 1 functional channel for suite2p
        for j = 1:nPlanes % number of planes
            roisCoord{1,j} = cell(1,length(tempRoi{j}));
            for k = 1:length(tempRoi) % neuron in each plane 
                bound = boundary(double(tempRoi{j,k}.xpix)', double(tempRoi{j,k}.ypix)',1); % restricted bound
                tempCoord = [tempRoi{j,k}.xpix(bound)' tempRoi{j,k}.ypix(bound)'];
                roisCoord{1,j}{k} = tempCoord;
            end
        end
    roisCoord=num2cell(roisCoord);
%     newtonecolormap=colormaptonotopy(9);
    %%
    col = [0.5,0.5,0.5];
    h=figure(4); 
        x = roisCoord{1}{1}{1}(:,1); %freehand rois have the outlines in x-y coordinates
    y = roisCoord{1}{1}{1}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
    hold on;
    imagesc(1:size(ops.meanImg,2),1:size(ops.meanImg,1),imadjust(int16(ops.meanImg))); colormap gray
%     cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',col),...
%         ((roisCoord{1}{1})))
%     cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',...
%         newtonecolormap(y,:)),roisCoord{1}{1}(tunedcells),...
%         mat2cell(bftunedcellsind,1,ones(1,length(bftunedcellsind))))

%     colorbar
%     patch(x,y,C(colormapIndex(bftunedcellsind,:)),'EdgeColor','none');
    xlim([1,size(ops.meanImg,2)]); ylim([1,size(ops.meanImg,1)])
%     set(gca,'XDir','reverse')
%         set(gca,'YDir','normal')
    
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    %%
%     saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_Tuningmap.pdf'])
    
%     if strcmp(animal, 'RoofBuddy2') || strcmp(animal, 'RoofBuddy1')
%         load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\z-stack\site' num2str(siteind(analysefilelistind(fn))) '_depth.mat'])
        %     load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\z-stack\coordpixsite.mat'])
        %     load([datadrivecolumn{analysefilelistind(fn)} ':\Jenni\' animal '\z-stack\site' num2str(siteind(analysefilelistind(fn))) '_vecstartpy.mat'])
        %
        for cc = 1:length(roisCoord{1})
            for px = 1:length(roisCoord{1}{cc})
                depthroi{cc}(px) = depth(floor(roisCoord{1}{cc}(px,2)/2),floor(roisCoord{1}{cc}(px,1)/2));
            end
        end
        mediandeppthrois = cellfun(@median,depthroi,'UniformOutput',1);
        % correct for depth 0
        mediandeppthrois(mediandeppthrois==0) = 1;

    
    %     roimicroncoord =
    if strcmp(animal, 'RoofBuddy2') || strcmp(animal, 'RoofBuddy1')
        
        colorredwhite = flipud([ones(200,1),(log10([1.01:0.05:11])/log10(11))',(log10([1.01:0.05:11])/log10(11))']);%([0:0.01:1-0.01].^1.2)',([0:0.01:1-0.01].^1.05)'];
        h=figure(); hold on; subplot(1,3,1)
        cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',newtonecolormap(y,:)),roisCoord{1}(tunedcells),mat2cell(bftunedcellsind,1,ones(1,length(bftunedcellsind))))
        xlim([1,size(ops.meanImg,2)]); ylim([1,size(ops.meanImg,1)])
        set(gca,'XDir','reverse')
        set(gca,'YDir','normal')
        
        set(gca,'xtick',[]); set(gca,'ytick',[]);
        hold on; subplot(1,3,2)
        cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',colorredwhite(round(y),:)),roisCoord{1},mat2cell(mediandeppthrois,1,ones(1,length(mediandeppthrois))))
        xlim([1,size(ops.meanImg,2)]); ylim([1,size(ops.meanImg,1)])
        set(gca,'XDir','reverse')
        set(gca,'YDir','normal')
        
        set(gca,'xtick',[]); set(gca,'ytick',[]);
        hold on; subplot(1,3,3)
        imagesc(1:size(ops.meanImg,2),1:size(ops.meanImg,1),imadjust(int16(ops.meanImg))); colormap gray
        set(gca,'XDir','reverse')
        set(gca,'YDir','normal')
        
        set(gca,'xtick',[]); set(gca,'ytick',[]);
        % depth scatter plot
        [depthvalue,depthind]=sort(mediandeppthrois(tunedcells));
        figure();
        plot(uniquetone(bftunedcellsind(depthind))/1000,depthvalue,'.');
        xticks(uniquetone/1000)
    end
    
    Tonotopystruct(siteind(analysefilelistind(fn))).cellstuning(1,:) = 1:length(freqmaxresp);
    Tonotopystruct(siteind(analysefilelistind(fn))).cellstuning(2,:) = zeros(1,length(freqmaxresp));
    Tonotopystruct(siteind(analysefilelistind(fn))).cellstuning(2,tunedcells) = 1;
    Tonotopystruct(siteind(analysefilelistind(fn))).cellstuning(3,:) =  zeros(1,length(freqmaxresp));
    Tonotopystruct(siteind(analysefilelistind(fn))).cellstuning(3,tunedcells) =  bftunedcellsind;
    Tonotopystruct(siteind(analysefilelistind(fn))).meantuningcurve =  meantuningdf;
    Tonotopystruct(siteind(analysefilelistind(fn))).mediantuningcurve =  mediantuningdf;
    Tonotopystruct(siteind(analysefilelistind(fn))).uniquetone =  uniquetone;
    Tonotopystruct(siteind(analysefilelistind(fn))).uniqueduration =  uniqueduration;
    Tonotopystruct(siteind(analysefilelistind(fn))).durationevoked =  durationevoked;
    Tonotopystruct(siteind(analysefilelistind(fn))).nbrep = nbrep;
    if strcmp(animal, 'RoofBuddy2') || strcmp(animal, 'RoofBuddy1')
        Tonotopystruct(siteind(analysefilelistind(fn))).depth = mediandeppthrois;
    end
    %     end
% end
% save(['N:\Jenni\' animal '\z-stack\Tonotopystruct.mat'],'Tonotopystruct')

% disp('')
% tuning duration % carefull with nb of repetition
%         signicelldurationtuning = cellfun(@(x) friedman(x',nbrep,'off'),meantuningdf);
%     tunedcells = find(signicelltuning<0.05);

%%

    for cc =  1:size(TC,1)
        for ii = 1:length(uniqueduration)
            for tid = 1:length(uniquetone)
                meanpsthdurationtuning{cc}{ii}(tid,:) = mean(squeeze(roiMatevoked(cc,intersect(find(seqid == uniqueduration(ii)),find(toneid == uniquetone(tid))),:)));
            end
            meantuningdfduration{cc}{ii} = squeeze(mean(roiMatevoked(cc,find(seqid == uniqueduration(ii)),responsewindow),3));
            plotmeantuningdur{cc}(ii,:) = mean((meanpsthdurationtuning{cc}{ii}(:,responsewindow)));
            plotmediantuningdur{cc}(ii,:) = median((meanpsthdurationtuning{cc}{ii}(:,responsewindow)));

        end
        meantuningdurationtuning{cc} = cell2mat(cellfun(@(x) mean(x,2),meantuningdfduration{cc},'UniformOutput',0));
    end

    [dfmaxvaluedur]= (cell2mat(cellfun(@(x) max(x),meantuningdurationtuning(tunedcells),'UniformOutput',false)')');
    dfmaxinddur = (cell2mat(cellfun(@(x) find(x==max(x)),meantuningdurationtuning(tunedcells),'UniformOutput',false)')');
    [sortdur,sortdurind] = sort(dfmaxinddur);
    bestdurationcellsind = cellfun(@(x) (find(x == max(x))),meantuningdurationtuning(tunedcells));




    h = figure(5); hold on; subplot(3,1,[1,2]); imagesc(uniqueduration,1:length(tunedcells),cell2mat(meantuningdurationtuning(tunedcells)')); colormap gray; %caxis([1,1.1]) % cell2mat(cellfun(@(x) mean(x),meantuningdurationtuning(tunedcells(sortdurind)),'UniformOutput',0)')
    %     xticklabels(mat2cell(uniqueduration,1,ones(1,length(uniqueduration))));
    ylabel('Significantly tuned ROIs');
%%
    subplot(3,1,[3]);plot(uniqueduration,cell2mat(meantuningdurationtuning(tunedcells)'),'k');
xticks([1:10]);
    xticklabels(mat2cell(uniquetone/1000,1,ones(1,10)))
    ylabel('df/f');
    xlabel('Tone duration [ms]')
    saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_DurationTuning.pdf'])


    figure(4);
    imagesc(1:10,1:5,meantuningdurationtuning{299}'); colormap gray;
    xticklabels(mat2cell(uniquetone/1000,1,ones(1,10)));
    yticks([1:5])
    yticklabels(mat2cell(uniqueduration,1,ones(1,5)));
    ylabel('durations [ms]')
    xlabel('Tone freq.')
    title('Example cell 299')
load(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\TonotopyMeanEvokedMovie.mat'])

single neuron pdf analysis output

%%

% load behavior file for upsweep and down sweep % Npo add indice of frames
% for the startstim
%     frameaddindice = nbframestmp(find(ismember(behaviorfileind{fn},sweepind(fn))));
%     load(['K:\Jenni\' animal '\behavior\' behaviorfile{sweepind(fn)} '.mat']); %RoofBuddy1tonotopy02-03-2021-12-50
%     sweepstartstim = batimagingdata.frametone+frameaddindice;%1:400:44000;
%
%
%     durid = batimagingdata.durationseq;%batimagingdata.;%repmat(1:11,1,10);
%     sweepid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
%     uniquesweep = unique(sweepid);
%     uniquedurationsweep = unique(durid); % also order of presentation
%     % nb of repetition
%     nbrep = (length(durid)/length(uniquesweep))/length(uniquedurationsweep);
%
%     for tid = 1:length(uniquesweep)
%         for dd = 1:length(uniquedurationsweep)
%             sweepordertone(tid,dd,:) = intersect(find(sweepid==tid),find(durid==uniquedurationsweep(dd))); % tone id duration and trial rep
%         end
%     end
%
%     for cc = 1:size(TC,1)
%         for ii = 1:length(sweepstartstim)
%             sweeproiMatevoked(cc,ii,:) = TC(cc,sweepstartstim(ii)-10:sweepstartstim(ii)+40)./median(TC(cc,sweepstartstim(ii)-10:sweepstartstim(ii)-1),2);
%         end
%     end
%     % cc is roi, ii is tone id, dr is duration, 4th dimension is trial nb, 5th dimension is
%     % time
%     sweeproiMatevokedstimdur = [];
%     sweeproiMatevokedstim = [];
%     for cc = 1:size(TC,1)
%         sweeproiMatevokedtest{cc} = [];
%
%         for ii = 1:length(uniquesweep)
%             for dr = 1:length(uniquedurationsweep)
%                 sweeproiMatevokedstimdur(cc,ii,dr,:,:) = squeeze(sweeproiMatevoked(cc,intersect(find(sweepid == uniquesweep(ii)),find(durid == uniquedurationsweep(dr))),:));%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
%             end
%             sweeproiMatevokedtest{cc}(ii,:) = mean(squeeze(sweeproiMatevoked(cc,reshape(sweepordertone(ii,:,:),1,length(uniquedurationsweep)*nbrep),responsewindow)),2)'; % % test repmat([1:5]',1,10)' according to this first 10
%             sweeproiMatevokedstim(cc,ii,:,:) = squeeze(sweeproiMatevoked(cc,find(sweepid == uniquesweep(ii)),:));
%         end
%     end
%
%     % test tuning WN etc, ...
%     %     sweeptesttuningdur = cellfun(@(x) anova2(x,nbrep,'off'),sweeproiMatevokedtest,'UniformOutput',0)
%     for cc =  1:size(TC,1)
%         [sweepPval(cc,:),~,sweepstats(cc)] = anova2(sweeproiMatevokedtest{cc}',nbrep,'off'); % tones are columns and duration is repeated in 10 trials
%     end
%     sweepdurrespcells = find(sweepPval(:,2)<0.05);
%     sweeptuningrespcells = find(sweepPval(:,1)<0.05);
%     sweepinteractioncells = find(sweepPval(:,3)<0.05);
%
%
%
%
%     % load echo delay
%     frameaddindice = nbframestmp(find(ismember(behaviorfileind{fn},echoind(fn))));
%     load(['N:\Jenni\' animal '\behavior\' behaviorfile{echoind(fn)} '.mat']); %RoofBuddy1tonotopy02-03-2021-12-50
%     echostartstim = batimagingdata.frametone+frameaddindice;%1:400:44000;
%
%
%     delayid = [10,12,14,16,18,20,2,4,6,8];%batimagingdata.;%repmat(1:11,1,10);
%     echoid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
%     uniqueecho = unique(echoid);
%     uniquedurationecho = unique(delayid); % also order of presentation
%     % nb of repetition
%
%     for cc = 1:size(TC,1)
%         for ii = 1:length(echostartstim)
%             echoroiMatevoked(cc,ii,:) = TC(cc,echostartstim(ii)-10:echostartstim(ii)+40)./median(TC(cc,echostartstim(ii)-10:echostartstim(ii)-1),2);
%         end
%     end
%     % cc is roi, ii is tone id, dr is duration, 4th dimension is trial nb, 5th dimension is
%     % time
%     for cc = 1:size(TC,1)
%         for ii = 1:length(uniqueecho)
%             echoroiMatevokedstim(cc,ii,:,:) = squeeze(echoroiMatevoked(cc,find(echoid == uniqueecho(ii)),:));%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
%         end
%     end
%
%
%     % do summary statistics % 1st output is colums duration % 2nd is rows
%     % tuning % 3rd is interaction
%     tonerespcells = numel(tunedcells)/length(signibaseline)*100;
%     prcttuningcells = numel(tuningrespcells)/length(signibaseline)*100;
%     prctdurationcells = numel(durrespcells)/length(signibaseline)*100;
%     prctofboth = 100*numel(intersect(tuningrespcells,durrespcells))/length(signibaseline);
%     prctofintraction = 100*(numel(interactioncells)/length(signibaseline));
%     prcttuningsweep = 100*numel(sweeptuningrespcells)/length(signibaseline);
%     prcttuningduration = 100*numel(sweepdurrespcells)/length(signibaseline);
%     %     labelsprcttuning = {''};
%     figure(50); hold on;
%     subplot(2,3,1); pie([tonerespcells,(100-tonerespcells)]);
%     title('% of tone responsive cells');
%     subplot(2,3,2); pie([prcttuningcells,(100-prcttuningcells)]);
%     title('% of frequency tuned cells');
%     subplot(2,3,3); pie([prctdurationcells,(100-prctdurationcells)]);
%     title('% of duration tuned cells');
%     subplot(2,3,4); pie([prctofintraction,(100-prctofintraction)]);
%     title('% of tone/duration tuned cells');
%     subplot(2,3,5); pie([prcttuningsweep,(100-prcttuningsweep)]);
%     title('% of complex sounds tuned cells');
%     subplot(2,3,6); pie([prcttuningduration,(100-prcttuningsweep)]);
%     title('% of duration for complex sounds tuned cells');
%     saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_SummaryStats.pdf'])
%
%
%     %     subplot(3,3,2); pie([prctdurationcells*100,(1-prctdurationcells)*100]);
%     % %     legend('Responsive','Non-responsive','Location','SouthWest'); legend box off
%     %     title('% of frequency duration cells');
%     %
%     %     % plot median and mean tuning pop
%     %         subplot(3,2,3); pl
%
%
%     % plot cell count per frequency
%
%     % plot cell count per duration
%
%
%     % mapping of duration tunin
%     colduration = autumn(5);
%     h=figure(65); hold on;
%     imagesc(1:size(ops.meanImg,2),1:size(ops.meanImg,1),imadjust(int16(ops.meanImg))); colormap gray
%     cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',col),roisCoord{1})
%     cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',colduration(y,:)),roisCoord{1}(tunedcells),mat2cell(bestdurationcellsind,1,ones(1,length(bestdurationcellsind))))
%     xlim([1,size(ops.meanImg,2)]); ylim([1,size(ops.meanImg,1)])
%     set(gca,'YDir','reverse')
%     set(gca,'xtick',[]); set(gca,'ytick',[]);
%     saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Population\' path{analysefilelistind(fn)}(end-4:end) '_DurationTuningmap.pdf'])
%     ylimval = [0.8,1.6];
%
%     uniquedurationechoandwn = unique([uniquedurationsweep,uniquedurationecho]);
%     colormapcomplex = cool(length(uniquedurationechoandwn));




%     load('K:\Jenni\RoofBuddy2\z-stack\site21\depthvalue_RoofBuddy2_013_011.mat')

%     for cc =[40,41,61,403,188,192,194,137]% 1:size(TC,1)
%         % plot example neurons for the poster
%         close all
%         h= figure(cc+6); hold on;
%         for tid = 1:length(uniquetone)
%             shadedErrorBar(((tid-1)*71)+1:((tid)*71),[nan(1,10),smooth(squeeze(mean(roiMatevokedstim(cc,tid,:,:),3)),3)',nan(1,10)],[nan(1,10),std(squeeze(roiMatevokedstim(cc,tid,:,:)))./sqrt(nbrep),nan(1,10)],{'Color',newtonecolormap(tid,:),'LineWidth',2},1);
%             plot([((tid-1)*71)+20,((tid-1)*71)+20],[0,max(max(squeeze(mean(roiMatevokedstim(cc,:,:,11:15),3))))],'Color',[0.8,0.8,0.8],'LineWidth',1.5)%;plot([((tid-1)*51),((tid-1)*51)],[0,2],'k','LineWidth',2)
%         end
%         ylim([0.95,2.1]); axis off;
%         set(gca,'fontsize',12); set(gca,'fontname','arial');
%         saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Examplecells\' path{analysefilelistind(fn)}(end-4:end) '_Roi' num2str(cc) '_MeanEvokedResponse.svg'])
%
%         h=figure(cc+7); hold on; colormap(gray);
%         imagesc(1:length(uniquetone),1:length(uniqueduration),cell2mat(cellfun(@(x) mean(x(:,responsewindow),2),meanpsthdurationtuning{cc},'UniformOutput',0))');axis tight; caxis([0.95,1.3]); %colorbar;
%         set(gca,'xtick',[1:length(uniquetone)]); set(gca,'xticklabels',cellfun(@(x) num2str(x),mat2cell(round(uniquetone/100)/10,1,ones(1,10)),'UniformOutput',0));
%         xlabel('Tone frequency [kHz]')
%         set(gca,'ytick',[1:length(uniqueduration)]); set(gca,'yticklabels',cellfun(@(x) num2str(x),mat2cell(uniqueduration,1,ones(1,length(uniqueduration))),'UniformOutput',0));
%         ylabel('Duration [ms]'); caxis([0.95,1.2]); colorbar;
%                set(gca,'fontsize',12); set(gca,'fontname','arial');
%         saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Examplecells\' path{analysefilelistind(fn)}(end-4:end) '_Roi' num2str(cc) '_ToneDurationMatrix.svg'])
%
%         h=figure(cc+8); hold on;  colormap(gray);
%         imagesc(1:length(uniquesweep)-1,1:length(uniquedurationsweep),squeeze(mean(mean(sweeproiMatevokedstimdur(cc,1:2,:,:,responsewindow),4),5))');axis tight; caxis([0.95,1.3]); %colorbar;
%         set(gca,'xtick',[1:length(uniquesweep)-1]); set(gca,'xticklabels',{'White Noise','Upsweep'});
%                 set(gca,'ytick',[1:length(uniquedurationsweep)]); set(gca,'yticklabels',cellfun(@(x) num2str(x),mat2cell(uniquedurationsweep,1,ones(1,length(uniquedurationsweep))),'UniformOutput',0));
%         set(gca,'fontsize',12); set(gca,'fontname','arial'); ylabel('Duration [ms]'); caxis([0.95,1.2]);
%         saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Examplecells\' path{analysefilelistind(fn)}(end-4:end) '_Roi' num2str(cc) '_ComplexDurationMatrix.svg'])
%
%         h=figure(cc+10); hold on;
%         sweepdurationcolor = winter(length(uniquedurationsweep));
%         for sid = 1:length(uniquesweep)-1
%             for durid = 1:length(uniquedurationsweep)
%                 clear meansweepcurve
%                 meansweepcurve = [nan(1,10),smooth(squeeze(mean(sweeproiMatevokedstimdur(cc,sid,durid,:,:),4))',3)',nan(1,10)];
%                 plot(((sid-1)*71)+1:((sid)*71),meansweepcurve(1,1:end),'Color',colormapcomplex(find(uniquedurationechoandwn == uniquedurationsweep(durid)),:),'LineWidth',2);
%             end
%             plot([((sid-1)*71)+20,((sid-1)*71)+20],[0,2],'Color',[0.8,0.8,0.8]);%plot([((sid-1)*51),((sid-1)*51)],[0,2],'k','LineWidth',2)
%         end
%         ylim([0.9,1.6]); axis off
%         legend(cellfun(@(x) num2str(x),mat2cell(uniquedurationsweep,1,ones(1,length(uniquedurationsweep))),'UniformOutput',0)); legend box off
%         set(gca,'fontsize',12); set(gca,'fontname','arial'); ylabel('Mean df/f'); caxis([0.95,1.2]);
%         saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Examplecells\' path{analysefilelistind(fn)}(end-4:end) '_Roi' num2str(cc) '_PSTHComplexDuration.svg'])
%
%
%         h= figure(cc+9); hold on; delaycolor=cool(length(uniquedurationecho));
%           sid = 1;
%         for durid = 1:length(uniquedurationecho)
%             plot(((sid-1)*51)+1:((sid)*51),smooth(squeeze(mean(echoroiMatevokedstim(cc,durid,:,:),3)),4),'Color',colormapcomplex(find(uniquedurationechoandwn == uniquedurationecho(durid)),:),'LineWidth',2)
%         end
%         legend(cellfun(@(x) num2str(x),mat2cell(uniquedurationecho,1,ones(1,length(uniquedurationecho))),'UniformOutput',0)); legend box off
%         plot([((sid-1)*51)+18,((sid-1)*51)+18],[0,2],'Color',[0.8,0.8,0.8]);%plot([((sid-1)*51),((sid-1)*51)],[0,2],'k','LineWidth',2)
%         ylim([0.9,1.6]); axis off
%         xlim([1,((sid)*51)]);
%                        set(gca,'fontsize',12); set(gca,'fontname','arial');
%
%            saveas(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Examplecells\' path{analysefilelistind(fn)}(end-4:end) '_Roi' num2str(cc) '_PSTHDelayTuning.svg'])
%
%
%     end
%
%         h=figure(6); %set(h, 'Visible', 'off');
% %         hold on; subplot(3,4,1); hold on;
%         imagesc(1:size(ops.meanImg,2),1:size(ops.meanImg,1),imadjust(int16(ops.meanImg))); colormap gray
%         xlim([1,size(ops.meanImg,2)]); ylim([1,size(ops.meanImg,1)])
%         set(gca,'YDir','reverse')
%         set(gca,'xtick',[]); set(gca,'ytick',[]);
%         %     cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',col),roisCoord{1}(cc))
%         cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',newtonecolormap(y,:)),roisCoord{1}(cc),mat2cell(freqmaxresp(cc),1,ones(1,length(freqmaxresp(cc))))) % modify to color cell with peak tuning (even if not signi)
%         % plot mean psth per tone
%         subplot(3,4,[2,3,4]); hold on;
%         for tid = 1:length(uniquetone)
%             shadedErrorBar(((tid-1)*51)+1:((tid)*51),squeeze(mean(roiMatevokedstim(cc,tid,:,:),3)),std(squeeze(roiMatevokedstim(cc,tid,:,:))),{'Color',newtonecolormap(tid,:)},1);
%             plot([((tid-1)*51)+10,((tid-1)*51)+10],[0,2],'k--');plot([((tid-1)*51),((tid-1)*51)],[0,2],'k','LineWidth',2)
%         end
%         ylim(ylimval); ylabel('mean df/f'); xlabel(['PSTH per tone' ' (' num2str(xvalues(1)) '-'  num2str(xvalues(end)) 's)'])
%         set(gca,'xtick',[((0:9)*51)+10]); set(gca,'xticklabels',cellfun(@(x) num2str(x),mat2cell([zeros(1,10)],1,ones(1,10)),'UniformOutput',0));
%         xlim([1,(tid)*51]);
%         % plot mean tuning curve / plot median tuning curve window (peak amp for 5 first frames)
%         subplot(3,4,[5]); hold on;
%         plot(uniquetone,max(squeeze(mean(roiMatevokedstim(cc,:,:,responsewindow),3))'),'k','Linewidth',2);
%         plot(uniquetone,max(squeeze(median(roiMatevokedstim(cc,:,:,responsewindow),3))'),'Color',[0.5,0.5,0.5],'Linewidth',2);
%         ylim(ylimval); ylabel('Peak amp over 5 first frames');xlim([min(uniquetone),max(uniquetone)]);
%         set(gca,'xtick',uniquetone); set(gca,'xticklabels',cellfun(@(x) num2str(x),mat2cell(round(uniquetone/100)/10,1,ones(1,10)),'UniformOutput',0));
%         xlabel('Tone frequency [kHz]'); legend('Mean','Median'); legend box off;
%
%         subplot(3,4,[6,7]); hold on; % that's not the correct matrice for this
%         imagesc(1:length(uniquetone),1:length(uniqueduration),cell2mat(cellfun(@(x) mean(x(:,responsewindow),2),meanpsthdurationtuning{cc},'UniformOutput',0))');axis tight; caxis([0.95,1.3]); %colorbar;
%         set(gca,'xtick',[1:length(uniquetone)]); set(gca,'xticklabels',cellfun(@(x) num2str(x),mat2cell(round(uniquetone/100)/10,1,ones(1,10)),'UniformOutput',0));
%         xlabel('Tone frequency [kHz]')
%         set(gca,'ytick',[1:length(uniqueduration)]); set(gca,'yticklabels',cellfun(@(x) num2str(x),mat2cell(uniqueduration,1,ones(1,length(uniqueduration))),'UniformOutput',0));
%         ylabel('Duration [ms]')
%
%         subplot(3,4,[8]); hold on;
%         plot(max(plotmeantuningdur{cc}'),uniqueduration,'k','Linewidth',2);
%         plot(max(plotmediantuningdur{cc}'),uniqueduration,'Color',[0.5,0.5,0.5],'Linewidth',2);
%         ylim([min(uniqueduration),max(uniqueduration)]); xlabel('Peak amp over 5 first frames');xlim(ylimval);
%         set(gca,'ytick',uniqueduration); set(gca,'yticklabels',cellfun(@(x) num2str(x),mat2cell(uniqueduration,1,ones(1,length(uniqueduration))),'UniformOutput',0));
%         ylabel('Duration [ms]'); legend('Mean','Median'); legend box off;
%
%         % plot cells response as function of sweep id per time
%         subplot(3,4,[9]); hold on;
%         for sid = 1:length(uniquesweep)
%             for durid = 1:length(uniquedurationsweep)
%                 clear meansweepcurve
%                 meansweepcurve = squeeze(mean(sweeproiMatevokedstimdur(cc,sid,durid,:,:),4))';
%                 plot(((sid-1)*51)+1:((sid)*51),meansweepcurve(1,1:end),'Color',[0.8,0.8,0.8]-0.1*(durid-1));
%             end
%             plot([((sid-1)*51)+10,((sid-1)*51)+10],[0,2],'k--');plot([((sid-1)*51),((sid-1)*51)],[0,2],'k','LineWidth',2)
%         end
%         ylim(ylimval); ylabel('mean df/f'); xlabel(['PSTH per sweep' ' (' num2str(xvalues(1)) '-'  num2str(xvalues(end)) 's)'])
%         set(gca,'xtick',[((0:2)*50)+10]); set(gca,'xticklabels',{'WN','US','DS'});
%
%         subplot(3,4,[10]); hold on;
%         imagesc(1:length(uniquesweep),1:length(uniquedurationsweep),squeeze(mean(mean(sweeproiMatevokedstimdur(cc,:,:,:,responsewindow),4),5))');axis tight; caxis([0.95,1.3]); %colorbar;
%         set(gca,'xtick',[1:length(uniquesweep)]); set(gca,'xticklabels',{'WN','US','DS'});
%         xlabel('Stimuli type')
%         set(gca,'ytick',[1:length(uniquedurationsweep)]); set(gca,'yticklabels',cellfun(@(x) num2str(x),mat2cell(uniquedurationsweep,1,ones(1,length(uniquedurationsweep))),'UniformOutput',0));
%         ylabel('Duration [ms]')
%
%
%         % Plot echo delay tuning
%         subplot(3,4,[11]);  hold on;
%         sid = 1;
%         for durid = 1:length(uniquedurationecho)
%             plot(((sid-1)*51)+1:((sid)*51),squeeze(mean(echoroiMatevokedstim(cc,durid,:,:),3)),'Color',[0.8,0.8,0.8]-0.06*(durid-1))
%         end
%         plot([((sid-1)*51)+10,((sid-1)*51)+10],[0,2],'k--');%plot([((sid-1)*51),((sid-1)*51)],[0,2],'k','LineWidth',2)
%         ylim(ylimval); ylabel('mean df/f'); xlabel(['PSTH per echo delay' ' (' num2str(xvalues(1)) '-'  num2str(xvalues(end)) 's)'])
%         set(gca,'xtick',[10]); set(gca,'xticklabels',{'HypSweep+Echo'});xlim([1,((sid)*51)]);
%
%         subplot(3,4,[12]); hold on;
%         plot(uniqueduration,max(squeeze(mean(echoroiMatevokedstim(cc,:,:,responsewindow+7),3))),'k','LineWidth',2)
%         plot(uniqueduration,max(squeeze(median(echoroiMatevokedstim(cc,:,:,responsewindow+7),3))),'Color',[0.5,0.5,0.5],'LineWidth',2)
%         xlabel('Delay [ms]'); legend('Mean','Median'); legend box off; xlim([uniquedurationecho(1),uniquedurationecho(end)])
%
%
%         suptitle(['ROI # ' num2str(cc)]);
%         set(h, 'Visible', 'off');
%         print(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Neurons\ROI_' num2str(cc)],'-dpdf')
% %         print(h,['N:\Jenni\' animal '\figure\tonotopy\site' num2str(siteind(analysefilelistind(fn))) '\Neurons\ROI_' num2str(cc)],'-dsvg')
%
%         %    saveas(h,['N:\Jenni\' animal '\figure\site1\' path{analysefilelistind(fn)}(end-4:end) '_ROI_' num2str(cc) '.pdf'])
%
%
%     end


% end


