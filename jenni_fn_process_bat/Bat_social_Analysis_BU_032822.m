function Bat_social_Analysis(analysis)
animal = 'Gray2';
animaldrive = 'K';
datatable = readtable([ animaldrive ':\Jenni\' animal '\ProcessingProgress.xlsx']);

manualroisdone = find(table2array(datatable(:,13))==1);%find(~cellfun(@isempty,cellfun(@(x) find(x=='1'),table2array(datatable(:,13)),'UniformOutput',0)));%find(table2array(datatable(:,13))==1);
daytraining = table2array(datatable(:,1));
siteind = table2array(datatable(:,5));
preprocessedrois = unique(siteind(manualroisdone));
framerate = table2array(datatable(:,8));
protocol = table2array(datatable(:,7));
path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,10));
pclampfile = table2array(datatable(:,11));
stimselection = find(table2array(datatable(:,14))==1);
datadrivecolumn = (table2array(datatable(:,15)));

tonotopyfilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'tonotopy')>0,behaviorfile,'UniformOutput',0))));
sweepfilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'LinearSweep')>0,behaviorfile,'UniformOutput',0))));
angiefilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Angie')>0,behaviorfile,'UniformOutput',0))));
echofilesind = find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Echo')>0,behaviorfile,'UniformOutput',0))));
melfilesind =  find(~cellfun(@isempty,(cellfun(@(x) strfind(x,'Mel')>0,behaviorfile,'UniformOutput',0))));

lowframeind = find(framerate==31.25);
%
% load('N:\Jenni\RoofBuddy2\figure\correctangiestim.mat')

for ii = 1:length(preprocessedrois)
    behaviorfileind{ii} = intersect(find(siteind==preprocessedrois(ii)),lowframeind);
    analysefilelistind(ii) = find(siteind==preprocessedrois(ii),1);
    if length(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),tonotopyfilesind))>1 % for RoofBuddy1 1st recording linear sweep was called tonotopy
        tonotopyind(ii) = max(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),tonotopyfilesind));%find(daytraining==preprocessedrois(ii),1),find(behaviorfile)); % change this to compare the name of the file and 'Mel' strcmp()
    else
        tonotopyind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),tonotopyfilesind);%find(daytraining==preprocessedrois(ii),1),find(behaviorfile)); % change this to compare the name of the file and 'Mel' strcmp()
    end
    if length(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),tonotopyfilesind))>1
        sweepind(ii) = min(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),tonotopyfilesind));
    else
        sweepind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),sweepfilesind);
    end
    %     if numel(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind))==0
    %         angieind(ii) = nan;
    %     else
    if length(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),angiefilesind))>1
        angieind(ii) = max(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),angiefilesind));
    else
        angieind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),angiefilesind);
    end
    %     end
    echoind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),echofilesind);
    if numel(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind))>1 % case of when there was 2 mel files
        melind(ii) = max(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind));
    elseif numel(intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind))==0
        melind(ii) = nan;
    else
        melind(ii) = intersect(intersect(lowframeind,find(siteind==preprocessedrois(ii))),melfilesind);
    end
end

if ~exist([animaldrive ':/Jenni/' animal '/roistructure/socialroisdata.mat'])
    for fn = 1:length(analysefilelistind)
        close all
        clearvars -except animal manualroisdone datatable angieind fn analysefilelistind behavindnb behaviorfileind behaviorfile lowframeind preprocessedrois ...
            socialroisdata pclampfile datadrivecolumn stimselection path framerate protocol siteind angiefilesind DataDrive
        
        DataDrive = datadrivecolumn{analysefilelistind(fn)}; % data drive for data, excel file is in N and stim too
        
        clear TC roiMatevoked tcnorm tunedcells signicelltuning meantuningdf roiMatevokedstim centroidcells dist
        % read frames from h5
        h5list = dir([DataDrive ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' '*.h5']);
        nbframestmp(1) = 0;
        for filenb = 1:length(h5list)
            clear hinfo
            hinfo = hdf5info([DataDrive ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' h5list(filenb).name]);
            nbframestmp(filenb+1) = hinfo.GroupHierarchy.Datasets.Dims(3);
        end
        
        % load behavior file
        load([DataDrive ':\Jenni\' animal '\behavior\' behaviorfile{angieind(fn)} '.mat']); %RoofBuddy1tonotopy02-03-201-12-50
        % get the nb of frames for the h5 placement
        behavindnb = (angieind(fn)-analysefilelistind(fn))+1;
        cumsumframes = cumsum(nbframestmp);
        startstim = batimagingdata.frametone+cumsumframes(behavindnb);
        if strcmp(animal,'RoofBuddy1') && siteind(analysefilelistind(fn))==9
            startstim = startstim-5;
        elseif strcmp(animal,'RoofBuddy1') && siteind(analysefilelistind(fn))==8
            startstim = startstim-2;
        end
        toneid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
        uniquetone = unique(toneid);
        for tid = 1:length(uniquetone)
            ordertone(tid,:) = find(toneid==uniquetone(tid)); % tone id duration and trial rep
        end
        % nb of repetition
        nbrep = (length(toneid)/length(uniquetone));
        
        % load ROIs
        load([DataDrive ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_TC_plane0.mat']);
        load([DataDrive ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
        load([DataDrive ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\31.25_rois_coord_plane0.mat']);
        % load mean image
        load([DataDrive ':\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\suite2p\plane0\fall.mat'])
        
        
        % get abf file for recording session
        %         [abfdata,samplingrate,header] = abfload([DataDrive ':\Jenni\' animal '\pclamp\',[pclampfile{angieind(fn)}, '.abf']]);
        %         Fs = 1/(samplingrate*(10^-6)); % capture the sampling rate in Hz
        %
        %         % Identify channel with all potential channel types
        %         photogate = find(strcmp(header.recChNames, 'Photogate'));
        %         sound_channel = find(strcmp(header.recChNames, 'Sound'));
        %         framechannel = find(strcmp(header.recChNames, 'Frame'));
        %
        %         % Identify the samples for each frame, frame data is always channel 1
        %         frametransitions = diff(abfdata(:,framechannel));%abs(diff(abfdata(:,framechannel))); %calculate derivative of frame signal to identify transitions
        %         ind = find(frametransitions>0.4); %identify when there is a positive change in the signal (frame start)
        %         ind2 = diff(ind); %the frame signal transition sometimes takes more than one sample, so need to identify the big transition points
        %         ind3 = find(ind2>80); %length of ind3 should now be the total number of frames acquired
        %         numsamples_perframe = round(mean(ind2(ind3))); %average number of samples per frame, e.g. 25600 for 100KhZ at 3.91 frames per second
        %
        %         ind_end = ind(ind3); %now just identify the exact samples where frame transition occurs [[this will be "last sample in frame n"]]
        %         ind_end = [ind_end;ind_end(end)+numsamples_perframe];
        %         ind_start = [ind_end(1:(length(ind_end)-1))]; %create index of frame start
        %
        %         numframes = length(ind_end); % tot nb of frames +1
        %
        %         filectr = 1;
        %         frame_recording = abfdata(:,framechannel);
        %         %     lickingchannel = find(strcmp(header.recChNames, 'Dig_licks'));
        %         %     rewardchannel = find(strcmp(header.recChNames, 'Water'));
        %         mic_recording{filectr} = abfdata(:,photogate);
        %         sound_recording{filectr} = abfdata(:,sound_channel);
        %         % load corresponding behavrio file and get when the stim where played +
        %         % length stim duration
        %
        %         stimframeind = batimagingdata.frametone;
        %         lengthstimincalamp = (unique(batimagingdata.durationseq)/1000+0.2)*Fs;% add 200 ms of buffe 100ms before and 100 ms after
        %         correspondingframe = ind_start(stimframeind);
        %
        stimwavname = dir(['N:\Jenni\wavfiles\Angie\' '*.wav']);
        for si = 1:length({stimwavname.name})
            wav_tmp{si} = audioread(['N:\Jenni\wavfiles\Angie\' stimwavname(si).name]);
        end
        %
        if ~isempty(intersect(stimselection,analysefilelistind(fn)))
            correctedstimseq = repmat([2:17,1],1,10);
            correctedstim = [2:6,8:12,15,17];
        else
            correctedstimseq = repmat([1:17],1,10);
            correctedstim = [1:17];
        end
        
        %
        %     clear Matmicstim Matsoundstim
        %         for si = 1:length(stimframeind)
        %                 Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):((correspondingframe(si))+(lengthstimincalamp)));
        %                 Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+(lengthstimincalamp)));
        %         end
        %
        %         for si = 1:length({stimwavname.name})
        %             clear h
        %             h =figure(si); hold on;
        %                 subplot(1,3,1);spectrogram(wav_tmp{si},400,300,[],200000,'yaxis'); % plot theorithecal wav file
        %                 subplot(1,3,2);spectrogram(reshape(Matsoundstim(1:round((length(wav_tmp{si})/200000)*250000),find(correctedstimseq==si)),[((length(wav_tmp{si})/200000)*250000)*10,1]),400,300,[],250000,'yaxis'); ylim([0 100]); %(0.1*250000):(0.15*250000),si)),400,300,[],250000,'yaxis'); ylim([0 100]); %ylim([0,100000]); % plot sound input
        %                 subplot(1,3,3);spectrogram(mean(Matsoundstim(1:round((length(wav_tmp{si})/200000)*250000),find(correctedstimseq==si)),2),400,300,[],250000,'yaxis'); ylim([0 100]); %(0.1*250000):(0.15*250000),si)),400,300,[],250000,'yaxis'); ylim([0 100]); %ylim([0,100000]); % plot sound input
        %                 saveas(h,['N:\Jenni\' animal '\figure\social\Stim\site' num2str(siteind(analysefilelistind(fn))) '\' path{analysefilelistind(fn)}(end-4:end) '_' [stimwavname(si).name] '.pdf'])
        %         end
        
        %TC = TC;%./median(TC')';
        
        responsewindow = 10:18;
        baselinewindow = 1:9;
        
        
        for cc = 1:size(TC,1)
            for ii = 1:length(startstim)
                roiMatevoked{cc}(ii,:) = TC(cc,startstim(ii)-9:startstim(ii)+41);%./median(TC(cc,startstim(ii)-10:startstim(ii)),2);
                % roiMatevokednorm(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40)./median(TC(cc,startstim(ii)-10:startstim(ii)),2);
            end
        end
        %         temporal shiuffle
        %         for cc = 1:size(roiMatevoked,1)
        %             for trialnb = 1:size(roiMatevoked,2)
        %                 for ii =1:10
        %                     shuffledmat(cc,trialnb,:,ii) = roiMatevoked(cc,trialnb,randperm(51));
        %                 end
        %             end
        %         end
        
        % cc is roi, ii is tone id, 3rd dimension is trial nb, 4th dimension is
        % time
        for cc = 1:size(TC,1)
            %         cellfun(@(x) x(find(correctedstimseq==1),:),roiMatevoked,'UniformOutput',0)
            for ii = 1:length(correctedstim)
                %                 bla(ii,:,:) = cellfun(@(x) x(find(correctedstimseq==ii),:),roiMatevoked,'UniformOutput',0);
                roiMatevokedstim(cc,ii,:,:) = squeeze(roiMatevoked{cc}(find(correctedstimseq == correctedstim(ii)),:))./median(squeeze(roiMatevoked{cc}(find(correctedstimseq == correctedstim(ii)),1:9)),2);%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
                clear correctedstimseqshuffle
                correctedstimseqshuffle = (cellfun(@(x) datasample(x,length(x),'Replace',false),repmat({correctedstimseq},100,1),'UniformOutput',0));
                % %                     [stimindx,stimindy] = find(correctedstimseqshuffle == correctedstim(ii));
                roiMatevokedstimshuffletmp{cc}(ii,:,:,:) = reshape(cell2mat(cellfun(@(x,y) x(find(y==ii),:),repmat({squeeze(roiMatevoked{cc}(:,:))},100,1),correctedstimseqshuffle,'UniformOutput',0)),10,100,51);
                %                     roiMatevokedstimshuffle(cc,ii,:,:,:) = squeeze(roiMatevoked(cc,correctedstimseqshuffle(stimindx,stimindy),:));%squeeze(shuffledmat(cc,find(correctedstimseq == correctedstim(ii)),:,:));%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
            end
        end
        nbcells = size(roiMatevokedstim,1);
        nbstim = size(roiMatevokedstim,2);
        nbrep = size(roiMatevokedstim,3);
        
        roiMatevokedstimshuffle = reshape(cat(1,roiMatevokedstimshuffletmp{:}),nbcells,nbstim,nbrep,100,size(roiMatevokedstim,4));
        roiMatevokedstimshuffle = permute(roiMatevokedstimshuffle,[1,2,3,5,4]);
        %             bli = permute(squeeze(bla),[2,1]);
        
        % analysis per category
        
        matstimevoked = roiMatevokedstim(:,:,:,responsewindow);
        matstimevokedshuffle = roiMatevokedstimshuffle(:,:,:,responsewindow,:);
        
        % ttest per stim
        % Do paired tests per stim per cells
        PeakAmp = squeeze(max(permute(matstimevoked(:,:,:,:),[4,1,2,3])));
        MedianAmp = squeeze(median(permute(matstimevoked,[4,1,2,3])));
        MedianAmpshuffle = squeeze(median(permute(matstimevokedshuffle,[4,1,2,3,5])));
        AreaAmp = squeeze(trapz(permute(matstimevoked,[4,1,2,3])))./length(responsewindow);
        
        % test for category analysis
        if nbstim == 12
            Category = [ones(1,5),ones(1,5)*2,ones(1,2)*3];
        elseif nbstim == 17
            Category = [ones(1,6),ones(1,6)*2,4,ones(1,4)*3];
        end
        % do anova per cells over stims + baseline
        for cc = 1:nbcells
            % category test 1 vs 2 (FMBs vs Echo)
            peaktestcat{cc} = [reshape(PeakAmp(cc,find(Category==1),:),nbrep*max(find(Category==1)),1),reshape(PeakAmp(cc,find(Category==2),:),nbrep*max(find(Category==1)),1)];
            medtestcat{cc} = [reshape(MedianAmp(cc,find(Category==1),:),nbrep*max(find(Category==1)),1),reshape(MedianAmp(cc,find(Category==2),:),nbrep*max(find(Category==1)),1)];
            areatestcat{cc} = [reshape(AreaAmp(cc,find(Category==1),:),nbrep*max(find(Category==1)),1),reshape(AreaAmp(cc,find(Category==2),:),nbrep*max(find(Category==1)),1)];
        end
        % identify cvells that are category specific with anova1
        peakanovacells = find(cellfun(@(x) anova1(x,[],'off'),peaktestcat)<=0.05);
        medanovacells = find(cellfun(@(x) anova1(x,[],'off'),medtestcat)<=0.05);
        areaanovacells = find(cellfun(@(x) anova1(x,[],'off'),areatestcat)<=0.05);
        % compute selectivity index per cell
        tuningpeak = mean(PeakAmp,3);
        tuningmed = mean(MedianAmp,3);
        tuningarea = mean(AreaAmp,3);
        dfmed = tuningmed-1; % in the sojner paper they remove the baseline, 1 is an approximation
        dfcat1med = mat2cell(dfmed(:,find(Category==1)),ones(1,nbcells),max(find(Category==1)));
        dfcat2med = mat2cell(dfmed(:,find(Category==2)),ones(1,nbcells),max(find(Category==1)));
        dfpeak = tuningpeak-1; % in the sojner paper they remove the baseline, 1 is an approximation
        dfcat1peak = mat2cell(dfpeak(:,find(Category==1)),ones(1,nbcells),max(find(Category==1)));
        dfcat2peak = mat2cell(dfpeak(:,find(Category==2)),ones(1,nbcells),max(find(Category==1)));
        dfarea = tuningarea-1; % in the sojner paper they remove the baseline, 1 is an approximation
        dfcat1area = mat2cell(dfarea(:,find(Category==1)),ones(1,nbcells),max(find(Category==1)));
        dfcat2area = mat2cell(dfarea(:,find(Category==2)),ones(1,nbcells),max(find(Category==1)));
        %             for stim = 1:max(max(find(Category==1)))
        selectivityindexperstimmed = cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),dfcat1med,dfcat2med,'UniformOutput',0);
        selectivityindexperstimpeak = cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),dfcat1peak,dfcat2peak,'UniformOutput',0);
        selectivityindexperstimarea = cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),dfcat1area,dfcat2area,'UniformOutput',0);
        
        %shuffle data mat for selectivity significance
        tuningmedshuffle = squeeze(mean(MedianAmpshuffle,3));
        dfmedshuffle = tuningmedshuffle-1; % in the sojner paper they remove the baseline, 1 is an approximation
        dfcat1medshuffle = mat2cell(dfmedshuffle(:,find(Category==1),:),ones(1,nbcells));%find(Category==1)),ones(1,nbcells),max(find(Category==1)),10);
        dfcat2medshuffle = mat2cell(dfmedshuffle(:,find(Category==2),:),ones(1,nbcells));%,ones(1,nbcells),max(find(Category==1)))
        selectivityindexperstimmedshuffle = cellfun(@(x,y) (abs(squeeze(x))-abs(squeeze(y)))./(abs(squeeze(x))+abs(squeeze(y))),dfcat1medshuffle,dfcat2medshuffle,'UniformOutput',0);
        
        % save ROIs mat for further analysis
        socialroisdata(fn).site = siteind(analysefilelistind(fn));
        socialroisdata(fn).roiMatevoked = roiMatevoked;
        socialroisdata(fn).roiMatevokedstim = roiMatevokedstim;
        %         socialroisdata(fn).roiMatevokedstimshuffle = roiMatevokedstimshuffle;
        socialroisdata(fn).correctedstim = correctedstim;
        socialroisdata(fn).correctedstimseq = correctedstimseq;
        socialroisdata(fn).framerate = framerate(analysefilelistind(fn));
        socialroisdata(fn).roisCoord = roisCoord;
        socialroisdata(fn).meanImg = ops.meanImg;
        socialroisdata(fn).selectivityindexperstimmed = selectivityindexperstimmed;
        socialroisdata(fn).selectivityindexperstimmedshuffle = selectivityindexperstimmedshuffle;
        socialroisdata(fn).selectivityindexperstimpeak = selectivityindexperstimpeak;
        socialroisdata(fn).selectivityindexperstimarea = selectivityindexperstimarea;
        socialroisdata(fn).peakanovacells = peakanovacells;
        socialroisdata(fn).medanovacells = medanovacells;
        socialroisdata(fn).areaanovacells = areaanovacells;
        % compute in microns for that 2X and 4X are different
        
        if strcmp(animal,'RoofBuddy1')
            centroidcells = cellfun(@(x) median((x./[2.55,2.66])),roisCoord{1},'UniformOutput',0);
        elseif strcmp(animal,'RoofBuddy2') && siteind(analysefilelistind(fn))==1 || siteind(analysefilelistind(fn))==2
            centroidcells = cellfun(@(x) median((x./[2.55,2.66])),roisCoord{1},'UniformOutput',0);
        else
            centroidcells = cellfun(@(x) median((x./[1.32,1.37])),roisCoord{1},'UniformOutput',0);
        end
        
        socialroisdata(fn).centroidcells =centroidcells;
        
        for cc = 1:length(centroidcells)
            for cci = length(centroidcells):-1:cc
                
                dist(cc,cci) = pdist([centroidcells{cc};centroidcells{cci}],'euclidean');
            end
        end
        socialroisdata(fn).dist = dist;
        
    end
%     if strcmp(animal,'RoofBuddy2')
%     save([DataDrive ':/Jenni/' animal '/roistructure/socialroisdata.mat'],'socialroisdata','-v7.3') % not for RoofBuddy 2 (split in between 2 drives)
%     else
    save([DataDrive ':/Jenni/' animal '/roistructure/socialroisdata.mat'],'socialroisdata','-v7.3') % not for RoofBuddy 2 (split in between 2 drives)
%     end
end

% load roi data structure
DataDrive = datadrivecolumn{analysefilelistind(1)}; % data drive for data, excel file is in N and stim too
datarois = load([DataDrive ':/Jenni/' animal '/roistructure/socialroisdata.mat']);
datarois = datarois.socialroisdata;

% colormap stim FMBs red, echo blue (because of the selectivity index,
% social towars 1, so red), mouse voc if not fucked up green
stimcolormap = zeros(17,3);
for stimind = 1:6
    stimcolormap(stimind,:) = [1-(0.1*(stimind)),0.3,0.2]; % red for FMBs 
end
for stimind = 7:12
    stimcolormap(stimind,:) = [0.2,0.3,1-(0.1*(stimind-6))]; % Blue for echolocation 
end
stimcolormap(13,:) = [0,0,0];
for stimind = 14:17
    stimcolormap(stimind,:) = [0.2,1-(stimind-13)/10,0.2];
end
xvalues = [-10:40]/31.25;
stimwavname = dir(['N:\Jenni\wavfiles\Angie\' '*.wav']);
for si = 1:length({stimwavname.name})
    wav_tmp{si} = audioread(['N:\Jenni\wavfiles\Angie\' stimwavname(si).name]);
end
Lengthstim = max(cellfun(@(x) size(x,1)/200000,wav_tmp)); % stim are 100ms apart fron the mouse calls
colorredwhite = flipud([ones(100,1),(log10([1.01:0.1:11])/log10(11))',(log10([1.01:0.1:11])/log10(11))']);%([0:0.01:1-0.01].^1.2)',([0:0.01:1-0.01].^1.05)'];
Socialind = find(cellfun(@(x) sum(strfind(x,'FMB')), {stimwavname.name})>0);
Echoind = find(cellfun(@(x) sum(strfind(x,'Echo')), {stimwavname.name})>0);
Noise = find(cellfun(@(x) sum(strfind(x,'noise')), {stimwavname.name})>0);
Mouseind = find(cellfun(@(x) sum(strfind(x,'voc')), {stimwavname.name})>0);

for fn = 1:length(analysefilelistind)
    switch analysis
        case 'population'
            if strcmp(animal,'RoofBuddy1')
                datarois = datarois(2:end);
            end
            % get the stim to compare % also adapt the SocialInd and Echo
            % Ind to the corresponding corrected stim sequences
            responsewindow = 11:19; % from sound onset
            baselinewindow = 1:10; % 10 frames before the sound
            if size(datarois(1).correctedstim,2) == 17
            else
                %             concatmatevokedFMB = cell2mat(cellfun(@(x) x(:,1:5,:,:),{datarois(:).roiMatevokedstim},'UniformOutput',0)');
                %             indEcho = cellfun(@(x) size(x,2) ,{datarois(:).roiMatevokedstim},'UniformOutput',1);
                %             indEcho(indEcho==12)= 6;
                %             indEcho(indEcho==17)= 7;
                %             concatmatevokedEcho = cell2mat(cellfun(@(x,y) x(:,y:y+4,:,:),{datarois(:).roiMatevokedstim},num2cell(indEcho),'UniformOutput',0)'); % select right category for corrected stim
                Socialind = 1:5; % double check this
                Echoind = 6:10; % double check this
                Noise = [];
                Mouseind = 15:16;
            end
            concatmatevokedFMB = cell2mat(cellfun(@(x) x(:,Socialind,:,:),{datarois(:).roiMatevokedstim},'UniformOutput',0)');
            concatmatevokedEcho = cell2mat(cellfun(@(x) x(:,Echoind,:,:),{datarois(:).roiMatevokedstim},'UniformOutput',0)');
            nbcellssite = cellfun(@(x) size(x,1), {datarois(:).roiMatevokedstim});
            vectcellssitestmp = cellfun(@(x,y) ones(1,size(x,1))*y, {datarois(:).roiMatevokedstim},mat2cell([1:length({datarois(:).roiMatevokedstim})],1,ones(1,length({datarois(:).roiMatevokedstim}))),'UniformOutput',0);
            vectcellssites = cat(2,vectcellssitestmp{:});
            medianconcatmatevokedFMB = squeeze(median(concatmatevokedFMB,3));
            medianconcatmatevokedEcho = squeeze(median(concatmatevokedEcho,3));
            tmpconcamat = cat(2,concatmatevokedFMB,concatmatevokedEcho);
            cumnbcellspersite = [0,cumsum(cellfun(@(x) size(x,1),{datarois(:).roiMatevokedstim}))];
            siteid = [datarois(:).site];
            for cc = 1:size(concatmatevokedFMB,1)
                clear tmppsth
                for ll = 1:10
                    tmppsth{ll} = tmpconcamat(cc,datasample([1:max(Echoind)],max(Echoind),'Replace',false),datasample([1:size(tmpconcamat,3)],10,'Replace',false),:);
                end
                shuffleconcatmatevokedFMB(cc,:,:,:) = mean(cell2mat(cellfun(@(x) x(:,Socialind,:,:),tmppsth,'UniformOutput',0)'));
                shuffleconcatmatevokedEcho(cc,:,:,:) = mean(cell2mat(cellfun(@(x) x(:,Echoind,:,:),tmppsth,'UniformOutput',0)'));
            end
            shufflemedianconcatmatevokedFMB = squeeze(median(shuffleconcatmatevokedFMB,3));
            shufflemedianconcatmatevokedEcho = squeeze(median(shuffleconcatmatevokedEcho,3));
            
            
            [sortconcat,sortconcatind] = sort(squeeze(max(permute(squeeze(mean(medianconcatmatevokedFMB(:,:,responsewindow),2)),[2,1]))));
            concatmatevokedFMBBcorr = concatmatevokedFMB-mean(concatmatevokedFMB(:,:,:,baselinewindow),4);
            concatmatevokedEchoBcorr = concatmatevokedEcho-mean(concatmatevokedEcho(:,:,:,baselinewindow),4);
            medianconcatmatevokedFMBBcorr = squeeze(median(concatmatevokedFMBBcorr,3));
            medianconcatmatevokedEchoBcorr = squeeze(median(concatmatevokedEchoBcorr,3));
            shuffleconcatmatevokedFMBBcorr = shuffleconcatmatevokedFMB-mean(shuffleconcatmatevokedFMB(:,:,:,1:10),4);
            shuffleconcatmatevokedEchoBcorr = shuffleconcatmatevokedEcho-mean(shuffleconcatmatevokedEcho(:,:,:,1:10),4);
            shufflemedianconcatmatevokedFMBBcorr = squeeze(median(shuffleconcatmatevokedFMBBcorr,3));
            shufflemedianconcatmatevokedEchoBcorr = squeeze(median(shuffleconcatmatevokedEchoBcorr,3));
            
            medrespFMB = median(medianconcatmatevokedFMBBcorr(:,:,responsewindow),3);
            medrespEcho =  median(medianconcatmatevokedEchoBcorr(:,:,responsewindow),3);
            shufflemedrespFMB = median(shufflemedianconcatmatevokedFMBBcorr(:,:,responsewindow),3);
            shufflemedrespEcho =  median(shufflemedianconcatmatevokedEchoBcorr(:,:,responsewindow),3);
            MedindexSelect = (abs(medrespFMB)-abs(medrespEcho))./(abs(medrespFMB)+abs(medrespEcho));
            shuffleMedindexSelect = (abs(shufflemedrespFMB)-abs(shufflemedrespEcho))./(abs(shufflemedrespFMB)+abs(shufflemedrespEcho));
            [sortmedindval,sortmedindind] = sort(mean(MedindexSelect,2));
            
            figure(); hold on; ax(1)=subplot(2,3,1); hold on; imagesc(xvalues,1:length(sortmedindind),squeeze(mean(medianconcatmatevokedFMB(fliplr(sortmedindind'),:,:),2))); colormap(colorredwhite); caxis([1,1.1]);
            plot([0,0],[1,length(sortmedindind)],'k'); axis tight
            title('Social mean'); colorbar; xlabel('Time from onset [s]')
            
            ax(2)=subplot(2,3,2); hold on; imagesc(xvalues,1:length(sortmedindind),squeeze(mean(medianconcatmatevokedEcho(fliplr(sortmedindind'),:,:),2))); colormap(colorredwhite); caxis([1,1.1]);
            plot([0,0],[1,length(sortmedindind)],'k'); axis tight
            title('Echo mean'); colorbar; xlabel('Time from onset [s]')
            
            %             subplot(1,5,3);imagesc(squeeze(mean(medianconcatmatevokedFMB(fliplr(sortmedindind'),:,:),2))-squeeze(mean(medianconcatmatevokedEcho(fliplr(sortmedindind'),:,:),2))); colormap(redblue(100)); caxis([-0.2,0.2]);
            %             title('Mean evoked response FMBs- Echo'); colorbar
            % check your selectivity index
            %             subplot(1,5,4); imagesc(squeeze(mean(medianconcatmatevokedFMBBcorr(fliplr(sortmedindind'),:,:),2))-squeeze(mean(medianconcatmatevokedEchoBcorr(fliplr(sortmedindind'),:,:),2))); colormap(redblue(100)); caxis([-0.1,0.1]);
            %             title('Mean evoked response baseline corrected FMBs- Echo'); colorbar
            ax(3)=subplot(2,3,3); hold on; imagesc(xvalues,1:length(sortmedindind),squeeze(mean(medianconcatmatevokedFMBBcorr(fliplr(sortmedindind'),:,:),2)-mean(medianconcatmatevokedEchoBcorr(fliplr(sortmedindind'),:,:),2))); colormap(redblue(100)); caxis([-0.1,0.1]);
            plot([0,0],[1,length(sortmedindind)],'k'); axis tight
            title('FMBs - Echo'); colorbar; xlabel('Time from onset [s]')
            xlabel('Time from onset [s]')
            
            a(4)= subplot(2,3,[4]); hold on;
            plot(mean(shuffleMedindexSelect(sortmedindind),2),1:size(sortmedindind,1),'.','Color',[0.5,0.5,0.5]);
            plot(sortmedindval,1:size(sortmedindind,1),'.','Color',[0,0,0]); ylim([1, length(sortmedindind)]);%colormap(redblue(100)); caxis([-0.2,0.2]);
            xlabel('SI value')
            title(sprintf(['selectivity index','\n','(median response window)']))
            
            [histshufflebin,histshufflevalue] = hist(mean(shuffleMedindexSelect,2),100);
            boundarymedselect = prctile(mean(shuffleMedindexSelect,2)',[2.5,97.5]);
            [histbin,histsvalue] = hist(sortmedindval,100);
            ax(5)=subplot(2,3,5); hold on;
            plot(histshufflevalue,histshufflebin,'Color',[0.5,0.5,0.5],'LineWidth',1.5);
            plot(histsvalue,histbin,'k','LineWidth',1.5); xlabel('SI value')
            plot([boundarymedselect(1),boundarymedselect(1)],[0,max(histshufflebin)],'--','Color',[0.5,0.5,0.5])
            plot([boundarymedselect(2),boundarymedselect(2)],[0,max(histshufflebin)],'--','Color',[0.5,0.5,0.5]); ylim([0,max(histshufflebin)])
            xlim([-1,1]);legend('shuffle','data'); legend boxoff
            title('Distribution SI')
            colormap(ax(1),colorredwhite)
            colormap(ax(2),colorredwhite)
            % add % of selective cells per site
            ax(6)=subplot(2,3,6); hold on;
            echoselectind = sortmedindind(find(sortmedindval<=boundarymedselect(1))); 
            socialselectind = sortmedindind(find(sortmedindval>=boundarymedselect(2))); 
            xbar = 4:4:(length({datarois(:).roiMatevokedstim})+1)*4;
            prtcechotot = (numel(echoselectind)./numel(sortmedindval))*100;
            prtcsocialtot = (numel(socialselectind)./numel(sortmedindval))*100;            
            for st = 1:length({datarois(:).roiMatevokedstim})
              prctecho(st) = (sum(ismember(echoselectind,find(vectcellssites == st)))./numel(find(vectcellssites == st)))*100;
              prtcsocial(st) = (sum(ismember(socialselectind,find(vectcellssites == st)))./numel(find(vectcellssites == st)))*100;
              bar(xbar(st)-2,prctecho(st),'FaceColor',[0 .2 .8])
              bar(xbar(st)-1,prtcsocial(st),'FaceColor',[0.8 .2 0])
            end
            a = bar(xbar(st+1)-2,prtcechotot,'FaceColor',[0 .2 0.8])
            b = bar(xbar(st+1)-1,prtcsocialtot,'FaceColor',[0.8 .2 0])
            xticks([4:4:((length({datarois(:).roiMatevokedstim})+1)*4)]-1.5)
            xticklabels([cellfun(@num2str,num2cell(1:length({datarois(:).roiMatevokedstim})),'UniformOutput',0),{'Total'}])
            xlabel('site');ylabel('% selective cells'); ylim([0,50]);
            title('% selective per site')
            legend([a,b],{'Echo','Social'}); legend boxoff
 
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\SummarySelectivityAllSites.pdf'],'-dpdf','-fillpage')
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\SummarySelectivityAllSites.svg'],'-dsvg')
            
            % distance and selectivity index % do it per site
            nbsite = length([datarois(:).site]);
            siteveccellstmp = cellfun(@(x,y) ones(1,size(x,1))*y,{datarois(:).roiMatevokedstim},num2cell(1:nbsite),'UniformOutput',0);
            siteveccells = [siteveccellstmp{:}];
            for st = 1:nbsite
                for cc = find(siteveccells==st)
                    for cci = max(find(siteveccells==st)):-1:cc
                        indcc = (cc-find(siteveccells==st,1))+1;
                        indcci = (cci-find(siteveccells==st,1))+1;
                        %                     disp(num2str(cc))
                        if mean(MedindexSelect(cc,:))<boundarymedselect(1) || mean(MedindexSelect(cc,:))>boundarymedselect(2) || mean(MedindexSelect(cci,:))>boundarymedselect(2) || mean(MedindexSelect(cci,:))<boundarymedselect(1)
                            pariselect{st}(indcc,indcci) =  mean(MedindexSelect(cc,:))+ mean(MedindexSelect(cci,:));
                        else
                            pariselect{st}(indcc,indcci) = nan;
                        end
                    end
                end
            end
            dist = {datarois(:).dist};
            for st = 1:nbsite
                vecdist{st} = reshape(dist{st},1,size(dist{st},1)*size(dist{st},2));
                vecsi{st} = reshape(abs(pariselect{st})/2,1,size(dist{st},1)*size(dist{st},2));
                indzeronan = unique([find(vecdist{st}==0),find(isnan(vecsi{st}))]);
                vecdist{st}(indzeronan) = [];
                vecsi{st}(indzeronan) = [];
                [discretebins,discreteedge] = discretize(vecdist{st},10);
                for i = unique(discretebins)
                    meansi(i) = mean(vecsi{st}(find(discretebins==i)));
                    semsi(i) = std(vecsi{st}(find(discretebins==i)));%/(sqrt(numel(find(discretebins==i))));
                end
                figure(st); hold on; subplot(1,2,1); hold on
%                 plot(discretebins,vecsi{st} ,'.','Color',[0.5,0.5,0.5]);
                shadedErrorBar(unique(discretebins),meansi,semsi,{'Color',[0,0,0],'LineWidth',2},1);
                xticks(unique(discretebins));
                xticklabels(num2cell(discreteedge));
                ylabel('Similarity Index'); xlabel('Distance pair [/mum]');
%                 linfit = fitlm(discretebins,vecsi{st},'linear');
%                 refline(linfit.Coefficients.Estimate(2),linfit.Coefficients.Estimate(1));
%                 dim = [.15 .65 .1 .1]; ylim([0,0.4]); xlim([1,10]); axis square
%                 annotation('textbox',dim,'String',['p = ' num2str(round(linfit.Coefficients.pValue(2)*10000)/10000) ' r =' num2str(round(linfit.Rsquared.Ordinary*10000)/10000)],'EdgeColor','none')
                %                 legend('data','fit'); legend box off
                % do equally sized bins
                [a,b,c] = histcounts(vecdist{st},prctile(vecdist{st},[0:10:100]));
                for i = unique(c)
                    meansieq(i) = mean(vecsi{st}(find(c==i)));
                    semsieq(i) = std(vecsi{st}(find(c==i)));%/(sqrt(numel(find(c==i))));
                end
                subplot(1,2,2); hold on;
%                 plot(c,vecsi{st} ,'.','Color',[0.5,0.5,0.5]);
                shadedErrorBar(unique(c),meansieq,semsieq,{'Color',[0,0,0],'LineWidth',2},1);
                
                %                 plot(c,vecsi{st} ,'.','Color',[0.5,0.5,0.5]);
                xticks(unique(c));
                xticklabels(num2cell(round(b(2:end))));
                linfit = fitlm(c,vecsi{st},'linear');
                refline(linfit.Coefficients.Estimate(2),linfit.Coefficients.Estimate(1))
                ylabel('Similarity Index'); xlabel('Distance pair [microns]');
                dim = [.6 .65 .1 .1];  ylim([0,0.4]); xlim([1,10]); axis square
                annotation('textbox',dim,'String',['p = ' num2str(round(linfit.Coefficients.pValue(2)*100)/100) ' r =' num2str(round(linfit.Rsquared.Ordinary*10000)/10000)],'EdgeColor','none')
%                 print(['N:\Jenni\' animal '\figure\social\Population\site' num2str(datarois(st).site) '\SimilaritySleectivityDistance.pdf'],'-dpdf','-fillpage')
%                 saveas(gcf,['N:\Jenni\' animal '\figure\social\Population\site' num2str(datarois(st).site) '\SimilaritySleectivityDistance.svg'])
%                 
            end
            
            % population anlaysis
            vecdistall = [vecdist{:}];
            vecsiall = [vecsi{:}];
            [discretebins,discreteedge] = discretize(vecdistall,10);
            figure(st); hold on; subplot(1,2,1);
            plot(discretebins,vecsiall ,'.','Color',[0.5,0.5,0.5]);
            xticks(unique(discretebins));
            xticklabels(num2cell(discreteedge));
            ylabel('Similarity Index'); xlabel('Distance pair [microns]');
            %             linfit = fitlm(discretebins,vecsi,'linear'); ylim([0,1]);
            %             refline(linfit.Coefficients.Estimate(2),linfit.Coefficients.Estimate(1));
            %             dim = [.15 .8 .1 .1];
            %             annotation('textbox',dim,'String',['p = ' num2str(round(linfit.Coefficients.pValue(2)*10000)/10000) ' r =' num2str(round(linfit.Rsquared.Ordinary*10000)/10000)])
            % do equally sized bins
            
            [a,b,c] = histcounts(vecdistall,prctile(vecdistall,[0:10:100]));
            for i = unique(c)
                meansieqall(i) = mean(vecsiall(find(c==i)));
                semsieqall(i) = std(vecsiall(find(c==i)));%(sqrt(numel(find(c==i))));
            end
            figure(); hold on; subplot(1,2,2); hold on; 
            shadedErrorBar(unique(c),meansieqall,semsieqall,{'Color',[0,0,0],'LineWidth',2});
            xticks(unique(c)); ylim([0,0.4]); xlim([1,10]);
            xticklabels(num2cell(round(b(2:end))));
            linfit = fitlm(c,vecsiall,'linear');
            refline(linfit.Coefficients.Estimate(2),linfit.Coefficients.Estimate(1))
            ylabel('Similarity Index'); xlabel('Distance pair [microns]');
            dim = [.6 .8 .1 .1]; ylim([0,0.4]);
            annotation('textbox',dim,'String',['p = ' num2str(round(linfit.Coefficients.pValue(2)*100)/100) ' r =' num2str(round(linfit.Rsquared.Ordinary*10000)/10000)],'EdgeColor','none')
            
            
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\DistSelectivityAllSites.pdf'],'-dpdf','-bestfit')
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\DistSelectivityAllSites.png'],'-dpng')

            print([DataDrive ':\Jenni\' animal '\figure\social\Population\DistSelectivityAllSites.svg'],'-dsvg')

            % decoding for all sites % tmp concamat only contains soxial
            % and echo stim
            loopvalue = 100;%0;
            correctedstim = datarois(fn).correctedstim;
            for ii = 1:size(tmpconcamat,2)
                for iif = size(tmpconcamat,2):-1:1
                    for lp = 1:loopvalue
                        clear RadomSelectedTarget roiMatevokedstim concaactivity RadomSelectedFoil testtriallength vec1 vecsi b TestFoil TestTarget VectId classif
                        roiMatevokedstim = tmpconcamat;
                        if ii == iif
                            traintriallength = (size(roiMatevokedstim,3)/2)-1;
                            RadomSelectedTarget = datasample(1:size(roiMatevokedstim,3),traintriallength,'Replace',false);
                            RadomSelectedFoil = datasample(find(~ismember(1:size(roiMatevokedstim,3),RadomSelectedTarget)),traintriallength,'Replace',false);
                            TestTarget = (find(~ismember(1:size(roiMatevokedstim,3),[RadomSelectedTarget,RadomSelectedFoil])));
                            TestFoil = find(~ismember(1:size(roiMatevokedstim,3),[RadomSelectedTarget,RadomSelectedFoil,TestTarget]));
                        else
                            traintriallength = size(roiMatevokedstim,3)-1;
                            RadomSelectedTarget = datasample(1:size(roiMatevokedstim,3),traintriallength,'Replace',false);
                            RadomSelectedFoil = datasample(1:size(roiMatevokedstim,3),traintriallength,'Replace',false);
                            TestTarget = 1:size(roiMatevokedstim,3);
                            TestTarget(ismember(TestTarget,RadomSelectedTarget)) = [];
                            
                            TestFoil = 1:size(roiMatevokedstim,3);
                            TestFoil(ismember(TestFoil,RadomSelectedFoil)) = [];
                        end
                        
                        
                        TestTrials = [TestTarget,TestFoil];
                        VectId = [ones(1,length(TestTarget)),ones(1,length(TestFoil))*2];
                        
                        vec1 = mean(squeeze(mean(roiMatevokedstim(:,ii,RadomSelectedTarget,responsewindow),4)),2);
                        vecsi = mean(squeeze(mean(roiMatevokedstim(:,iif,RadomSelectedFoil,responsewindow),4)),2);
                        w{lp} = vec1-vecsi;
                        b = -(dot(w{lp},vec1) + dot(w{lp},vecsi))/2;
                        
                        %                                 shufflevec1{lp} = squeeze(mean(mean(Shuffletonetc{ii}(RadomSelectedTarget,targetwindow+10,:),2)));
                        %                                 shufflevec2{lp} = squeeze(mean(mean(Shuffletonetc{ii}(RadomSelectedFoil,baselinewindow+10,:),2)));
                        %                                 shufflew{lp} = vec1{lp}-vec2{lp};
                        %                                 shuffleb = -(dot(shufflew{lp},shufflevec1{lp}) + dot(shufflew{lp},shufflevec2{lp}));
                        concaactivity = [squeeze(mean(mean(roiMatevokedstim(:,ii,TestTarget,responsewindow),4),2)),squeeze(mean(mean(roiMatevokedstim(:,iif,TestFoil,responsewindow),4),2))];
                        %                                 shuffleconcaactivity = cat(1,mean(Shuffletonetc{ii}(TestTarget,targetwindow+10,:),2),mean(Shuffletonetc{ii}(TestFoil,baselinewindow+10,:),2));
                        
                        for i = 1:length(TestTrials)
                            classif(i) = dot(squeeze(concaactivity(:,i)),w{lp})+b;
                            %                                         shuffleclassif(i) = dot(squeeze(shuffleconcaactivity(i,1,:)),shufflew{lp})+shuffleb;
                        end
                        
                        ClassifiedTarget = find(classif>0);
                        ClassifiedFoil = find(classif<0);
                        %                                 ShuffleClassifiedTarget = find(shuffleclassif>0);
                        %                                 ShuffleClassifiedFoil = find(shuffleclassif<0);
                        %
                        DecodingAccuracytonetc(ii,iif,lp) = (numel(intersect(ClassifiedTarget,find(VectId==1)))+ numel(intersect(ClassifiedFoil,find(VectId==2))))/(length(TestTarget)+length(TestFoil));
                        %                                 ShuffleDecodingAccuracytonetc{ii}(lp,1) = (numel(intersect(ShuffleClassifiedTarget,find(VectId==1)))+ numel(intersect(ShuffleClassifiedFoil,find(VectId==2))))/(length(TestTarget)+length(TestFoil));
                        DecodingWeight(ii,iif,:,lp) = w{lp};
                    end
                    InterestingNeurons{ii,iif} = [find(mean(DecodingWeight(ii,iif,:,:),4)>0.1);find(mean(DecodingWeight(ii,iif,:,:),4)<-0.1)];
                    InterestingNeuronsvalue{ii,iif} = squeeze(squeeze(mean(DecodingWeight(ii,iif,InterestingNeurons{ii,iif}),4)));
                end
            end
            %             colorredwhite = flipud([ones(100,1),(log10([1.01:0.1:11])/log10(11))',(log10([1.01:0.1:11])/log10(11))']);%([0:0.01:1-0.01].^1.2)',([0:0.01:1-0.01].^1.05)'];
            t = num2cell((squeeze(mean(DecodingAccuracytonetc*100,3)))); % extact values into cells
            x = repmat(1:10,10,1); % generate x-coordinates
            y = x'; % generate y-coordinates
            t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
            % Draw Image and Label Pixels
            h = figure(); hold on; subplot(3,3,[1:2,4:5,7:8]); hold on; 
            imagesc((squeeze(mean(DecodingAccuracytonetc*100,3)))); caxis([50,100]); colormap(colorredwhite); axis tight;
            %             text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
            ylabel('stim Id'); xlabel('stim Id');
            set(gca,'YDir','reverse')
            xticks([1:size(tmpconcamat,2)]) 
            yticks((([1:length(correctedstim(1:size(tmpconcamat,2)))])))
            if size(tmpconcamat,2)==12
                xticklabels([cellfun(@(x) x(5:8),{stimwavname(correctedstim(1:6)).name},'UniformOutput',0),cellfun(@(x) x(9:13),{stimwavname(correctedstim(7:12)).name},'UniformOutput',0)])%xticklabels({stimwavname(correctedstim(1:size(tmpconcamat,2))).name}) % modify to get correct names for when the sequence was out of order
                yticklabels([cellfun(@(x) x(5:8),{stimwavname(correctedstim(1:6)).name},'UniformOutput',0),cellfun(@(x) x(9:13),{stimwavname(correctedstim(7:12)).name},'UniformOutput',0)])%xticklabels({stimwavname(correctedstim(1:size(tmpconcamat,2))).name}) % modify to get correct names for when the sequence was out of order
            elseif size(tmpconcamat,2)==10
                xticklabels([cellfun(@(x) x(5:8),{stimwavname(correctedstim(1:5)).name},'UniformOutput',0),cellfun(@(x) x(9:13),{stimwavname(correctedstim(6:10)).name},'UniformOutput',0)])%xticklabels({stimwavname(correctedstim(1:size(tmpconcamat,2))).name}) % modify to get correct names for when the sequence was out of order
                yticklabels([cellfun(@(x) x(5:8),{stimwavname(correctedstim(1:5)).name},'UniformOutput',0),cellfun(@(x) x(9:13),{stimwavname(correctedstim(6:10)).name},'UniformOutput',0)])%xticklabels({stimwavname(correctedstim(1:size(tmpconcamat,2))).name}) % modify to get correct names for when the sequence was out of order
            end
            xtickangle(45)
            ytickangle(45)
            colorbar;
            axis square
            title('Pairwise decoding accuracy')
            % decode simply echo vs social 
            
            % plot neural trajectory (stim is like 3 frames so let's do first 10 frames)
            % compare weight decoder and selectivity
            samestim = [ones(size(tmpconcamat,2),2).*(1:size(tmpconcamat,2))'];
            for dci = 1:size(tmpconcamat,2)
            samestimdecodeweight(:,dci) = squeeze(mean(DecodingWeight(samestim(dci,1),samestim(dci,2),:,:),4));
            end
            for dci = 1:size(Socialind,2)
            diffstimdecodeweight(:,:,dci) = squeeze(mean(DecodingWeight(dci,Echoind,:,:),4));
            end
            subplot(3,3,3); hold on;
            for dci = 1:size(Socialind,2)
            plot(1:size(samestimdecodeweight,1),abs(squeeze(diffstimdecodeweight(:,:,dci))'),'.','Color',[0,0,0])
            end
            plot(1:size(samestimdecodeweight,1),abs(samestimdecodeweight),'.','Color',[0.5,0.5,0.5])
            xlim([1,size(samestimdecodeweight,1)]);
            xlabel('cell nb')
            ylabel('Decoding weights')
            [weigthssamehist,weigthssamehistval] = hist(reshape((samestimdecodeweight),1,size(samestimdecodeweight,1)*size(samestimdecodeweight,2)),40);
            [weigthsdiffhist,weigthsdiffhistval] = hist(reshape((diffstimdecodeweight),1,size(diffstimdecodeweight,1)*size(diffstimdecodeweight,2)*size(diffstimdecodeweight,3)),40);
            subplot(3,3,6); hold on;
            plot(weigthssamehistval,weigthssamehist./max(weigthssamehist),'Color',[0.5,0.5,0.5],'LineWidth',2);
            plot(weigthsdiffhistval,weigthsdiffhist./max(weigthsdiffhist),'Color',[0,0,0],'LineWidth',2);
            xlabel('Weights pairwise decoding over loop');
            legend('Same stim','Across cat'); legend box off
%             subplot(2,3,3); hold on;
%             plot(median(abs(samestimdecodeweight)'),'.','Color',[0.5,0.5,0.5])
%             for dci = 1:size(Socialind,2)
%             plot(median(abs(squeeze(diffstimdecodeweight(:,:,1)))),'.','Color',[0,0,0])
%             end
            subplot(3,3,9); hold on;
            plot(squeeze(mean(mean((diffstimdecodeweight),3))),(mean((shuffleMedindexSelect),2)),'.','Color',[0.5,0.5,0.5])
            plot(squeeze(mean(mean((diffstimdecodeweight),3))),(mean((MedindexSelect),2)),'.','Color',[0,0,0])
            axis square
            xlabel(['Mean pairwise decoding weigths'])
            ylabel('Selectivity index')
            legend('Shuffle','Data'); legend boxoff
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\AllSites_PairwiseDecodingAccuracy.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\AllSites_PairwiseDecodingAccuracy.svg'],'-dsvg')
            
            
            % decoding as a function of time % multiclass decoder %
            % category decoder to start for category 1 vs 2
            nbcells = size(tmpconcamat,1);
            nbtrials = (size(tmpconcamat,2)/2)*size(tmpconcamat,3);
            socialtmptime = reshape(tmpconcamat(:,1:size(tmpconcamat,2)/2,:,:),nbcells,nbtrials,length(xvalues));
            echotmptime = reshape(tmpconcamat(:,(size(tmpconcamat,2)/2)+1:end,:,:),nbcells,nbtrials,length(xvalues));
            for tt = 1:length(xvalues)
                for lp = 1:loopvalue
                    clear RadomSelectedTarget tmptrial roiMatevokedstim concaactivity SOconcaactivity RadomSelectedFoil testtriallength vecso vecec b TestFoil TestTarget VectId classif
                    RadomSelectedTarget = datasample(1:nbtrials,nbtrials/2,'Replace',false);
                    RadomSelectedFoil = datasample(1:nbtrials,nbtrials/2,'Replace',false);
                    TestTarget = 1:nbtrials;
                    TestTarget(ismember(TestTarget,RadomSelectedTarget)) = [];
                    TestFoil = 1:nbtrials;
                    TestFoil(ismember(TestFoil,RadomSelectedFoil)) = [];
                    TestTrials = [TestTarget,TestFoil];
                    VectId = [ones(1,length(TestTarget)),ones(1,length(TestFoil))*2];
                    
                    vecso = squeeze(mean(socialtmptime(:,RadomSelectedTarget,tt),2));
                    vecec = squeeze(mean(echotmptime(:,RadomSelectedFoil,tt),2));
                    wcattime{lp}(tt,:) = vecso-vecec;
                    b = -(dot(wcattime{lp}(tt,:),vecso) + dot(wcattime{lp}(tt,:),vecec))/2;
                    concaactivity = [squeeze(socialtmptime(:,TestTarget,tt)),squeeze(echotmptime(:,TestFoil,tt))];
                    
                    for i = 1:length(TestTrials)
                        classif(i) = dot(squeeze(concaactivity(:,i)), wcattime{lp}(tt,:))+b;
                    end
                    
                    ClassifiedTarget = find(classif>0);
                    ClassifiedFoil = find(classif<0);
                    DecodingAccuracytimecat(tt,lp) = (numel(intersect(ClassifiedTarget,find(VectId==1)))+ numel(intersect(ClassifiedFoil,find(VectId==2))))/(length(TestTarget)+length(TestFoil));
                    DecodingWeighttimecat(:,tt,lp) = wcattime{lp}(tt,:);

%                     % do a within cat to check for so 1st
                    RadomSelectedTargetSO =  datasample(1:nbtrials,floor(nbtrials/4),'Replace',false);
                    tmptrial = 1:nbtrials;
                    tmptrial(ismember(tmptrial,RadomSelectedTargetSO)) = [];                     
                    RadomSelectedFoilSO =  datasample(tmptrial,floor(nbtrials/4),'Replace',false);
                    tmptrial(ismember(tmptrial,RadomSelectedFoilSO)) = [];   
                    SOTestTarget = datasample(tmptrial,floor(nbtrials/4),'Replace',false);
                    tmptrial(ismember(tmptrial,SOTestTarget)) = [];   
                    SOTestFoil = tmptrial;
                    vec1so = squeeze(mean(socialtmptime(:,RadomSelectedTargetSO,tt),2));
                    vec2so = squeeze(mean(socialtmptime(:,RadomSelectedFoilSO,tt),2));
                    wcattimeso{lp}(tt,:) = vec1so-vec2so;
                    b = -(dot(wcattimeso{lp}(tt,:),vec1so) + dot(wcattimeso{lp}(tt,:),vec2so))/2;
                    SOconcaactivity = [squeeze(socialtmptime(:,SOTestTarget,tt)),squeeze(socialtmptime(:,SOTestFoil,tt))];
                    SOVectId = [ones(1,floor(nbtrials/4)),ones(1,floor(nbtrials/4))*2];
                    for i = 1:length(nbtrials/2)
                        classifwithin(i) = dot(squeeze(SOconcaactivity(:,i)), wcattimeso{lp}(tt,:))+b;
                    end
                    ClassifiedTarget = find(classif>0);
                    ClassifiedFoil = find(classif<0);
                    SODecodingAccuracytimecat(tt,lp) = (numel(intersect(ClassifiedTarget,find(SOVectId==1)))+ numel(intersect(ClassifiedFoil,find(SOVectId==2))))/(length(SOTestTarget)+length(SOTestFoil));
                    SODecodingWeighttimecat(:,tt,lp) = wcattimeso{lp}(tt,:);
                end
            end
            plot decoding accuracy
            h = figure(); hold on; subplot(3,1,[1,2]); hold on;
            a = shadedErrorBar(xvalues,mean(DecodingAccuracytimecat*100,2)',std((DecodingAccuracytimecat*100)'),{'k'},1);
            xlim([xvalues(1),xvalues(end)]); ylim([20,100]); 
            xlabel('Time from sound onset [s]'); ylabel('Decoding accuracy [%]')
            b = shadedErrorBar(xvalues,mean(SODecodingAccuracytimecat*100,2)',std((SODecodingAccuracytimecat*100)'),{'Color',[0.8,0.3,0.2]},1);
            % plot stimboundaries
            c = plot([xvalues(11),xvalues(11)],[0,100],'k--'); plot([xvalues(14),xvalues(14)],[0,100],'k--');
            legend([a.mainLine,b.mainLine],{'Social vs. Echo','Social vs. Social'}); legend boxoff
            % plot PSTH to compare the rise of the fluorescence indicator
            subplot(3,1,3); hold on; 
            plot(xvalues,median(squeeze(mean(socialtmptime))),'r');
            plot(xvalues,median(squeeze(mean(echotmptime))),'b');
            xlim([xvalues(1),xvalues(end)]); ylim([0.95,1.3]); 
            xlabel('Time from sound onset [s]'); ylabel('df/f')
            plot([xvalues(11),xvalues(11)],[0,100],'k--'); plot([xvalues(14),xvalues(14)],[0,100],'k--');
            legend('social','echo','start stim','end stim'); legend box off
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\AllSites_timedecoding.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\AllSites_timedecoding.svg'],'-dsvg')
%             
            % load abffile to compare the time elapse between sound onset
            % and frame % Identify number of calls per frames
            %  maybe add the tonotopy curve just to check the onset
            %  response the trace needs to be correlated with nb of pulse
            
            % plot tonotopy against SI
            notselectiveind = 1:nbcells;
            notselectiveind([echoselectind;socialselectind]) = [];
            load([animaldrive ':\Jenni\' animal '\z-stack\Tonotopystruct.mat']) % first 5 sites for RFB2 for the moment
            alltuning = cat(2,Tonotopystruct.cellstuning);
            mediandftuning = (cell2mat(cellfun(@(x) max(x'),[Tonotopystruct.mediantuningcurve],'UniformOutput',0)'));
            correchoind = echoselectind;%intersect(echoselectind,find(ismember(vectcellssites,[4:8])))-1192; % 4 to 8 for RoofBuddy1 -1192
            corrsocialind = socialselectind;%intersect(socialselectind,find(ismember(vectcellssites,[4:8])))-1192;
            for nbtones = 1:max(alltuning(3,correchoind))
            histechotuning(nbtones) = sum(alltuning(3,correchoind)==nbtones);
            histsocialtuning(nbtones) = sum(alltuning(3,corrsocialind)==nbtones);
            histtuning(nbtones) = sum(alltuning(3,:)==nbtones);
            end
            histechotuning(nbtones+1)=sum(alltuning(3,correchoind)==0);
            histsocialtuning(nbtones+1)= sum(alltuning(3,corrsocialind)==0);
            histtuning(nbtones+1)= sum(alltuning(3,:)==0);

            h=figure(); hold on; subplot(1,2,1); hold on;
            plot([histechotuning(end),histechotuning(1:end-1)],'Color',[0.2,0.2,0.8],'LineWidth',2);
            plot([histsocialtuning(end),histsocialtuning(1:end-1)],'Color',[0.8,0.2,0.2],'LineWidth',2);
%             plot([histtuning(end),histtuning(1:end-1)],'Color',[0,0,0],'LineWidth',2);
            ylabel('Nb cells')
            xlabel('Best freq.')
            xticks(1:nbtones+1)
            xticklabels({'No bf','4','5.7','8','11.3','16','22.7','32','45.2','64','90.5'})
            xlim([1,11]);
            legend('Echo','Social'); legend boxoff
            title('Hist tuning for selective cells')
            subplot(1,2,2); hold on;
            plot([histechotuning(end),histechotuning(1:end-1)]./max(histechotuning),'Color',[0.2,0.2,0.8],'LineWidth',2);
            plot([histsocialtuning(end),histsocialtuning(1:end-1)]./max(histsocialtuning),'Color',[0.8,0.2,0.2],'LineWidth',2);
            plot([histtuning(end),histtuning(1:end-1)]./max(histtuning),'Color',[0,0,0],'LineWidth',2);
            ylabel('Nb cells [norm.]')
            xlabel('Best freq.');
            xlim([1,11]);
            xticks(1:nbtones+1)
            title('Hist tuning for selective cells [norm.]')
            xticklabels({'No bf','4','5.7','8','11.3','16','22.7','32','45.2','64','90.5'})
            legend('Echo','Social','All'); legend boxoff
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_Tuning.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_Tuning.svg'],'-dsvg')
        
            % due the tuning selectivity but with df/f instead of bf
            mediantuningPSTH = [Tonotopystruct.mediantuningcurve];
            mediantuningPSTHnorm = cellfun(@(x) (x-min(x))./(max(x)-min(x)),mediantuningPSTH,'UniformOutput',0);
            maxmediantuning = cell2mat(cellfun(@(x) max(x(:,responsewindow)'),mediantuningPSTH,'UniformOutput',0)')
            maxmediantuningnorm = (maxmediantuning'-min(maxmediantuning'))./(max(maxmediantuning')-min(maxmediantuning'));
            maxmediantuningnorm = maxmediantuningnorm';
            h=figure(); hold on; subplot(1,2,1); hold on;
%             plot(mean(maxmediantuning(echoselectind,:)),'Color',[0.2,0.2,0.8],'LineWidth',2);
%             plot(mean(maxmediantuning(socialselectind,:)),'Color',[0.8,0.2,0.2],'LineWidth',2);
%             plot([histtuning(end),histtuning(1:end-1)],'Color',[0,0,0],'LineWidth',2);
            shadedErrorBar(1:nbtones,mean(maxmediantuning(echoselectind,:)),std(maxmediantuning(echoselectind,:))./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            shadedErrorBar(1:nbtones,mean(maxmediantuning(socialselectind,:)),std(maxmediantuning(socialselectind,:))./sqrt(numel(socialselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            shadedErrorBar(1:nbtones,mean(maxmediantuning(notselectiveind,:)),std(maxmediantuning(notselectiveind,:))./sqrt(numel(notselectiveind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            ylabel('Mean peak df/f')
            xlabel('Tones [kHz]')
            xticks(1:nbtones)
            xticklabels({'4','5.7','8','11.3','16','22.7','32','45.2','64','90.5'})
            xlim([1,10]);
%             legend('Echo','Social'); legend boxoff
            title('Mean peak df/f tuning for selective cells')
            subplot(1,2,2); hold on;
%             plot(mean(maxmediantuningnorm(echoselectind,:)),'Color',[0.2,0.2,0.8],'LineWidth',2);
%             plot(mean(maxmediantuningnorm(socialselectind,:)),'Color',[0.8,0.2,0.2],'LineWidth',2);
%             plot(mean(maxmediantuningnorm(notselectiveind,:)),'Color',[0,0,0],'LineWidth',2);
            a = shadedErrorBar(1:nbtones,mean(maxmediantuningnorm(echoselectind,:)),std(maxmediantuningnorm(echoselectind,:))./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            b = shadedErrorBar(1:nbtones,mean(maxmediantuningnorm(socialselectind,:)),std(maxmediantuningnorm(socialselectind,:))./sqrt(numel(socialselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            c = shadedErrorBar(1:nbtones,mean(maxmediantuningnorm(notselectiveind,:)),std(maxmediantuningnorm(notselectiveind,:))./sqrt(numel(notselectiveind))*2,{'Color',[0,0,0],'LineWidth',2},1);

            ylabel('Mean peak df/f')
            xlabel('Tones [kHz]')
            xlim([1,10]);
            xticks(1:nbtones)
            title('Mean peak df/f tuning for selective cells [norm.]')
            xticklabels({'4','5.7','8','11.3','16','22.7','32','45.2','64','90.5'})
%             legend('Echo','Social','Not selective'); legend boxoff
            legend([a.mainLine,b.mainLine,c.mainLine],{'Echo','Social','Not selective'}); legend boxoff
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_Tuningdf.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_Tuningdf.svg'],'-dsvg')
        
            load([animaldrive ':\Jenni\' animal '\z-stack\Sweepstruct.mat']) % first 5 sites for RFB2 for the moment
            % response to white noise 
            tmpsweeppeak = {Sweepstruct.ampcells};
            sitenb = unique(vectcellssites);
            whpeak = cell2mat(cellfun(@(x) squeeze(x(1,:,:)),tmpsweeppeak(~cellfun(@isempty,tmpsweeppeak)),'UniformOutput',0));
            whpeaknorm = (whpeak-min(whpeak))./(max(whpeak)-min(whpeak));
            uppeak = cell2mat(cellfun(@(x) squeeze(x(2,:,:)),tmpsweeppeak(~cellfun(@isempty,tmpsweeppeak)),'UniformOutput',0));
            uppeaknorm = (uppeak-min(uppeak))./(max(uppeak)-min(uppeak));
            dspeak = cell2mat(cellfun(@(x) squeeze(x(3,:,:)),tmpsweeppeak(~cellfun(@isempty,tmpsweeppeak)),'UniformOutput',0));
            dspeaknorm = (dspeak-min(dspeak))./(max(dspeak)-min(dspeak));
            if strcmp(animal,'RoofBuddy2') % restrict cells to last 2 sites (6 and 7 corresponfing to 1 and 22)
            echoselectind = echoselectind(find(echoselectind>max(cumsum(nbcellssite(1:end-2)))))-max(cumsum(nbcellssite(1:end-2)));
            socialselectind = socialselectind(find(socialselectind>max(cumsum(nbcellssite(1:end-2)))))-max(cumsum(nbcellssite(1:end-2)));
            notselectiveind = notselectiveind(find(notselectiveind>max(cumsum(nbcellssite(1:end-2)))))-max(cumsum(nbcellssite(1:end-2)));
            end
            h=figure(); hold on; subplot(1,2,1); hold on;
            shadedErrorBar(1:8,mean(whpeak(:,echoselectind),2),std(whpeak(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            shadedErrorBar(1:8,mean(whpeak(:,socialselectind),2),std(whpeak(:,socialselectind)')./sqrt(numel(socialselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            shadedErrorBar(1:8,mean(whpeak(:,notselectiveind),2),std(whpeak(:,notselectiveind)')./sqrt(numel(notselectiveind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            ylabel('Peak df/f')
            xlabel('Duration [ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueduration})
            xlim([1,8]);
%             legend('Echo','Social'); legend boxoff
            title('WH resp. for selective cells')
            subplot(1,2,2); hold on;
            a = shadedErrorBar(1:8,mean(whpeaknorm(:,echoselectind),2),std(whpeaknorm(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            b = shadedErrorBar(1:8,mean(whpeaknorm(:,socialselectind),2),std(whpeaknorm(:,socialselectind)')./sqrt(numel(socialselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            c = shadedErrorBar(1:8,mean(whpeaknorm(:,notselectiveind),2),std(whpeaknorm(:,notselectiveind)')./sqrt(numel(notselectiveind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            ylabel('Peak df/f')
            xlabel('Duration [ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueduration})
            xlim([1,8]);
            legend([a.mainLine,b.mainLine,c.mainLine],{'Echo','Social','Not selective'}); legend boxoff
            title('WH resp. for selective cells [norm.]')
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_WH.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_WH.svg'],'-dsvg')
 
            h=figure(); hold on; subplot(1,2,1); hold on;
            shadedErrorBar(1:8,mean(uppeak(:,echoselectind),2),std(uppeak(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            shadedErrorBar(1:8,mean(uppeak(:,socialselectind),2),std(uppeak(:,socialselectind)')./sqrt(numel(socialselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            shadedErrorBar(1:8,mean(uppeak(:,notselectiveind),2),std(uppeak(:,notselectiveind)')./sqrt(numel(notselectiveind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            ylabel('Peak df/f')
            xlabel('Rate [kHz/ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueratekhz})
            xlim([1,8]);
%             legend('Echo','Social'); legend boxoff
            title('US resp. for selective cells')
            subplot(1,2,2); hold on;
            a = shadedErrorBar(1:8,mean(uppeaknorm(:,echoselectind),2),std(uppeaknorm(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            b = shadedErrorBar(1:8,mean(uppeaknorm(:,socialselectind),2),std(uppeaknorm(:,socialselectind)')./sqrt(numel(socialselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            c = shadedErrorBar(1:8,mean(uppeaknorm(:,notselectiveind),2),std(uppeaknorm(:,notselectiveind)')./sqrt(numel(notselectiveind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            ylabel('Peak df/f')
            xlabel('Rate [kHz/ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueratekhz})
            xlim([1,8]);
            legend([a.mainLine,b.mainLine,c.mainLine],{'Echo','Social','Not selective'}); legend boxoff
            title('US resp. for selective cells [norm.]')
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_US.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_US.svg'],'-dsvg')
            
            h=figure(); hold on; subplot(1,2,1); hold on;
            shadedErrorBar(1:8,mean(dspeak(:,echoselectind),2),std(dspeak(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            shadedErrorBar(1:8,mean(dspeak(:,socialselectind),2),std(dspeak(:,socialselectind)')./sqrt(numel(socialselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            shadedErrorBar(1:8,mean(dspeak(:,notselectiveind),2),std(dspeak(:,notselectiveind)')./sqrt(numel(notselectiveind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            ylabel('Peak df/f')
            xlabel('Rate [kHz/ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueratekhz})
            xlim([1,8]);
%             legend('Echo','Social'); legend boxoff
            title('DS resp. for selective cells')
            subplot(1,2,2); hold on;
            a = shadedErrorBar(1:8,mean(dspeaknorm(:,echoselectind),2),std(dspeaknorm(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            b = shadedErrorBar(1:8,mean(dspeaknorm(:,socialselectind),2),std(dspeaknorm(:,socialselectind)')./sqrt(numel(socialselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            c = shadedErrorBar(1:8,mean(dspeaknorm(:,notselectiveind),2),std(dspeaknorm(:,notselectiveind)')./sqrt(numel(notselectiveind))*2,{'Color',[0,0,0],'LineWidth',2},1);

            ylabel( 'Peak df/f')
            xlabel('Rate [kHz/ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueratekhz})
            xlim([1,8]);
            legend([a.mainLine,b.mainLine,c.mainLine],{'Echo','Social','Not selective'}); legend boxoff
            title('DS resp. for selective cells [norm.]')
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_DS.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_DS.svg'],'-dsvg')
            
            % Peak latency for all ds us and all
            
            latencypeaktmp = {Sweepstruct.peaklatency};
            whlatency = cell2mat(cellfun(@(x) squeeze(x(1,:,:)),latencypeaktmp(~cellfun(@isempty,latencypeaktmp)),'UniformOutput',0))/31.25;
            uslatency = cell2mat(cellfun(@(x) squeeze(x(2,:,:)),latencypeaktmp(~cellfun(@isempty,latencypeaktmp)),'UniformOutput',0))/31.25;
            dslatency = cell2mat(cellfun(@(x) squeeze(x(3,:,:)),latencypeaktmp(~cellfun(@isempty,latencypeaktmp)),'UniformOutput',0))/31.25;
            h = figure(); hold on; subplot(1,3,1); hold on; 
            a = shadedErrorBar(1:8,mean(whlatency(:,echoselectind),2),std(whlatency(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            b = shadedErrorBar(1:8,mean(whlatency(:,socialselectind),2),std(whlatency(:,socialselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            c = shadedErrorBar(1:8,mean(whlatency(:,notselectiveind),2),std(whlatency(:,notselectiveind)')./sqrt(numel(echoselectind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            title('Latency WH') 
            ylabel('Median Latency [s]')
            xlabel('Rate [kHz/ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueduration}); xlim([1,8])
            subplot(1,3,2); hold on; 
            a = shadedErrorBar(1:8,mean(uslatency(:,echoselectind),2),std(uslatency(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            b = shadedErrorBar(1:8,mean(uslatency(:,socialselectind),2),std(uslatency(:,socialselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            c = shadedErrorBar(1:8,mean(uslatency(:,notselectiveind),2),std(uslatency(:,notselectiveind)')./sqrt(numel(echoselectind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            title('Latency US') 
            ylabel('Median Latency [s]')
            xlabel('Rate [kHz/ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueratekhz}); xlim([1,8])
                  subplot(1,3,3); hold on; 
            a = shadedErrorBar(1:8,mean(dslatency(:,echoselectind),2),std(dslatency(:,echoselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.2,0.2,0.8],'LineWidth',2},1);
            b = shadedErrorBar(1:8,mean(dslatency(:,socialselectind),2),std(dslatency(:,socialselectind)')./sqrt(numel(echoselectind))*2,{'Color',[0.8,0.2,0.2],'LineWidth',2},1);
            c = shadedErrorBar(1:8,mean(dslatency(:,notselectiveind),2),std(dslatency(:,notselectiveind)')./sqrt(numel(echoselectind))*2,{'Color',[0,0,0],'LineWidth',2},1);
            title('Latency DS')
            ylabel('Median Latency [s]')
            xlabel('Rate [kHz/ms]')
            xticks(1:8)
            xticklabels({Sweepstruct(21).uniqueratekhz}); xlim([1,8])
            legend([a.mainLine,b.mainLine,c.mainLine],{'Echo','Social','Not selective'}); legend boxoff
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_Latency.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_Latency.svg'],'-dsvg')
      
            % Peak amplitude generally for each category per stim type
            [tuninghistecho,tuninghistechoind] = histcounts(mean(mediandftuning(echoselectind,:),2),10);
            [tuninghistsocial,tuninghistsocialind] = histcounts(mean(mediandftuning(socialselectind,:),2),10);
            [tuninghistother,tuninghistotherind] = histcounts(mean(mediandftuning(notselectiveind,:),2),10);
            [whhistecho,whhistechoind] = histcounts(mean(whpeak(:,echoselectind)),10);
            [whhistsocial,whhistsocialind] = histcounts(mean(whpeak(:,socialselectind)),10);
            [whhistother,whistotherind] = histcounts(mean(whpeak(:,notselectiveind)),10);
            [ushistecho,ushistechoind] = histcounts(mean(uppeak(:,echoselectind)),10);
            [ushistsocial,ushistsocialind] = histcounts(mean(uppeak(:,socialselectind)),10);
            [ushistother,usistotherind] = histcounts(mean(uppeak(:,notselectiveind)),10);
            [dshistecho,dshistechoind] = histcounts(mean(dspeak(:,echoselectind)),10);
            [dshistsocial,dshistsocialind] = histcounts(mean(dspeak(:,socialselectind)),10);
            [dshistother,dsistotherind] = histcounts(mean(dspeak(:,notselectiveind)),10);
 
            h = figure; hold on; subplot(1,4,1); hold on;
            plot(tuninghistechoind(1:end-1)+diff(tuninghistechoind(1:2))/2,tuninghistecho/sum(tuninghistecho),'Color',[0.2,0.2,0.8],'LineWidth',2)
            plot(tuninghistsocialind(1:end-1)+diff(tuninghistsocialind(1:2))/2,tuninghistsocial/sum(tuninghistsocial),'Color',[0.8,0.2,0.2],'LineWidth',2)
            plot(tuninghistotherind(1:end-1)+diff(tuninghistotherind(1:2))/2,tuninghistother/sum(tuninghistother),'Color',[0,0,0],'LineWidth',2)
%             legend('Echo','Social','Not selective'); legend boxoff
            ylabel('% cells'); xlabel('Mean peak amplitude'); 
            title('Tones'); ylim([0,0.4]); xlim([min(tuninghistechoind(1:end-1)+diff(tuninghistechoind(1:2))/2) max(tuninghistechoind(1:end-1)+diff(tuninghistechoind(1:2))/2)])
            subplot(1,4,2); hold on;
            plot(whhistechoind(1:end-1)+diff(whhistechoind(1:2))/2,whhistecho/sum(whhistecho),'Color',[0.2,0.2,0.8],'LineWidth',2)
            plot(whhistsocialind(1:end-1)+diff(whhistsocialind(1:2))/2,whhistsocial/sum(whhistsocial),'Color',[0.8,0.2,0.2],'LineWidth',2)
            plot(whistotherind(1:end-1)+diff(whistotherind(1:2))/2,whhistother/sum(whhistother),'Color',[0,0,0],'LineWidth',2)
%             legend('Echo','Social','Not selective'); legend boxoff
            ylabel('% cells'); xlabel('Mean peak amplitude'); 
            title('White noise');ylim([0,0.4]);  xlim([min(whhistechoind(1:end-1)+diff(whhistechoind(1:2))/2) max(whhistechoind(1:end-1)+diff(whhistechoind(1:2))/2)])
            subplot(1,4,3); hold on;
            plot(ushistechoind(1:end-1)+diff(ushistechoind(1:2))/2,ushistecho/sum(ushistecho),'Color',[0.2,0.2,0.8],'LineWidth',2)
            plot(ushistsocialind(1:end-1)+diff(ushistsocialind(1:2))/2,ushistsocial/sum(ushistsocial),'Color',[0.8,0.2,0.2],'LineWidth',2)
            plot(usistotherind(1:end-1)+diff(usistotherind(1:2))/2,ushistother/sum(ushistother),'Color',[0,0,0],'LineWidth',2)
%             legend('Echo','Social','Not selective'); legend boxoff
            ylabel('% cells'); xlabel('Mean peak amplitude'); xlim([min(ushistechoind(1:end-1)+diff(ushistechoind(1:2))/2) max(ushistechoind(1:end-1)+diff(ushistechoind(1:2))/2)])
            title('US');ylim([0,0.4])
            subplot(1,4,4); hold on;
            plot(dshistechoind(1:end-1)+diff(dshistechoind(1:2))/2,dshistecho/sum(dshistecho),'Color',[0.2,0.2,0.8],'LineWidth',2)
            plot(dshistsocialind(1:end-1)+diff(dshistsocialind(1:2))/2,dshistsocial/sum(dshistsocial),'Color',[0.8,0.2,0.2],'LineWidth',2)
            plot(dsistotherind(1:end-1)+diff(dsistotherind(1:2))/2,dshistother/sum(dshistother),'Color',[0,0,0],'LineWidth',2)
            legend('Echo','Social','Not selective'); legend boxoff
            ylabel('% cells'); xlabel('Mean peak amplitude'); 
            title('DS');ylim([0,0.4]); xlim([min(dshistechoind(1:end-1)+diff(dshistechoind(1:2))/2) max(dshistechoind(1:end-1)+diff(dshistechoind(1:2))/2)])
            saveas(h,[DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_amplitude.pdf'])
            print([DataDrive ':\Jenni\' animal '\figure\social\Population\Selectivity_amplitude.svg'],'-dsvg')
            

     
            % duration tuning
            
            
            % 
        case 'category'
            clear selectivityind selectivityind12 dfcat2 dfcat1
            close all
            
            responsewindow = 10:18;
            baselinewindow = 1:9;
            nbcells = size(datarois(fn).roiMatevokedstim,1);
            nbstim = size(datarois(fn).roiMatevokedstim,2);
            nbrep = size(datarois(fn).roiMatevokedstim,3);
            matstimevoked = datarois(fn).roiMatevokedstim(:,:,:,responsewindow);
            matbaseline = datarois(fn).roiMatevokedstim(:,:,:,baselinewindow);
            testmat = squeeze(median(datarois(fn).roiMatevokedstim(:,:,:,responsewindow),4));
            baselinemat = (squeeze(median(datarois(fn).roiMatevokedstim(:,:,:,baselinewindow),4)));
            meantestmat = mean(testmat,3);
            vecttest = reshape(meantestmat,[1,nbstim*nbcells]);
            teststruct.cells = reshape(repmat([1:nbcells],nbstim,1),[1,nbstim*nbcells]);%repmat([1:nbcells]',nbstim,1);%repmat(mat2cell(1:nbcells,1,ones(1,nbcells))',nbstim,1);
            teststruct.stim = repmat([1:nbstim]',nbcells,1)';
            %             baselineresp = squeeze(median(datarois(fn).roiMatevokedstim(:,:,:,baselineresponse),4));
            % ttest per stim
            % Do paired tests per stim per cells
            PeakAmp = squeeze(max(permute(matstimevoked,[4,1,2,3])));
            BaseAmp = squeeze(max(permute(matbaseline,[4,1,2,3])));
            MedianAmp = squeeze(median(permute(matstimevoked,[4,1,2,3])));
            MedianBase = squeeze(median(permute(matbaseline,[4,1,2,3])));
            AreaAmp = squeeze(trapz(permute(matstimevoked,[4,1,2,3])));%./length(responsewindow);
            AreaBase = squeeze(trapz(permute(matbaseline,[4,1,2,3])));%./length(responsewindow);
            ttestAlpha = 0.05;
            for cc = 1:nbcells
                peakAct = squeeze(PeakAmp(cc,:,:));
                baseAct = squeeze(BaseAmp(cc,:,:));
                medAct = squeeze(MedianAmp(cc,:,:));
                medBase = squeeze(MedianBase(cc,:,:));
                areaAct = squeeze(AreaAmp(cc,:,:));
                areaBase = squeeze(AreaBase(cc,:,:));
                for j = 1:nbstim
                    %                     try
                    [h,p] = ttest(peakAct(j,:),baseAct(j,:),'alpha',ttestAlpha,'tail','right');
                    ttestToneP(j,cc) = p;
                    ttestToneH(j,cc) = h;
                    [h,p] = ttest(medAct(j,:),medBase(j,:),'alpha',ttestAlpha,'tail','right');
                    medttestToneP(j,cc) = p;
                    medttestToneH(j,cc) = h;
                    [h,p] = ttest(areaAct(j,:),areaBase(j,:),'alpha',ttestAlpha,'tail','right');
                    areattestToneP(j,cc) = p;
                    areattestToneH(j,cc) = h;
                end
            end
            stimevokedcells = find(sum(ttestToneH)>0);
            notevokedcells = find(sum(ttestToneH)==0);
            prctevokedcells = (numel(stimevokedcells)/nbcells)*100;
            medstimevokedcells = find(sum(medttestToneH)>0);
            mednotevokedcells = find(sum(medttestToneH)==0);
            medprctevokedcells = (numel(medstimevokedcells)/nbcells)*100;
            areastimevokedcells = find(sum(areattestToneH)>0);
            areanotevokedcells = find(sum(areattestToneH)==0);
            areaprctevokedcells = (numel(areastimevokedcells)/nbcells)*100;
            
            % test for category analysis
            if nbstim == 12
                Category = [ones(1,5),ones(1,5)*2,ones(1,2)*3];
            elseif nbstim == 17
                Category = [ones(1,6),ones(1,6)*2,4,ones(1,4)*3];
            end
            % simple ttest for average per stim per cat for cat 1 & 2
            %             for cc = 1:nbcells
            %                 peakAct = reshape(MedianAmp(cc,find(Category==1),:),[1,numel(find(Category==1))*nbrep]);
            %                 baseAct = reshape(MedianAmp(cc,find(Category==2),:),[1,numel(find(Category==1))*nbrep]);
            %                 %                     try
            %                 [h,p] = ttest(peakAct,baseAct,'alpha',ttestAlpha,'tail','right');
            %                 catttestToneP(cc) = p;
            %                 catttestToneH(cc) = h;
            %             end
            %             selectivecells = find(catttestToneH);
            % do anova per cells over stims + baseline
            for cc = 1:nbcells
                groupNames = reshape(repmat([1:(nbstim*2)],nbrep,1),[1,(nbstim*2)*nbrep]);%reshape(repmat([1:(nbstim+1)],nbrep,1),[1,(nbstim+1)*nbrep]);
                mattemp = cat(1,squeeze(MedianAmp(cc,:,:)),squeeze(MedianBase(cc,:,:)));
                [p,h,stats] = anova1(reshape(mattemp',[(nbstim*2)*nbrep,1])',groupNames,'off');
                anovaPeak (1,cc) = p;
                [results,~,~,~] = multcompare(stats,'Display','off');
                for stimid = 1:nbstim
                    results_reduced(stimid) = results(intersect(find(round(results(:,1))==stimid),find(round(results(:,2))==(stimid+nbstim))),6);%results(find(results(:,1)==1),6); % select comparison with baseline
                end
                anovaSignifTone{cc} = (results_reduced<0.05);
                %                 anovaPeak (2,cc) = p<0.05 && sum(results<0.05)>0;
                %                 [results,~,~,~] = multcompare(stats,'Display','off');
                %                 results = results(results(:,2)==(nbstim+1),6);
                %                 anovaSignifTone(:,cc) = (results<0.05);
                %                 anovaPeak (2,cc) = p<0.05 && sum(results<0.05)>0;
                %                 [p,h,stats] = anova1(peakANOVACorr(:,:,i),groupNames,'off');
                %                 anovaPeakCorr (1,i) = p;
                %
                %                 [results,~,~,~] = multcompare(stats,'Display','off');
                %                 results = results(results(:,2)==(nTones+1),6);
                %                 anovaSignifToneCorr(:,i) = (results<0.05);
                %                 anovaPeakCorr (2,i) = p<0.05 && sum(results<0.05)>0;
                %
                % category test 1 vs 2 (FMBs vs Echo)
                peaktestcat{cc} = [reshape(PeakAmp(cc,find(Category==1),:),nbrep*max(find(Category==1)),1),reshape(PeakAmp(cc,find(Category==2),:),nbrep*max(find(Category==1)),1)];
                medtestcat{cc} = [reshape(MedianAmp(cc,find(Category==1),:),nbrep*max(find(Category==1)),1),reshape(MedianAmp(cc,find(Category==2),:),nbrep*max(find(Category==1)),1)];
                areatestcat{cc} = [reshape(AreaAmp(cc,find(Category==1),:),nbrep*max(find(Category==1)),1),reshape(AreaAmp(cc,find(Category==2),:),nbrep*max(find(Category==1)),1)];
            end
            % identify cvells that are category specific with anova1
            peakanovacells = find(cellfun(@(x) anova1(x,[],'off'),peaktestcat)<=0.05)
            medanovacells = find(cellfun(@(x) anova1(x,[],'off'),medtestcat)<=0.05)
            areaanovacells = find(cellfun(@(x) anova1(x,[],'off'),areatestcat)<=0.05)
            
            
            
            % not great
            % selectivity index
            
            %             for cc = 1:nbcells
            %                 mediancat1{cc} = [];
            %                 mediancat2{cc} = [];
            %                 mediancat3{cc} = [];
            %                 for j = find(Category==1)
            %                     mediancat1{cc} = [mediancat1{cc};squeeze(datarois(fn).roiMatevokedstim(cc,j,:,:))];
            %                 end
            %                 for j = find(Category==2)
            %                     mediancat2{cc} = [mediancat2{cc};squeeze(datarois(fn).roiMatevokedstim(cc,j,:,:))];
            %                 end
            %                 for j = find(Category==3)
            %                     mediancat3{cc} = [mediancat3{cc};squeeze(datarois(fn).roiMatevokedstim(cc,j,:,:))];
            %                 end
            %                 selectivityind12(cc,:) = (median(mediancat1{cc})-median(mediancat2{cc}))./(median(mediancat1{cc})+median(mediancat2{cc}));
            %             end
            %             % selectivity index sorted by min if negative and max if
            %             % positive
            % %                        negcells = find(median(selectivityind12(:,responsewindow),2)<=0);
            % %                        poscells = find(median(selectivityind12(:,responsewindow),2)>0);
            % %                        [selectvityvalneg,selectvityindneg] = sort(min(selectivityind12(negcells,responsewindow)'));
            %             %            [selectvityvalpos,selectvityindpos] = sort(max(selectivityind12(poscells,responsewindow)'));
            %             medselect = median(selectivityind12(:,responsewindow)');
            %             [selectivityval,selectivityind] =  sort(medselect);
            %             selectivecat1 = selectivityind(selectivityval>=prctile(selectivityval,95));
            %             selectivecat2 = selectivityind(selectivityval<=prctile(selectivityval,5));
            %              % tuning curve analysis
            %             figure(); hold on;
            %             subplot(3,1,[1,2]); hold on; imagesc(xvalues,1:nbcells,selectivityind12(selectivityind,:)); caxis([prctile(selectivityval,5),prctile(selectivityval,95)]); colormap(flipud(redblue))
            %             axis tight; colorbar
            %             subplot(3,1,3); hold on;  shadedErrorBar(xvalues,median(selectivityind12(selectivecat1,:)),std(selectivityind12(selectivecat1,:)),{'Color',stimcolormap(1,:),'LineWidth',1.5},1);
            %             subplot(3,1,3); hold on;  shadedErrorBar(xvalues,median(selectivityind12(selectivecat2,:)),std(selectivityind12(selectivecat1,:)),{'Color',stimcolormap(7,:),'LineWidth',1.5},1);
            %             ylim([-0.1,0.1]); ylabel('Median selectivity top 5%');
            %             xlim([xvalues(1),xvalues(end)]); xlabel('Time from stim onset [s]')
            %             print(['N:\Jenni\' animal '\figure\social\Population\site' num2str(datarois(fn).site) '\site' num2str(datarois(fn).site) '_HeatMapSelectivity.pdf'],'-dpdf','-fillpage')
            %            close
            %             % tuning curve sorted by selectivity index
            %             tuningpeaknorm = mean(PeakAmp,3)./max(mean(PeakAmp,3)')';
            %             tuningmednorm = mean(MedianAmp,3)./max(mean(MedianAmp,3)')';
            %             tuningareanorm = mean(AreaAmp,3)./max(mean(AreaAmp,3)')';
            %             tuningpeak = mean(PeakAmp,3);
            %             tuningmed = mean(MedianAmp,3);
            %             tuningarea = mean(AreaAmp,3);
            %
            %             figure(); hold on; subplot(1,3,1); hold on; imagesc(tuningpeaknorm(selectivityind,:)); caxis([0.8,1]);axis tight
            %             title('peak')
            %             subplot(1,3,2); hold on; imagesc(tuningmednorm(selectivityind,:)); caxis([0.8,1]);axis tight
            %             title('median')
            %             subplot(1,3,3); hold on; imagesc(tuningareanorm(selectivityind,:)); caxis([0.8,1]);axis tight
            %             title('area')
            %             print(['N:\Jenni\' animal '\figure\social\Population\site' num2str(datarois(fn).site) '\site' num2str(datarois(fn).site) '_HeatMapSocialTuning.pdf'],'-dpdf','-fillpage')
            %            close
            %             % map activity per stim
            %             col = [0.5,0.5,0.5];
            %             meanimg = datarois(fn).meanImg;
            % %             bla = colormap(redblue(1000));
            %             colormapselect = colormap(flipud(redblue(100)));
            %             selectcorrespondance = linspace(-1,1,100);
            %             colorredwhite = flipud([ones(100,1),(log10([1.01:0.1:11])/log10(11))',(log10([1.01:0.1:11])/log10(11))']);%([0:0.01:1-0.01].^1.2)',([0:0.01:1-0.01].^1.05)'];
            % %             [histpeak,histpeakvalue] = hist(reshape(tuningpeak,1,size(tuningpeak,1)*size(tuningpeak,2)),100);
            % %             [histmed,histmedvalue] = hist(reshape(tuningmed,1,size(tuningpeak,1)*size(tuningpeak,2)),100);
            % %             [histarea,histareavalue] = hist(reshape(tuningarea,1,size(tuningarea,1)*size(tuningarea,2)),100);
            %
            %             tuningpeakmin = (mean(PeakAmp,3)-min(mean(PeakAmp,3)')')./((max(mean(PeakAmp,3)')')-min(mean(PeakAmp,3)')');
            %             difftuningpeak = mean(PeakAmp(:,find(Category==1),:),3)-mean(PeakAmp(:,find(Category==2),:),3);
            %             difftuningpeaknorm = (difftuningpeak-min(difftuningpeak))./(max(difftuningpeak)-min(difftuningpeak))
            
            % caxis max 1.1 median value from 0.95 min
            
            colormapselect = colormap(flipud(redblue(100)));
            selectcorrespondance = linspace(-1,1,100);
            colorredwhite = flipud([ones(100,1),(log10([1.01:0.1:11])/log10(11))',(log10([1.01:0.1:11])/log10(11))']);%([0:0.01:1-0.01].^1.2)',([0:0.01:1-0.01].^1.05)'];
            
            tuningpeak = mean(PeakAmp,3);
            tuningmed = mean(MedianAmp,3);
            tuningarea = mean(AreaAmp,3);
            minb = 0.9;
            maxb = 1.2;
            vectmed = reshape(tuningmed,1,size(tuningmed,1)*size(tuningmed,2));
            boundedvect = (vectmed-minb)./(maxb-minb);
            correspondancevectmed = nan(1,size(vectmed,2));
            correspondancevectmed = ceil(boundedvect*100);
            correspondancevectmed(find(correspondancevectmed>=100)) = length(colorredwhite);
            correspondancevectmed(find(correspondancevectmed<=0)) = 1;
            correspondancematmed = reshape(correspondancevectmed,nbcells,nbstim);
            dfmed = tuningmed-1; % in the sojner paper they remove the baseline, 1 is an approximation
            dfcat1med = mat2cell(dfmed(:,find(Category==1)),ones(1,nbcells),max(find(Category==1)));
            dfcat2med = mat2cell(dfmed(:,find(Category==2)),ones(1,nbcells),max(find(Category==1)));
            dfpeak = tuningpeak-1; % in the sojner paper they remove the baseline, 1 is an approximation
            dfcat1peak = mat2cell(dfpeak(:,find(Category==1)),ones(1,nbcells),max(find(Category==1)));
            dfcat2peak = mat2cell(dfpeak(:,find(Category==2)),ones(1,nbcells),max(find(Category==1)));
            dfarea = tuningarea-1; % in the sojner paper they remove the baseline, 1 is an approximation
            dfcat1area = mat2cell(dfarea(:,find(Category==1)),ones(1,nbcells),max(find(Category==1)));
            dfcat2area = mat2cell(dfarea(:,find(Category==2)),ones(1,nbcells),max(find(Category==1)));
            %             for stim = 1:max(max(find(Category==1)))
            selectivityindexperstimmed = cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),dfcat1med,dfcat2med,'UniformOutput',0);
            selectivityindexperstimpeak = cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),dfcat1peak,dfcat2peak,'UniformOutput',0);
            selectivityindexperstimarea = cellfun(@(x,y) (abs(x)-abs(y))./(abs(x)+abs(y)),dfcat1area,dfcat2area,'UniformOutput',0);
            
            %         save(['N:/Jenni/' animal '/roistructure/socialroisdata.mat'],'datarois')
            
            
            
        case 'plotselectivitymaps'
            if strcmp(animal,'RoofBuddy1')
                datarois = datarois(2:end);
            end
            clear selectcorrespondance colormapselect selectcorrespondance significantcells vectcorrespondanceselect
            
            % ttest category based
            
            
            %             end
            
            % plot histograms of selectivity: mean selectivity per stim
            % based on median response window
            % test for category analysis
            if size(datarois(fn).selectivityindexperstimmed{1},2)==5
                Category = [ones(1,5),ones(1,5)*2,ones(1,2)*3];
            elseif size(datarois(fn).selectivityindexperstimmed{1},2)==6
                Category = [ones(1,6),ones(1,6)*2,4,ones(1,4)*3];
            end
            meanselect = cellfun(@mean,datarois(fn).selectivityindexperstimmed);
            meanmeanselectshuffle = mean(cell2mat(cellfun(@(x) mean(x,1),datarois(fn).selectivityindexperstimmedshuffle,'UniformOutput',0)),2);
            boundarycells = cell2mat(cellfun(@(x) prctile(mean(x),[5,95]),datarois(fn).selectivityindexperstimmedshuffle,'UniformOutput',0));
            significantcells = [find(meanselect<(boundarycells(:,1)));find(meanselect>(boundarycells(:,2)))];
            meanselectshuffle = reshape(cell2mat(cellfun(@(x) mean(x,1),datarois(fn).selectivityindexperstimmedshuffle,'UniformOutput',0)),[size(datarois(fn).selectivityindexperstimmedshuffle,1)*size(datarois(fn).selectivityindexperstimmedshuffle{1},2),1]);
            [sortedmeanselect,sortedmeanselectid] = sort(meanselect);
            [sortedmeanselectshuffle,sortedmeanselectidshuffle] = sort(meanmeanselectshuffle);
            h= figure(); hold on; subplot(1,2,1); hold on;
            [histselectvalues,histselectbins]=hist(meanselect,20);
            [histselectvaluesshuffle,histselectbinsshuffle]=hist(meanselectshuffle,20);
            plot([histselectbins(1),histselectbins,histselectbins(end)],[0, histselectvalues./max(histselectvalues),0],'k','LineWidth',1.5);
            plot([histselectbinsshuffle(1),histselectbinsshuffle,histselectbinsshuffle(end)],[0, histselectvaluesshuffle./max(histselectvaluesshuffle),0],'Color',[0.5,0.5,0.5],'LineWidth',1.5);
            
            plot([0,0],[0,max(histselectvalues)],'--k'); xlim([min(histselectbinsshuffle),max(histselectbinsshuffle)]);
            ylim([0,1]);%max(histselectvaluesshuffle)]);
            subplot(1,2,2); hold on;
            plot(1:numel(find(sortedmeanselect<0)),sortedmeanselect(find(sortedmeanselect<0)),'.','MarkerSize',10,'MarkerFaceColor',[0.3,0.2,0.8]); xlim([1,length(sortedmeanselect)])
            plot(numel(find(sortedmeanselect<0))+1:numel(find(sortedmeanselect>=0))+numel(find(sortedmeanselect<0)),sortedmeanselect(find(sortedmeanselect>=0)),'.','MarkerSize',10,'MarkerFaceColor',[0.8,0.3,0.2]); xlim([1,length(sortedmeanselect)])
            plot(sortedmeanselectshuffle,'.','MarkerSize',10,'MarkerFaceColor',[0.7,0.7,0.7],'MarkerEdgeColor',[0.7,0.7,0.7]);
            
            
            
            % plot maps
            boundaryhist = prctile(meanselectshuffle,[5,95]);
            minselect = -0.3;%boundaryhist(1);
            maxselect = 0.3;%boundaryhist(2);
            colormapselect = colormap((redblue(100)));
            selectcorrespondance = linspace(-1,1,100);
            selectcorrespondance2 = linspace(minselect,maxselect,100);
            colorredwhite = flipud([ones(100,1),(log10([1.01:0.1:11])/log10(11))',(log10([1.01:0.1:11])/log10(11))']);%([0:0.01:1-0.01].^1.2)',([0:0.01:1-0.01].^1.05)'];
            %             minb = 0.9;
            %             maxb = 1.2;
            %             vectmed = reshape(tuningmed,1,size(tuningmed,1)*size(tuningmed,2));
            %             boundedvect = (vectmed-minb)./(maxb-minb);
            %
            %             correspondancevectmed = nan(1,size(vectmed,2));
            %             correspondancevectmed = ceil(boundedvect*100);
            %             correspondancevectmed(find(correspondancevectmed>=100)) = length(colorredwhite);
            %             correspondancevectmed(find(correspondancevectmed<=0)) = 1;
            
            
            posselect = find(meanselect>0);
            negselect = find(meanselect<=0);
            restselectneg = find(meanselect>minselect & meanselect<0);
            restselectpos = find(meanselect<maxselect & meanselect>=0);
            vectcorrespondanceselect(posselect(find(meanselect(posselect)>=maxselect))) = maxselect;
            vectcorrespondanceselect(negselect(find(meanselect(negselect)<=minselect))) = minselect;
            for i = 1:length(restselectneg)
                vectcorrespondanceselect(restselectneg(i)) = selectcorrespondance2(max(find(meanselect(restselectneg(i))>=selectcorrespondance2)));
            end
            for i = 1:length(restselectpos)
                vectcorrespondanceselect(restselectpos(i)) = selectcorrespondance2(max(find(meanselect(restselectpos(i))>=selectcorrespondance2)));
            end
            for ii = 1:length(vectcorrespondanceselect)
                selectcorrespondance(ii) = find(selectcorrespondance2==vectcorrespondanceselect(ii));
            end
            meanimg = datarois(fn).meanImg;
            %             h=figure(); hold on;
            %             for stim = 1:max(find(Category==2))
            %                 subplot(3,max(find(Category==1)),stim); hold on;
            %                 imagesc(1:size(meanimg,2),1:size(meanimg,1),imadjust(int16(meanimg))); colormap gray
            %                 %             cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',col),datarois(fn).roisCoord{1})
            %                 cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',colorredwhite(y,:)),datarois(fn).roisCoord{1},mat2cell(correspondancematmed(:,stim)',1,ones(1,length(correspondancematmed(:,stim)))))
            %                 %cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',bla(round(y*999)+1,:)),datarois(fn).roisCoord{1},mat2cell(tuningpeakmin(:,stim)',1,ones(1,length(tuningpeakmin(:,stim)))))
            %                 xlim([1,size(meanimg,2)]); ylim([1,size(meanimg,1)])
            %                 set(gca,'YDir','reverse')
            %                 set(gca,'xtick',[]); set(gca,'ytick',[]); axis square
            %                 if stim > max(find(Category==1))
            %                     subplot(3,max(find(Category==1)),max(find(Category==2)+(stim-max(find(Category==1))))); hold on;
            %                     imagesc(1:size(meanimg,2),1:size(meanimg,1),imadjust(int16(meanimg))); colormap gray
            %                     %             cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',col),datarois(fn).roisCoord{1})
            %                     cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',colormapselect(find(selectcorrespondance>=y(stim-max(find(Category==1))),1),:)),datarois(fn).roisCoord{1},selectivityindexperstim');
            %                     % cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',bli(round(y*999)+1,:)),datarois(fn).roisCoord{1},mat2cell(difftuningpeaknorm(:,stim-max(find(Category==1)))',1,ones(1,length(difftuningpeaknorm(:,stim-max(find(Category==1)))))));
            %
            %                     xlim([1,size(meanimg,2)]); ylim([1,size(meanimg,1)])
            %                     set(gca,'YDir','reverse')
            %                     set(gca,'xtick',[]); set(gca,'ytick',[]); axis square
            %                 end
            %             end
            %             print(['N:\Jenni\' animal '\figure\social\Population\site' num2str(datarois(fn).site) '\site' num2str(datarois(fn).site) '_MapPopulationActivity.pdf'],'-dpdf','-fillpage')
            % map selectivity index
            close
            %             medselectnorm = (medselect-min(medselect))./(max(medselect)-min(medselect));
            h=figure(); hold on;
            %             imagesc(1:size(meanimg,2),1:size(meanimg,1),imadjust(int16(meanimg))); colormap gray
            %             cellfun(@(x) patch('Faces',1:size(x,1),'Vertices',x,'EdgeColor','none','FaceColor',{}),datarois(fn).roisCoord{1})
            cellfun(@(x,y) patch('Faces',1:size(x,1),'Vertices',[(750-x(:,1)),x(:,2)],'EdgeColor','none','FaceColor',colormapselect(y,:)),datarois(fn).roisCoord{1},num2cell(selectcorrespondance))
            centroidcells = cellfun(@(x) median(x),datarois(fn).roisCoord{1},'UniformOutput',0)
            xlim([1,size(meanimg,2)]); ylim([1,size(meanimg,1)])
            set(gca,'xtick',[]); set(gca,'ytick',[]);
            print([animaldrive ':\Jenni\' animal '\figure\social\Population\site' num2str(datarois(fn).site) '\site' num2str(datarois(fn).site) '_MapPopulationSelectivity.pdf'],'-dpdf','-fillpage')
            close
            
            %             for cc = 1:length(centroidcells)
            %                 for cci = length(centroidcells):-1:1
            %                     dist(cc,cci) = pdist([centroidcells{cc};centroidcells{cci}],'euclidean');
            %                     if meanselect(cc)>0.4 || meanselect(cc)<-0.4 && meanselect(cci)>0.4 || meanselect(cci)<-0.4 && dist(cc,cci)>0
            %                         pairselect(cc,cci) = meanselect(cc)+meanselect(cci);
            %                     else
            %                         pairselect(cc,cci) = nan;
            %                     end
            %                 end
            %             end
            %             % plot selectivity index per dist
            %              figure(); hold on; plot(reshape(pairselect,1,size(dist,1)*size(dist,1)),reshape(dist,1,size(dist,1)*size(dist,1)),'.')
            %
            %
            %             p = anovan(vecttest,{teststruct.stim,teststruct.cells});
            %             [a,b,c] = anova2(meantestmat);
            %             t = table(teststruct.cells,meantestmat(:,1),meantestmat(:,2),meantestmat(:,3),meantestmat(:,4),meantestmat(:,5),meantestmat(:,6),meantestmat(:,7),meantestmat(:,8),meantestmat(:,9),meantestmat(:,10),meantestmat(:,11),meantestmat(:,12),'VariableNames',{'cells','meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10','meas11','meas12'});
            %             Meas = table([1:12]','VariableNames',{'Measurements'});
            %             rm = fitrm(t,'meas1-meas12','WithinDesign',Meas);
            
            %             teststruct.stim = repmat([1:nbstim]',nbcells,1); %mat2cell(reshape(repmat(1:nbstim,nbcells,1),nbstim*nbcells,1),ones(1,nbstim*nbcells),1);
            
            %             teststruct.cells = reshape(repmat(1:nbcells,nbstim,1),nbstim*nbcells,1);%repmat([1:nbcells]',nbstim,1);%repmat(mat2cell(1:nbcells,1,ones(1,nbcells))',nbstim,1);
            %             teststruct.stim = repmat([1:nbstim]',nbcells,1); %mat2cell(reshape(repmat(1:nbstim,nbcells,1),nbstim*nbcells,1),ones(1,nbstim*nbcells),1);
            % %             if nbstim == 12
            %                 teststruct.category = repmat([ones(1,5),ones(1,5)*2,ones(1,2)*3]',nbcells,1);
            %             else
            %                 teststruct.category = repmat([ones(1,6),ones(1,6)*2,4,ones(1,4)*3]',nbcells,1);
            %             end
            %             teststruct.meas = [];
            %              for cc = 1:nbcells
            %                  for stimid = 1:nbstim
            %                      teststruct.meas = [teststruct.meas;squeeze(testmat(cc,stimid,:))'];
            %                  end
            %              end
            %             % repated measure of anova within subject for cells
            %             t = table(teststruct.cells,teststruct.stim,teststruct.category,teststruct.meas(:,1),teststruct.meas(:,2),teststruct.meas(:,3),teststruct.meas(:,4),teststruct.meas(:,5),teststruct.meas(:,6),teststruct.meas(:,7),teststruct.meas(:,8),teststruct.meas(:,9),teststruct.meas(:,10),'VariableNames',{'cells','stim','category','meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10'});
            % %             anovan(teststruct.meas)
            %              Meas = table([3:12]','VariableNames',{'Measurements'});
            % %             stim = t(:,2);
            %             withinDesign = table([1:12]','VariableNames',{'Condition'});
            %             withinDesign.Stim = categorical(teststruct.stim);
            %             % create the repeated measures model and do the anova
            %                Meas = table([1:10]','VariableNames',{'Measurements'});
            %                rm = fitrm(t,'meas1-meas10 ~ stim','WithinDesign',Meas);
            % %             rm = fitrm(t,'meas1-meas10~stim','WithinDesign',Meas,'WithinModel','stim');
            %             [ranovatable,A,B] = ranova(rm);
            %             % anova1 stim per cell
            
        case 'psth'
            close all
            graycmapinv = flipud(colormap(gray));
            responsewindow = 11:19; % (baseline from 1:10, 10 is stim onset)
            figure(fn); hold on;
            [sortedcellsval,sortedcells] = sort(max(squeeze(mean(squeeze(mean(datarois(fn).roiMatevokedstim(:,:,:,responsewindow),3)),2))'));
            subplot(2,2,1); hold on; imagesc(xvalues,1:size(sortedcells,2),cell2mat(cellfun(@(x) median(x./median(x(:,1:10),2)),datarois(fn).roiMatevoked(sortedcells),'UniformOutput',false)')); caxis([0.98,1.1]); axis tight; 
            colormap(gca,'gray'); colorbar
            title('PSTH (median evoked response)')
            xlabel('Time fron stim onset [s]')
            for stimid = 1:length(datarois(fn).correctedstim)
                subplot(2,2,3); hold on; plot(xvalues,mean(squeeze(median(datarois(fn).roiMatevokedstim(sortedcells,datarois(fn).correctedstim(stimid),:,:),3))),'Color',stimcolormap(datarois(fn).correctedstim(stimid),:),'LineWidth',1.5); axis tight
            end
            ylim([0.95,1.25])
            legend({stimwavname(datarois(fn).correctedstim).name}); legend boxoff
            title('Population PSTH')
            xlabel('Time fron stim onset [s]')
            % plot onset stim and end of stimulus % stim onset line and
            % shaded area for stim length (light grey)
            plot([xvalues(responsewindow(1)),xvalues(responsewindow(1))],[0.95,1.25],'--k');
%             plot([xvalues(responsewindow(1)),xvalues(responsewindow(1)),xvalues(find(round((xvalues-Lengthstim)*100)==0)),xvalues(find(round((xvalues-Lengthstim)*100)==0))],[0.95 1.25 0.95 1.25],'Color',[.8,0.8,0.8])
            plot([xvalues(find(round((xvalues-Lengthstim)*100)==0)),xvalues(find(round((xvalues-Lengthstim)*100)==0))],[0.95,1.25],'--k');
            
            
            % plot 'tuning curve per cell' max value in range 10 to 24 frame
            stimtuningcurve = squeeze(max(permute(squeeze(median(datarois(fn).roiMatevokedstim(:,:,:,responsewindow),3)),[3,2,1])));%./max(squeeze(max(permute(squeeze(median(datarois(fn).roiMatevokedstim(:,:,:,responsewindow),3)),[3,2,1]))));
            [~,sortedstimtuning] = max(stimtuningcurve);
            [~,bla] = sort(sortedstimtuning);
            figure(fn); hold on;
            subplot(2,2,2);hold on; imagesc(1:length(datarois(fn).correctedstim),1:size(stimtuningcurve,2),stimtuningcurve(:,bla)'); axis tight; caxis([0.5,1]); axis tight; 
            colormap(gca,graycmapinv); colorbar
            title(['Peak Amplitude (df/f)'])
            xticks([1:length(datarois(fn).correctedstim)])
            xticklabels(cellfun(@(x) x(1:end-4),{stimwavname(datarois(fn).correctedstim).name},'UniformOutput',0))
            xtickangle(45)
            ylabel('cell nb')
            
            % plot latency per stim (max value for each stim)    
            for cc = 1:size(stimtuningcurve,1)
                for stimid = 1:length(datarois(fn).correctedstim)
                    stimlatency(cc,stimid) = find(squeeze(median(datarois(fn).roiMatevokedstim(cc,stimid,:,responsewindow),3))==stimtuningcurve(cc,stimid))/datarois(fn).framerate;
                end
            end
            figure(fn); hold on;
            subplot(2,2,4);hold on; imagesc(1:length(datarois(fn).correctedstim),size(stimtuningcurve,1),stimlatency(bla,:)); axis tight; caxis([1,11]./datarois(fn).framerate); axis tight; colormap(redblue); colorbar
            title('Latency for peak per stim')
            xticks([1:length(datarois(fn).correctedstim)])
            xticklabels(cellfun(@(x) x(1:end-4),{stimwavname(datarois(fn).correctedstim).name},'UniformOutput',0))
            xtickangle(45)
            print([Datadrive ':\Jenni\' animal '\figure\social\Population\site' num2str(datarois(fn).site) '\site' num2str(datarois(fn).site) '_SummaryEvokedResponse.pdf'],'-dpdf','-fillpage')
            
            bandwidthstim = round((sum((stimtuningcurve-min(stimtuningcurve')')./(max(stimtuningcurve')'-min(stimtuningcurve')'),2)/length(datarois(fn).correctedstim))*10)/10;
            
            %plot PSTHs for each cells (10 cells at a time)
            nbfig= ceil(length(sortedcells)/10);
            xvaluesstim = 1:(length(xvalues)-20)*(length(datarois(fn).correctedstim));%1:length([5:9,responsewindow])*(length(datarois(fn).correctedstim)+1);
            % get max or mean selectivity per cell
            %             maxpeaksel = cellfun(@max,datarois(fn).selectivityindexperstimpeak);
            %             maxmedsel = cellfun(@max,datarois(fn).selectivityindexperstimmed);
            %             maxareasel = cellfun(@max,datarois(fn).selectivityindexperstimarea);
            meanpeaksel = cellfun(@mean,datarois(fn).selectivityindexperstimpeak);
            meanmedsel = cellfun(@mean,datarois(fn).selectivityindexperstimmed);
            meanareasel = cellfun(@mean,datarois(fn).selectivityindexperstimarea);
            
            
            for nbfigind = 1:nbfig
                close all
                h = figure(fn+2); hold on;
                if nbfigind<nbfig
                    for cc = ((nbfigind-1)*10)+1:((nbfigind-1)*10)+10
                        for stimid = 1:length(datarois(fn).correctedstim)
                            subplot(5,2,cc-((nbfigind-1)*10)); hold on; plot(xvaluesstim(((length(xvalues)-20)*(stimid-1))+1):xvaluesstim(((length(xvalues)-20)*(stimid-1))+length(xvalues)-20),[smooth(squeeze(median(datarois(fn).roiMatevokedstim(cc,stimid,:,1:end-23),3)),3);nan(3,1)],'LineWidth',1.5,'Color',stimcolormap(datarois(fn).correctedstim(stimid),:));	%plot(xvaluesstim((length([5:9,responsewindow])*stimid-1)+1):xvaluesstim((length([5:9,responsewindow])*stimid-1)+16),smooth(squeeze(median(datarois(fn).roiMatevokedstim(cc,stimid,:,[5:9,responsewindow]),3)),3),'LineWidth',1.5,'Color',stimcolormap(datarois(fn).correctedstim(stimid),:)); %axis tight, caxis([0.8,1.2]); colormap(redblue)
                            plot([xvaluesstim(((length(xvalues)-20)*(stimid-1))+11) xvaluesstim(((length(xvalues)-20)*(stimid-1))+11)],[0.8, 1.3],'LineWidth',1,'Color',[0.5,0.5,0.5]); ylim([0.8 1.3])
                            %                 xlabel('Time from onset {s}'); %xlim([xvalues(5),xvalues(responsewindow(end))])
                            %                 if cc>((nbfigind-1)*10)+1
                            %                     yticks([])
                            %                 end
                            xticks([]);
                            if stimid == length(datarois(fn).correctedstim)
                                % add xaxis limits
                                xlim([xvaluesstim(1) xvaluesstim(end)]);
                            end
                        end
                        title({['Rois nb: ' num2str(cc)];['mean SI: med:' num2str(round(meanmedsel(cc)*100)/100) ' peak:' num2str(round(meanpeaksel(cc)*100)/100) ' area:' num2str(round(meanareasel(cc)*100)/100)]})%['max SI: med:' num2str(maxmedsel(cc)) ' peak:' num2str(maxpeaksel(cc)) ' area:' num2str(maxareasel(cc))];
                    end
                elseif nbfigind==nbfig
                    for cc = ((nbfigind-1)*10)+1:length(sortedcells)
                        for stimid = 1:length(datarois(fn).correctedstim)
                            subplot(5,2,cc-((nbfigind-1)*10)); hold on; plot(xvaluesstim(((length(xvalues)-20)*(stimid-1))+1):xvaluesstim(((length(xvalues)-20)*(stimid-1))+length(xvalues)-20),[smooth(squeeze(median(datarois(fn).roiMatevokedstim(cc,stimid,:,1:end-23),3)),3);nan(3,1)],'LineWidth',1.5,'Color',stimcolormap(datarois(fn).correctedstim(stimid),:));	%plot(xvaluesstim((length([5:9,responsewindow])*stimid-1)+1):xvaluesstim((length([5:9,responsewindow])*stimid-1)+16),smooth(squeeze(median(datarois(fn).roiMatevokedstim(cc,stimid,:,[5:9,responsewindow]),3)),3),'LineWidth',1.5,'Color',stimcolormap(datarois(fn).correctedstim(stimid),:)); %axis tight, caxis([0.8,1.2]); colormap(redblue)
                            plot([xvaluesstim(((length(xvalues)-20)*(stimid-1))+11) xvaluesstim(((length(xvalues)-20)*(stimid-1))+11)],[0.8, 1.3],'LineWidth',1,'Color',[0.5,0.5,0.5]); ylim([0.8 1.3])
                            %                 xlabel('Time from onset {s}'); %xlim([xvalues(5),xvalues(responsewindow(end))])
                            title({['Rois nb: ' num2str(cc)];['mean SI: med:' num2str(round(meanmedsel(cc)*100)/100) ' peak:' num2str(round(meanpeaksel(cc)*100)/100) ' area:' num2str(round(meanareasel(cc)*100)/100)]})%['max SI: med:' num2str(maxmedsel(cc)) ' peak:' num2str(maxpeaksel(cc)) ' area:' num2str(maxareasel(cc))];
                            %                 if cc>((nbfigind-1)*10)+1
                            %                     yticks([])
                            %                 end
                            xticks([]);
                            if stimid == length(datarois(fn).correctedstim)
                                % add xaxis limits
                                xlim([xvaluesstim(1) xvaluesstim(end)]);
                            end
                        end
                    end
                end
                suptitle(['Site: ' num2str(num2str(datarois(fn).site)), ' ind:' num2str(nbfigind)])
                print([DataDrive ':\Jenni\' animal '\figure\social\Cells\site' num2str(datarois(fn).site) '_StimPSTHs' num2str(nbfigind) '.pdf'],'-dpdf','-fillpage')
            end
            
            %     % check if cells have some specific stim they like % npo category
            %     %create mat pca
            %     if ~isempty(intersect(stimselection,analysefilelistind(fn)))
            %         numvoc = 2;
            %     else
            %         numvoc = 4;
            %     end
            %     for cc = 1:size(roiMatevoked,1)
            %         matpca(cc,:) = reshape(squeeze(median(roiMatevokedstim(cc,1:end-numvoc,:,:),3)),[1,(size(roiMatevokedstim,2)-numvoc)*size(roiMatevokedstim,4)]);
            %     end
            %     %     [a,b,c,explained]=pca(matpca');
            %     %     for pcind = 1:3
            %     %         for  ii = 1:length(uniquetone)
            %     %             for tt = 1:51
            %     %             projectionspc{pcind,ii}(:,tt) = dot(a(:,pcind)',(squeeze(median(roiMatevokedstim(:,ii,:,tt),3))));
            %     %             end
            %     %         end
            %     %     end
            %     %
            %     %     figure(3); hold on;
            %     %
            %     %     for  ii = 1%:length(uniquetone)
            %     %                     plot3(projectionspc{1,ii},projectionspc{2,ii},projectionspc{3,ii},'Color',colorcodeid(ii,:));
            %     %         for   tt = 1:22
            %     %             scatter3(projectionspc{1,ii}(:,tt),projectionspc{2,ii}(:,tt),projectionspc{3,ii}(:,tt),'Color',colorcodeid(ii,:),'MarkerSize',5*(tt/20));
            %     %         end
            %     %     end
            %     [coeff,score,latent,~,explained,~] = pca(matpca');
            %
            %     %%
            %     %            Xcentered = score*coeff';
            %     plot_range = 8:14;%51;
            %
            %     % build a colormap for each stim/category
            %
            %     % for stimind = 14:17
            %     % stimcolormap(14:17,:) = [0.9,0.2,0.2];
            %     % end
            %
            %     % else
            %
            %     % end
            %
            %
            %     figure(); hold on;
            %     for ii=1:length(correctedstim)-numvoc
            %         scorek{ii} = score(size(score,1)*(ii-1)/(length(correctedstim)-numvoc)+1:size(score,1)*ii/(length(correctedstim)-numvoc),1:3);
            %         %     if correctedstim(ii)<=6
            %         plot3(scorek{ii}(plot_range,1),scorek{ii}(plot_range,2),scorek{ii}(plot_range,3),'Color',stimcolormap(correctedstim(ii),:),'LineWidth',2); hold on;
            %         %         fnplt(cscvn(scorek{ii}(plot_range,:)'),1,'b'); hold on;
            %         %             plot3(scorek{ii}(plot_range,1),scorek{ii}(plot_range,2),scorek{ii}(plot_range,3),'b')
            %         %     elseif correctedstim(ii)>6 && correctedstim(ii)<13
            %         %         fnplt(cscvn(scorek{ii}(plot_range,:)'),1,'r'); hold on
            %         %     elseif correctedstim(ii) == 13
            %         %         fnplt(cscvn(scorek{ii}(plot_range,:)'),1,'g'); hold on
            %         %     end
            %     end
            %     set(gca,'ticklength',[0 0])
            %     set(gca,'XTickLabel',''); set(gca,'YTickLabel',''); set(gca,'ZTickLabel','');
            %     set(gca,'ticklength',[0 0])
            %     %     figure(1); hold on;
            %     % %     for ii=1:length(uniquetone)
            %     % %         for tt = 8:34
            %     % %             plot3(scorek{ii}(8:32,1),scorek{ii}(8:32,2),scorek{ii}(8:32,3),'Color',colorcodeid(ii,:))
            %     % %         end
            %     %     fnplt(cscvn(scorek(plot_range,:)'),1); hold on
            %     %     end
            %     %     score1 = score(1:size(score,1)/length(uniquetone),:);
            %     %     score2 = score(size(score,1)/length(uniquetone)+1:size(score,1)*2/length(uniquetone),:);
            %     %     score3 = score(size(score,1)*2/length(uniquetone)+1:size(score,1)*3/length(uniquetone),:);
            %
            %     %     figure;
            %     %     plot3(score1(:,1),score1(:,2),score1(:,3),'b'); hold on;
            %     %     for i=1:51%-ts_start+1:floor(num_fr_stim*4/3) - ts_start
            %     %         scatter3(score1(i,1),score1(i,2),score1(i,3), 3+0.5*i,'MarkerColor',colorcodeid)
            %     %         hold on;
            %     %     end
            %     %     scatter3(score1(1,1),score1(1,2),score1(1,3), 20,'b','filled'); hold on
            %     %     xlabel('PC1');ylabel('PC2'); zlabel('PC3');
            %     %     title('Transition of IC Population over Time Stim1 vs Stim 5')
        case 'decoding'
            clear DecodingWeight concaactivity DecodingAccuracytonetc InterestingNeuronsvalue InterestingNeurons
            % decoding pairwise
            loopvalue = 50;
            correctedstim = datarois(fn).correctedstim;
            responsewindow = 10:19;
            for ii = 1:length(correctedstim)
                for iif = length(correctedstim):-1:1
                    for lp = 1:loopvalue
                        clear RadomSelectedTarget roiMatevokedstim concaactivity RadomSelectedFoil testtriallength vec1 vec2 b TestFoil TestTarget VectId classif
                        roiMatevokedstim = datarois(fn).roiMatevokedstim;
                        if ii == iif
                            traintriallength = (size(roiMatevokedstim,3)/2)-1;
                            RadomSelectedTarget = datasample(1:size(roiMatevokedstim,3),traintriallength,'Replace',false);
                            RadomSelectedFoil = datasample(find(~ismember(1:size(roiMatevokedstim,3),RadomSelectedTarget)),traintriallength,'Replace',false);
                            TestTarget = (find(~ismember(1:size(roiMatevokedstim,3),[RadomSelectedTarget,RadomSelectedFoil])));
                            TestFoil = find(~ismember(1:size(roiMatevokedstim,3),[RadomSelectedTarget,RadomSelectedFoil,TestTarget]));
                        else
                            traintriallength = size(roiMatevokedstim,3)-1;
                            RadomSelectedTarget = datasample(1:size(roiMatevokedstim,3),traintriallength,'Replace',false);
                            RadomSelectedFoil = datasample(1:size(roiMatevokedstim,3),traintriallength,'Replace',false);
                            TestTarget = 1:size(roiMatevokedstim,3);
                            TestTarget(ismember(TestTarget,RadomSelectedTarget)) = [];
                            
                            TestFoil = 1:size(roiMatevokedstim,3);
                            TestFoil(ismember(TestFoil,RadomSelectedFoil)) = [];
                        end
                        
                        
                        TestTrials = [TestTarget,TestFoil];
                        VectId = [ones(1,length(TestTarget)),ones(1,length(TestFoil))*2];
                        
                        vec1 = mean(squeeze(mean(roiMatevokedstim(:,ii,RadomSelectedTarget,responsewindow),4)),2);
                        vecsi = mean(squeeze(mean(roiMatevokedstim(:,iif,RadomSelectedFoil,responsewindow),4)),2);
                        w{lp} = vec1-vecsi;
                        b = -(dot(w{lp},vec1) + dot(w{lp},vecsi))/2;
                        
                        %                                 shufflevec1{lp} = squeeze(mean(mean(Shuffletonetc{ii}(RadomSelectedTarget,targetwindow+10,:),2)));
                        %                                 shufflevec2{lp} = squeeze(mean(mean(Shuffletonetc{ii}(RadomSelectedFoil,baselinewindow+10,:),2)));
                        %                                 shufflew{lp} = vec1{lp}-vec2{lp};
                        %                                 shuffleb = -(dot(shufflew{lp},shufflevec1{lp}) + dot(shufflew{lp},shufflevec2{lp}));
                        concaactivity = [squeeze(mean(mean(roiMatevokedstim(:,ii,TestTarget,11:23),4),2)),squeeze(mean(mean(roiMatevokedstim(:,iif,TestFoil,11:23),4),2))];
                        %                                 shuffleconcaactivity = cat(1,mean(Shuffletonetc{ii}(TestTarget,targetwindow+10,:),2),mean(Shuffletonetc{ii}(TestFoil,baselinewindow+10,:),2));
                        
                        for i = 1:length(TestTrials)
                            classif(i) = dot(squeeze(concaactivity(:,i)),w{lp})+b;
                            %                                         shuffleclassif(i) = dot(squeeze(shuffleconcaactivity(i,1,:)),shufflew{lp})+shuffleb;
                        end
                        
                        ClassifiedTarget = find(classif>0);
                        ClassifiedFoil = find(classif<0);
                        %                                 ShuffleClassifiedTarget = find(shuffleclassif>0);
                        %                                 ShuffleClassifiedFoil = find(shuffleclassif<0);
                        %
                        DecodingAccuracytonetc(ii,iif,lp) = (numel(intersect(ClassifiedTarget,find(VectId==1)))+ numel(intersect(ClassifiedFoil,find(VectId==2))))/(length(TestTarget)+length(TestFoil));
                        %                                 ShuffleDecodingAccuracytonetc{ii}(lp,1) = (numel(intersect(ShuffleClassifiedTarget,find(VectId==1)))+ numel(intersect(ShuffleClassifiedFoil,find(VectId==2))))/(length(TestTarget)+length(TestFoil));
                        DecodingWeight(ii,iif,:,lp) = w{lp};
                    end
                    InterestingNeurons{ii,iif} = [find(mean(DecodingWeight(ii,iif,:,:),4)>0.1);find(mean(DecodingWeight(ii,iif,:,:),4)<-0.1)];
                    InterestingNeuronsvalue{ii,iif} = squeeze(squeeze(mean(DecodingWeight(ii,iif,InterestingNeurons{ii,iif}),4)));
                end
            end
            colorredwhite = flipud([ones(100,1),(log10([1.01:0.1:11])/log10(11))',(log10([1.01:0.1:11])/log10(11))']);%([0:0.01:1-0.01].^1.2)',([0:0.01:1-0.01].^1.05)'];
            
            h = figure(); hold on;
            imagesc(1:length(correctedstim),1:length(correctedstim),(squeeze(mean(DecodingAccuracytonetc*100,3)))); caxis([50,75]); colormap(colorredwhite); axis tight;
            ylabel('stim Id'); xlabel('stim Id');
            set(gca,'YDir','reverse')
            xticks([1:length(correctedstim)])
            xticklabels({stimwavname(correctedstim).name})
            xtickangle(45)
            yticks((([1:length(correctedstim)])))
            yticklabels(({stimwavname(correctedstim).name}))
            ytickangle(45)
            colorbar;
            saveas(h,['K:\Jenni\' animal '\figure\social\Population\site' num2str(siteind(analysefilelistind(fn))) '\' path{analysefilelistind(fn)}(end-4:end) '_PairwiseDecodingAccuracy.pdf'])
            %
            %
            %     %     figure(); hold on;
            %     %     for sti = 1:length(correctedstim)
            %     %         subplot(6,3,sti); hold on;
            %     %         plot(squeeze(mean(DecodingWeight(sti,:,:,:),4))');
            %     %         title(batimagingdata.wavnames{sti}(1:end-4))
            %     %     end
            %     %
            %     %     xpsth = ((1:size(roiMatevokedstim,4))-10)/31.25;
            %     %     nbcells = size(roiMatevoked,1);
            %     %     % plot pop imagesc
            %     %     figure(); hold on;
            %     %     for sti = 1:17
            %     %         subplot(6,3,sti); hold on;
            %     %         imagesc(xpsth,1:nbcells,(squeeze(median(roiMatevokedstim(:,sti,:,:),3)))); axis tight; caxis([0.8,1.2]); colormap(redblue(100))
            %     %         xlabel('Time from sound onset [s]'); ylabel('cell nb');
            %     %         title(batimagingdata.wavnames{sti}(1:end-4))
            %     %     end
            %     %
            %     %     subplot(6,3,sti+1); hold on;
            %     %     imagesc(xpsth,1:sti,squeeze(mean(squeeze(median(roiMatevokedstim(:,:,:,:),3)))));axis tight; caxis([0.8,1.2]); colormap(redblue(100))
            %     %     xlabel('Time from sound onset [s]'); ylabel('cell nb');
            %     %     title('Population PSTH')
            %     %
            %     %     figure(); hold on;
            %     %     imagesc(1:sti,1:nbcells,squeeze(mean(squeeze(median(roiMatevokedstim(:,:,:,10:24),3)),3)));axis tight; caxis([0.8,1.2]); colormap(redblue(100))
            %     %     xlabel('Stim nb'); ylabel('cell nb');
            %     %
            %     %     uniqueinterestingneurons = unique(vertcat(InterestingNeurons{:}));
            %     %     for cc = 1:nbcells
            %     %         rawtcstim(cc,:) = max(squeeze(median(roiMatevokedstim(cc,:,:,11:24),3))');
            %     %         tcstim(cc,:) = (rawtcstim(cc,:)-min(rawtcstim(cc,:)))./(max(rawtcstim(cc,:))-min(rawtcstim(cc,:)));
            %     %         bandwidthstim(cc) = sum(tcstim(cc,:));
            %     %     end
            %     % figure(); hold on;
            %     % for cc = 1:length(uniqueinterestingneurons(20:30,:))
            %     %     subplot(5,2,cc); hold on;
            %     %     plot(xpsth,squeeze(median(roiMatevokedstim(uniqueinterestingneurons(cc),:,:,:),3))')
            %     end
    end
    
end

end