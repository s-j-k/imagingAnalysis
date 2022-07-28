function Bat_movie_Analysis()

animal = 'RoofBuddy1';
datatable = readtable(['N:\Jenni\' animal '\ProcessingProgress.xlsx'],'Format','auto');

% manualroisdone = find(table2array(datatable(:,13))==1);
daytraining = table2array(datatable(:,1));
% preprocessedrois = unique(daytraining(manualroisdone));
framerate = table2array(datatable(:,8));
protocol = table2array(datatable(:,7));
path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,10));
uniquedaytraining = unique(daytraining);
site = table2array(datatable(:,5));
uniquesite = unique(site);
for ii = 1:length(uniquesite)
    analysefilelistind(ii) = find(site==uniquesite(ii),1);
end

% uniquedaytraining = unique(daytraining);
% for ii = 1:length(uniquedaytraining)
%     analysefilelistind(ii) = find(daytraining==uniquedaytraining(ii),1);
% end


%   for tid = 1:length(toneid)
%         clear bla bli
%     bla = max(max(max(MeanEvokedMovieTonotopy{tid},[],3)));
%     bli = MeanEvokedMovieTonotopy{tid}./bla;
%     figure(1); hold on; subplot(1,length(toneid),tid);imagesc(mean(bli(:,:,11:15),3));
%   end
    
for fn = 11%:length(analysefilelistind)
%for fn = 1:length(analysefilelistind(1:end))
    clearvars -except fn animal datatable manualroisdone daytraining preprocessedrois framerate protocol path behaviorfile uniquedaytraining analysefilelistind
    % read frames from h5
    h5list = dir(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' '*.h5']);
    nbframestmp(1) = 0;
    for filenb = 1:length(h5list)
        clear hinfo
        hinfo = hdf5info(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\' h5list(filenb).name]);
        nbframestmp(filenb+1) = hinfo.GroupHierarchy.Datasets.Dims(3);
    end
    cumsumframes = cumsum(nbframestmp);
    
    cd(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\'])
    load('EvokedMovie.mat')
%     stimnumber = [];
    behavtmp = behaviorfile(intersect(find(daytraining==uniquedaytraining(fn)),find(framerate==31.25)),:);
    for bf = 1:size(behavtmp,1)
        clear batimagingdata toneid tonedur MeanEvokedMovieTonotopy MeanEvokedMovieTonotopy MeanEvokedMovie MeanEvokedMovieDuration
        load(['N:\Jenni\' animal '\behavior\' behavtmp{bf} '.mat'])
%         stimnumber = [stimnumber, 1:length(batimagingdata.frametone)];%1:400:44000;
        if strfind(behavtmp{bf},'tonotopy')>1
            stimid{bf} = batimagingdata.toneseq;
            stimduration{bf} = batimagingdata.durationseq;
            toneid = unique(batimagingdata.toneseq);
            tonedur = unique(batimagingdata.durationseq);
            for tid = 1:length(toneid)
                MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,find(stimid{bf}==toneid(tid)),:),3));
                for tdu = 1:length(tonedur)
                    MeanEvokedMovieDuration{tid}(:,:,:,tdu) = mean(EvokedMovie(:,:,intersect(find(stimid{bf}==toneid(tid)),find(stimduration{bf}==tonedur(tdu))),:),3);
                end
            end
        end
        if strfind(behavtmp{bf},'Linear')>1
            stimid{bf} = batimagingdata.toneseq;
            stimduration{bf} = batimagingdata.durationseq;
            toneid = unique(batimagingdata.toneseq);
            tonedur = unique(batimagingdata.durationseq);
            for tid = 1:length(toneid)
                MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,find(stimid{bf}==toneid(tid)),:),3));
                for tdu = 1:length(tonedur)
                    MeanEvokedMovieDuration{tid}(:,:,:,tdu) = mean(EvokedMovie(:,:,intersect(find(stimid{bf}==toneid(tid)),find(stimduration{bf}==tonedur(tdu))),:),3);
                end
            end
        end
        if strfind(behavtmp{bf},'Angie')>1
            stimid{bf} = batimagingdata.toneseq;
            stimduration{bf} = batimagingdata.durationseq;
            toneid = unique(batimagingdata.toneseq);
            for tid = 1:length(toneid)
                MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,find(stimid{bf}==toneid(tid)),:),3));
            end
        end
        
        if strfind(behavtmp{bf},'Pulse')>1
            stimid{bf} = batimagingdata.toneseq;
            stimduration{bf} = batimagingdata.durationseq;
            toneid = unique(batimagingdata.toneseq);
            for tid = 1:length(toneid)
                MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,find(stimid{bf}==toneid(tid)),:),3));
            end
        end
        
        if strfind(behavtmp{bf},'Mel')>1
            stimid{bf} = batimagingdata.toneseq;
            stimduration{bf} = batimagingdata.durationseq;
            toneid = unique(batimagingdata.toneseq);
            for tid = 1:length(toneid)
                MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,find(stimid{bf}==toneid(tid)),:),3));
            end
        end
    end
    
    
end
end