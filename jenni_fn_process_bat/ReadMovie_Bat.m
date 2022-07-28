function ReadMovie_Bat()

animal = 'RoofBuddy2';
datatable = readtable(['N:\Jenni\' animal '\ProcessingProgress.xlsx']);

manualroisdone = find(table2array(datatable(:,13))==1);%find(~cellfun(@isempty,cellfun(@(x) find(x=='1'),table2array(datatable(:,13)),'UniformOutput',0)));%find(table2array(datatable(:,13))==1);
daytraining = table2array(datatable(:,1));
site = table2array(datatable(:,5));
siteind = table2array(datatable(:,5));
preprocessedrois = unique(daytraining(manualroisdone));
framerate = (table2array(datatable(:,8)));
protocol = table2array(datatable(:,7));
path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,10));
uniquedaytraining = unique(daytraining);
uniquesite = unique(site);
for ii = 1:length(uniquesite)
    analysefilelistind(ii) = find(site==uniquesite(ii),1);
end

uniquesite = unique(siteind);

for sitenb = 3%1:length(uniquesite)%1:length(analysefilelistind(1:end))
    frameratesite = unique(framerate(find(siteind==uniquesite(sitenb))));
    frameratesite(isnan(frameratesite)) = [];
    for fr = 2%1:length(frameratesite)
            clearvars -except siteind sitenb fr frameratesite fn framerate uniquesite site uniquesite animal datatable manualroisdone daytraining preprocessedrois protocol path behaviorfile uniquedaytraining analysefilelistind
        % read frames from h5
        h5list = dir(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\' '*.h5']);
        nbframestmp(1) = 0;
        for filenb = 1:length(h5list)
            clear hinfo
            hinfo = hdf5info(['N:\Jenni\' animal '\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\' h5list(filenb).name]);
            nbframestmp(filenb+1) = hinfo.GroupHierarchy.Datasets.Dims(3);
        end
        cumsumframes = cumsum(nbframestmp);
        
        % first get behavior file
        startstim = [];
        fileidstim = [];
        behavtmp = behaviorfile(intersect(find(siteind==uniquesite(sitenb)),find(framerate==frameratesite(fr))));%%find(cellfun(@(x) strcmp(x,'31.25'),framerate)==1)),:);
        for bf = 1:size(behavtmp,1)
            clear batimagingdata
            load(['N:\Jenni\' animal '\behavior\' behavtmp{bf} '.mat'])
            startstim = [startstim, batimagingdata.frametone+cumsumframes(bf)];%1:400:44000;
            fileidstim = [fileidstim; ones(size(batimagingdata.frametone,2),1)*bf];
        end
        startstim = startstim';
        for ii = 1:length(startstim)
            matframes(ii,:) = startstim(ii)-10:startstim(ii)+44;
        end
        
        
        %  Get raw aligned movie
        cd(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\suite2p\plane0\'])
        
        nimg = sum(nbframestmp);%23173;
        
        if ~exist(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\TonotopyMeanEvokedMovie.mat'])
            
            fileID = fopen('data.bin','r'); % open green channel binary file
            
            %    cd(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\31.25\suite2p\plane0'])
            data = load(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\31.25\suite2p\plane0\Fall.mat']);
            
            if frameratesite(fr) == 31.25
                ly = data.ops.Ly;
                lx = data.ops.Lx;
            elseif frameratesite(fr) == 62.5
                ly = data.ops.Ly/2;
                lx = data.ops.Lx;
            elseif frameratesite(fr) == 125
                ly = data.ops.Ly/4;
                lx = data.ops.Lx;
            end
            k=0; a=1;
            blksize = 2000;%2000; % nb of frames loaded at a time (depend on RAM)
            to_read = min(blksize,nimg-k);
            avgA = [];
            avgA = nan(lx,ly);
            while to_read>0
                clear A
                A = fread(fileID,ly*lx*to_read,'*int16');
                A = reshape(A,lx,ly,[]);
                datamovie{a} = A;
                avgA(:,:,a) = mean(A,3);
                a=a+1;
                k = k+to_read;
                to_read = min(blksize,nimg-k);
            end
            
            moviedata = [];
            for i = 1:length(datamovie)
                moviedata = cat(3,moviedata,datamovie{i});
            end
            meanimg = mean(moviedata,3);
            meanpoppsth = squeeze(mean(squeeze(mean(moviedata,1)),1));
            
            
            MatAllEvoked = moviedata(:,:,matframes);
            EvokedMovie = reshape(MatAllEvoked,[size(moviedata,1),size(moviedata,2),size(startstim,1),55]);
            
            
            %     save(['N:\Jenni\' animal '\' path{analysefilelistind(fn)} '\31.25\EvokedMovie.mat'],'EvokedMovie','-v7.3')
            for bf = 1:size(behavtmp,1)
                clear numstimind
                numstimind = find(fileidstim==bf);
                clear batimagingdata toneid tonedur MeanEvokedMovieTonotopy MeanEvokedMovieTonotopy MeanEvokedMovie MeanEvokedMovieDuration
                load(['N:\Jenni\' animal '\behavior\' behavtmp{bf} '.mat'])
                %         stimnumber = [stimnumber, 1:length(batimagingdata.frametone)];%1:400:44000;
                if strfind(behavtmp{bf},'tonotopy')>1
                    stimid{bf} = batimagingdata.toneseq;
                    stimduration{bf} = batimagingdata.durationseq;
                    toneid = unique(batimagingdata.toneseq);
                    tonedur = unique(batimagingdata.durationseq);
                    for tid = 1:length(toneid)
                        MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,numstimind(numstimind(find(stimid{bf}==toneid(tid)))),:),3));
                        for tdu = 1:length(tonedur)
                            MeanEvokedMovieDuration{tid}(:,:,:,tdu) = squeeze(mean(EvokedMovie(:,:,intersect(numstimind(find(stimid{bf}==toneid(tid))),numstimind(find(stimduration{bf}==tonedur(tdu)))),:),3));
                        end
                    end
                    save(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\TonotopyMeanEvokedMovie.mat'],'MeanEvokedMovie','-v7.3')
                    save(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\TonotopyDurationMeanEvokedMovie.mat'],'MeanEvokedMovieDuration','-v7.3')
                end
                
                if strfind(behavtmp{bf},'Linear')>1
                    stimid{bf} = batimagingdata.toneseq;
                    stimduration{bf} = batimagingdata.durationseq;
                    toneid = unique(batimagingdata.toneseq);
                    tonedur = unique(batimagingdata.durationseq);
                    for tid = 1:length(toneid)
                        MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,numstimind(find(stimid{bf}==toneid(tid))),:),3));
                        for tdu = 1:length(tonedur)
                            MeanEvokedMovieDuration{tid}(:,:,:,tdu) = squeeze(mean(EvokedMovie(:,:,intersect(numstimind(find(stimid{bf}==toneid(tid))),numstimind(find(stimduration{bf}==tonedur(tdu)))),:),3));
                        end
                    end
                    save(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\LinearMeanEvokedMovie.mat'],'MeanEvokedMovie','-v7.3')
                    save(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\LinearDurationMeanEvokedMovie.mat'],'MeanEvokedMovieDuration','-v7.3')
                end
                if strfind(behavtmp{bf},'Angie')>1
                    stimid{bf} = batimagingdata.toneseq;
                    stimduration{bf} = batimagingdata.durationseq;
                    toneid = unique(batimagingdata.toneseq);
                    for tid = 1:length(toneid)
                        MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,numstimind(find(stimid{bf}==toneid(tid))),:),3));
                    end
                    save(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\AngieMeanEvokedMovie.mat'],'MeanEvokedMovie','-v7.3')
                end
                
                if strfind(behavtmp{bf},'Pulse')>1
                    stimid{bf} = batimagingdata.toneseq;
                    stimduration{bf} = batimagingdata.durationseq;
                    toneid = unique(batimagingdata.toneseq);
                    for tid = 1:length(toneid)
                        MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,numstimind(find(stimid{bf}==toneid(tid))),:),3));
                    end
                    save(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\PulseMeanEvokedMovie.mat'],'MeanEvokedMovie','-v7.3')
                end
                
                if strfind(behavtmp{bf},'Mel')>1
                    stimid{bf} = batimagingdata.toneseq;
                    stimduration{bf} = batimagingdata.durationseq;
                    toneid = unique(batimagingdata.toneseq);
                    for tid = 1:length(toneid)
                        MeanEvokedMovie{tid} = squeeze(mean(EvokedMovie(:,:,numstimind(find(stimid{bf}==toneid(tid))),:),3));
                    end
                    save(['N:\Jenni\RoofBuddy2\' path{find(siteind==uniquesite(sitenb),1)} '\' num2str(frameratesite(fr)) '\MelMeanEvokedMovie.mat'],'MeanEvokedMovie','-v7.3')
                end
            end
        end
    end
end
% seqid = batimagingdata.durationseq;%batimagingdata.;%repmat(1:11,1,10);
% toneid = batimagingdata.toneseq;%[13454;4757;26907;53813;74264;76424;19026;38052;6727;9513;70000];
% uniquetone = unique(toneid);
% % [orderfreq,orderstim] = sort(toneid);
%
% % check evoked response over the course of the imaging session
% plot(mean(squeeze(mean(moviedata))));
% hold on; plot([startstim;startstim],[zeros(1,length(startstim));ones(1,length(startstim))*4000],'k');
%
% for ii = 1:length(startstim)
% Matevoked{ii} = moviedata(:,:,startstim(ii)-10:startstim(ii)+40);%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
% Matevokednorm{ii} = moviedata(:,:,startstim(ii)-10:startstim(ii)+40)./median(moviedata(:,:,startstim(ii)-10:startstim(ii)),3);
% end
% for ii = 1:length(uniquetone)
%     Matevokedstim{ii} = mean(cat(4,Matevoked{find(toneid == uniquetone(ii))}),4); %mean(cat(4,Matevoked{find(seqid == ii)}),4);
%     Matevokedstimnorm{ii} = mean(cat(4,Matevokednorm{find(toneid == uniquetone(ii))}),4);
% end
% for ii = 1:length(uniquetone)
%     Matevokedstimtrial{ii} = cat(4,Matevoked{find(toneid == uniquetone(ii))}); %mean(cat(4,Matevoked{find(seqid == ii)}),4);
%     Matevokedstimnormtrial{ii} = cat(4,Matevokednorm{find(toneid == uniquetone(ii))});
% end
%
% load('N:\Jenni\RoofBuddy1\012_h5\site2\31.25\31.25_TC_plane0.mat')
% TC = TC./median(TC')';
% for cc = 1:size(TC,1)
% for ii = 1:length(startstim)
% roiMatevoked(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40);%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
% % roiMatevokednorm(cc,ii,:) = TC(cc,startstim(ii)-10:startstim(ii)+40)./median(TC(cc,startstim(ii)-10:startstim(ii)),2);
% end
% for ii = 1:length(uniquetone)
%     roiMatevokedstim(cc,ii,:,:) = squeeze(roiMatevoked(cc,find(toneid == uniquetone(ii)),:));%./median(moviedata(:,:,startstim(ii):startstim(ii)+200),3);
% end
%
% % plot mean evoked response per roi
% figure(2); hold on; subplot(3,1,[1,2]); imagesc(squeeze((median(roiMatevoked,2)))); colormap gray;
% subplot(3,1,[3]); plot(median(squeeze((mean(roiMatevoked,2)))),'k'); xlim([1,50])
%
% % plot tuning 247
%
%
% for i = 1:11
% stimpopmat{i} = cell2mat(cellfun(@(x) mean(squeeze(mean(x,1)),1),Matevoked(find(seqid == i)),'UniformOutput',0)');%./median(median(cell2mat(cellfun(@(x) mean(squeeze(mean(x(:,:,1:20),1)),1),Matevoked(find(seqid == ii)),'UniformOutput',0)')));
% end
% % plot population resp
% cmap = [jet(8); [0 0 0];  [0.4 0.4 0.4]; [0.7 0.7 0.7]];
% figure(); hold on;
% for i = 1:11
%     bla(i) = plot([1:201]+(201*i-1),median(stimpopmat{orderstim(i)}),'Color',cmap(i,:));
% end
% % legend(bla(i).Line,'DS')
%
% % open IMJ
% disk = 'E:';
% javaaddpath ([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\mij.jar'])
% javaaddpath ([disk '\KishoreLab\Shared\Matlab\preprocessing\MIJI\ij-1.52i.jar'])
% MIJ.start([disk '\KishoreLab\Shared\Matlab\Fiji.app'])
% open movie for stim average
% MIJ.createImage(uint16(moviedata(:,:,1:200)))
% MIJ.createImage((uint16(Matevokedstimnorm{1}(:,:,:))))%mean(uint16(Matevokedstim{1}(:,:,:)),3))

end