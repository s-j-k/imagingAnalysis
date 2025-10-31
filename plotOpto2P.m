% load the data
animalpath='E:\sk274\';
session='10Hz_FOV1\';
datapath=[animalpath session session 'plane0\Fall.mat'];
load(datapath);
pclampPath = [animalpath session];
cd(pclampPath);
% pclampList = dir('*.abf');

files=dir('*.sbx');
for oo=1:length(files)
    file = files(oo).name;
%     pclampFile=pclampList(oo).name;
%     [frame_recording,opto_recording,sound_recording,numframes] = optoReadAbf(pclampFile,pclampPath);
    filename = file(1:length(file)-4);
    x=sbxread(filename,1,100);
    global info_loaded info
    files(oo).frameCount = info.max_idx;
    if oo==1
        files(oo).frameIdx = [1;files(oo).frameCount];
    else
        files(oo).frameIdx = [files(oo-1).frameIdx(2)+1; files(oo-1).frameIdx(2)+1+files(oo).frameCount];
    end
end

Fnorm=F;
for ii = 1:size(F,1)
    Fnorm(ii,:)=(F(ii,:)-median(F(ii,:)))./median(F(ii,:));
end

for tt=1:length(files)
    
    if tt==1
        numRep =(20:100:files(1).frameIdx(2)); % where the sound happens for a session
        files(tt).NumRep=numRep;
    else
        clear numRep
        numRep = (files(tt).frameIdx(1)+20:100:files(tt).frameIdx(2));
        files(tt).NumRep=numRep;
    end
end

for pp=1:length(files)
    for ll = 1:length(files(pp).NumRep)
        for ii=1:size(Fnorm,1)
            frameStart = files(pp).NumRep(ll);
            matEvoked{pp}(ii,ll,:) = Fnorm(ii,frameStart-10:frameStart+40);
            matEvokedTrialBaseline{pp}(ii,ll,:) = (F(ii,frameStart-10:frameStart+40)-mean(F(ii,frameStart-10:frameStart),2))./mean(F(ii,frameStart-10:frameStart),2);
        end
    end
end
% sess=2;
% figure;imagesc(squeeze(matEvoked{sess}(1,:,:))); % for one cell
% ylabel('trial');xlabel('frames');
% 
% figure;imagesc(squeeze(mean(matEvoked{sess}(:,:,:)))); % for all cells, file 1
% ylabel('trial');xlabel('frames');colorbar;caxis([-0.05 0.2]);
% colormap gray
% 
% figure;imagesc(squeeze(mean(matEvoked{sess+1}(:,:,:)))); % for all cells, file 2
% ylabel('trial');xlabel('frames');colorbar;caxis([-0.05 0.2]);
% colormap gray

cm=colormap(winter(10));
% psth
for qq=1:length(matEvokedTrialBaseline)
    f=figure(qq); 
    f.Position = [100 100 440 700];
    hold on; subplot(3,1,[1,2]);imagesc(squeeze(mean(matEvokedTrialBaseline{qq}(:,:,:))));colormap parula;
    caxis([-0.05 0.2]);xticks([1 12 30 45]);ylabel('trial');xlabel('seconds after tone');
        xticklabels({-0.5 0 0.5 1});
    if qq>=4
        ww2=10;
        for ww=1:ww2
            hold on;subplot(3,1,3);
            plot(1:51,squeeze(mean(matEvokedTrialBaseline{qq}(:,ww,:))),'Color',cm(ww,:)); % cell, trial, frame
            ylabel('dF/F');ylim([-0.1 0.2]); 
            xlim([1 50]);xticks([1 16 30 45]);
            xticklabels({-0.5 0 0.5 1});xlabel('seconds after tone');
        end
    else
        for ww=1:size(matEvokedTrialBaseline{qq},2)
            hold on;subplot(3,1,3);
            plot(1:51,squeeze(mean(matEvokedTrialBaseline{qq}(:,ww,:))),'Color',cm(ww,:)); % cell, trial, frame
            ylabel('dF/F');ylim([-0.1 0.2]); 
            xlim([1 50]);xticks([1 16 30 45]);
            xticklabels({-0.5 0 0.5 1});xlabel('seconds after tone');
        end
    end
    saveas(gcf,[char(session) 'population avg session ' num2str(qq) ' trial ' num2str(ww)]);
    saveas(gcf,[char(session) 'population avg session ' num2str(qq) ' trial ' num2str(ww) '.png']);
    
end
save([char(session) '_workspace.mat']);
clear all
close all
%%
% to quantify effect of opto, look at avg 10 frames before sound for each
% trial and compare baseline between light off vs. light on
% we think the baseline is going up
% plot is baseline avg on y axis, trial number on x axis
% scatter plot with each condition 
animalpath='E:\sk274\';
sessList=["20Hz_FOV1\","10Hz_FOV1\","20Hz_FOV3\","20Hz_FOV3\"];
bigMat={};
bigMat{1,1}='matEvoked Trial Baseline';bigMat{1,2}='matEvoked population Baseline';
bigMat{1,3}='file metadata';
for ee=1:length(sessList)
    datapath=[animalpath char(sessList(ee)) char(sessList(ee))];
    cd(datapath)
    load('_workspace.mat', 'matEvokedTrialBaseline','matEvoked','files')
    bigMat{ee+1,1} = matEvokedTrialBaseline;
    bigMat{ee+1,2} = matEvoked;
    bigMat{ee+1,3} = files;
end
save('allOptoData.mat','bigMat');
%%
clear baseAvg
% take each trial and average the baseline, using the session population
% baseline (bigMat column 2) or trial baseline (bigMat column 1)
for dd=2:5
    sess=bigMat{dd,2};
    for rr=1:size(sess,2)
        allTrials = sess(1,rr); % cell, trial, frame
        if size(allTrials,2)>20
            % this is a tuning curve file. ignore
        else
            baseAvg{dd,rr}=squeeze(mean(allTrials{:}(:,:,1:10),3));
            peakAvg{dd,rr}=squeeze(mean(allTrials{:}(:,:,16:26),3));
        end
    end
end

% 1 = 10Hz
% 2 = 20Hz
% 3 = no light
% 4 = baseline
conditions = [0 2 2 2 2 2 2 0; ...
              1 1 1 1 3 1 1 4; ...
              2 2 2 2 2 0 3 0; ...
              0 3 1 1 1 1 1 0];
          
save('allOptoData.mat','bigMat','conditions');
%% scatter plot for each condition, mean baseline and mean peak
% baseAvg{1,:}=[];
tenIdx=baseAvg(conditions==1);
for aa=1:length(tenIdx)
    tenh(aa) = cellfun(@(x) mean(mean(x(:,1:size(tenIdx{aa},2)))),tenIdx(aa));
end
twenIdx=baseAvg(conditions==2);
for aa=1:length(twenIdx)
    twenh(aa) = cellfun(@(x) mean(mean(x(:,1:size(twenIdx{aa},2)))),twenIdx(aa));
end
noLightIdx=baseAvg(conditions==3);
for aa=1:length(noLightIdx)
    nolight(aa) = cellfun(@(x) mean(mean(x(:,1:size(noLightIdx{aa},2)))),noLightIdx(aa));
end
figure; subplot(1,2,1); hold on; 
scatter(1,mean(nolight),55,[1 0 0],'filled'); 
scatter(2,mean(tenh),55,[0 0 1],'filled');
scatter(3,mean(twenh),55,[1 0 1],'filled');
plot(repmat(1,1,length(nolight)),nolight,'.','Color',[0.75,0,0],'MarkerSize',15); 
plot(repmat(2,1,length(tenh)),tenh,'.','Color',[0 0 0.75],'MarkerSize',15);
plot(repmat(3,1,length(twenh),1),twenh,'.','Color',[0.75 0 0.75],'MarkerSize',15);
% legend('no light','10Hz','20Hz'); legend box off
xlim([0.5 3.5]);ylim([-0.1 0.1]);
ylabel('Mean baseline dF/F');
xticks([1 2 3]);xticklabels({'no light','10Hz','20Hz'});xlabel('condition');
title('Mean baseline by condition');

% peak
tenIdx=peakAvg(conditions==1);
for aa=1:length(tenIdx)
    tenh(aa) = cellfun(@(x) mean(mean(x(:,1:size(tenIdx{aa},2)))),tenIdx(aa));
end
twenIdx=peakAvg(conditions==2);
for aa=1:length(twenIdx)
    twenh(aa) = cellfun(@(x) mean(mean(x(:,1:size(twenIdx{aa},2)))),twenIdx(aa));
end
noLightIdx=peakAvg(conditions==3);
for aa=1:length(noLightIdx)
    nolight(aa) = cellfun(@(x) mean(mean(x(:,1:size(noLightIdx{aa},2)))),noLightIdx(aa));
end

subplot(1,2,2); hold on; 
scatter([1],mean(nolight),55,[1 0 0],'filled'); 
scatter(2,nanmean(tenh),55,[0 0 1],'filled');
scatter(3,nanmean(twenh),55,[1 0 1],'filled');
plot(repmat(1,1,length(nolight)),nolight,'.','Color',[0.75,0,0],'MarkerSize',15); 
plot(repmat(2,1,length(tenh)),tenh,'.','Color',[0 0 0.75],'MarkerSize',15);
plot(repmat(3,1,length(twenh),1),twenh,'.','Color',[0.75 0 0.75],'MarkerSize',15);
% legend('no light','10Hz','20Hz'); legend box off
xlim([0.5 3.5]);ylim([-0.1 0.1]);
ylabel('Mean peak dF/F');
xticks([1 2 3]);xticklabels({'no light','10Hz','20Hz'});xlabel('condition');
title('Mean peak by opto condition');


%% 
%baseline
tenIdx=baseAvg(conditions==1);
twenIdx=baseAvg(conditions==2);
noLightIdx=baseAvg(conditions==3);
howmanyTrials = 5;
if howmanyTrials==0
    for aa=1:length(tenIdx)
        tenh(aa) = cellfun(@(x) mean(mean(x(:,1:size(tenIdx{aa},2)))),tenIdx(aa));
    end
    
    for aa=1:length(twenIdx)
        twenh(aa) = cellfun(@(x) mean(mean(x(:,1:size(twenIdx{aa},2)))),twenIdx(aa));
    end
    
    for aa=1:length(noLightIdx)
        nolight(aa) = cellfun(@(x) mean(mean(x(:,1:size(noLightIdx{aa},2)))),noLightIdx(aa));
    end
else
    for aa=1:length(tenIdx)
        tenh(aa) = cellfun(@(x) mean(mean(x(:,1:howmanyTrials))),tenIdx(aa));
    end
    twenIdx=baseAvg(conditions==2);
    for aa=1:length(twenIdx)
        twenh(aa) = cellfun(@(x) mean(mean(x(:,1:howmanyTrials))),twenIdx(aa));
    end
    noLightIdx=baseAvg(conditions==3);
    for aa=1:length(noLightIdx)
        nolight(aa) = cellfun(@(x) mean(mean(x(:,1:howmanyTrials))),noLightIdx(aa));
    end
end
conditionsMatrixB=padcat(nolight,tenh,twenh);
[~,~,aov]=anova1(conditionsMatrixB');
aov
title('anova comparing baseline, first 5 trials');
multcompare(aov)
%%
% peak
tenIdx=peakAvg(conditions==1);
twenIdx=peakAvg(conditions==2);
noLightIdx=peakAvg(conditions==3);
if howmanyTrials==0
    for aa=1:length(tenIdx)
        tenh(aa) = cellfun(@(x) mean(mean(x(:,1:size(tenIdx{aa},2)))),tenIdx(aa));
    end
    
    for aa=1:length(twenIdx)
        twenh(aa) = cellfun(@(x) mean(mean(x(:,1:size(twenIdx{aa},2)))),twenIdx(aa));
    end
    
    for aa=1:length(noLightIdx)
        nolight(aa) = cellfun(@(x) mean(mean(x(:,1:size(noLightIdx{aa},2)))),noLightIdx(aa));
    end
else
    for aa=1:length(tenIdx)
        tenh(aa) = cellfun(@(x) mean(mean(x(:,1:howmanyTrials))),tenIdx(aa));
    end
    for aa=1:length(twenIdx)
        twenh(aa) = cellfun(@(x) mean(mean(x(:,1:howmanyTrials))),twenIdx(aa));
    end
    for aa=1:length(noLightIdx)
        nolight(aa) = cellfun(@(x) mean(mean(x(:,1:howmanyTrials))),noLightIdx(aa));
    end    
end
conditionsMatrixP=padcat(nolight,tenh,twenh);
[~,~,aov]=anova1(conditionsMatrixP');title('anova comparing peak, first 5 trials');
aov
% multcompare(aov)

% baseline subtracted from peak
conditionsDiff=conditionsMatrixP-conditionsMatrixB;
[~,~,aov]=anova1(conditionsDiff');title('anova comparing peak-baseline');
multcompare(aov)

%%

