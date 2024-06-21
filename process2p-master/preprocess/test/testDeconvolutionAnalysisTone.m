clear;
datapath ='C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\smooth3000\behavior';
behavpath = 'C:\Users\zzhu34\Documents\tempdata\excitData\cd036\behavior';
%filenameList = {'cd036_confoo_smin.mat','cd036_foo_smin.mat','cd036_thre_smin.mat'};
%filenameList = {'cd036_calman_confoo_optb.mat','cd036_calman_foo90_optb_nosmin.mat',...
%    'cd036_calman_foo95_optb_nosmin.mat','cd036_calman_foo95_optb.mat','cd036_calman_foo99_optb.mat','cd036_calman_thre95_optb.mat'};
filenameList = {'cd036_calman_ar2_confoo_optb.mat','cd036_calman_ar2_foo90_optb_nosmin.mat','cd036_calman_ar2_foo95_optb_nosmin.mat',...
    'cd036_calman_ar2_foo90_optb.mat','cd036_calman_ar2_foo95_optb.mat','cd036_calman_ar2_thre95_optb.mat'};
%filenameList = {'cd036_foo_lambda0.mat','cd036_foo_lambda05.mat','cd036_foo_lambda10.mat'};
shortFilename = [];

peakDetection = false; peakFrames = 1:6;
for i = 1:length(filenameList); shortFilename{i} = filenameList{i}(14:end-4); shortFilename{i}(strfind(shortFilename{i},'_'))=' ';end
for i = 1:length(filenameList)
    disp([filenameList{i} ' started!'])
    load([datapath '\' filenameList{i}],'s','c','day');
    
    for j = 1:length(day)
        s{j} = [nan(size(s{j},1),1) s{j}];
        behavMatrix = importdata([behavpath '\cd036_' int2str(day(j)-1) 'v1.txt']);
        [tempTPM, tempFPM,tempT,tempF,tempTPT,tempFPT] = getPeakAct(behavMatrix,s{j},peakDetection,peakFrames);
        tPeakMean(:,j,i) = tempTPM;fPeakMean(:,j,i) = tempFPM;
        tPeakTrial(:,:,j,i) = tempTPT;fPeakTrial(:,:,j,i) = tempFPT;
        tPeak(:,:,:,j,i) = tempT;fPeak(:,:,:,j,i) = tempF;
    end
end
load([datapath '\selectDff.mat'],'selectDff');
disp('selected dff started!')
peakDetection = true; peakFrames = 1:12;
for j = 1:length(day)
    selectDff{j} = [nan(size(selectDff{j},1),1) selectDff{j}];
    behavMatrix = importdata([behavpath '\cd036_' int2str(day(j)-1) 'v1.txt']);
    [tempTP, tempFP, tempT, tempF,tempTPT,tempFPT] = getPeakAct(behavMatrix,selectDff{j},peakDetection,peakFrames);
    tPeakMeanDff(:,j) = tempTP;fPeakMeanDff(:,j) = tempFP;
    tPeakTrialDff(:,:,j) = tempTPT;fPeakTrialDff(:,:,j) = tempFPT;
    tPeakDff(:,:,:,j) = tempT;fPeakDff(:,:,:,j) = tempF;
    
end

%% load the data from the whole traces

clear;
mouse = 'cd017';
configpath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\allSessions\'];
behavpath = ['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\behavior'];
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\roi\ishere_plane0.mat'],'ishere');
if iscell(ishere); ishere = ishere{1}; end
peakDetection = false; peakFrames = 1:6;

[nFrames, ~] = func_readnFrames(mouse,'root',configpath);
nFrames = [0 0; nFrames];
configTable =  readtable([configpath '\' mouse '_config.csv']);

filenameList = {[mouse '_calman_ar2_foo90_pars_allday_s.mat'],[mouse '_calman_ar2_foo90_pars_day_s.mat'],...
    [mouse '_calman_ar2_foo95_optb_nosmin_s.mat']};

for i = 1:length(filenameList); shortFilename{i} = filenameList{i}(24:end-6); shortFilename{i}(strfind(shortFilename{i},'_'))=' ';end
day = 5;
selectNeuron = sum(ishere(:,day+1),2)==length(day);
tPeakMean = []; fPeakMean = [];
tPeakTrial = []; fPeakTrial = [];
tPeak =[]; fPeak = [];
for i = 1:length(filenameList)
    disp([filenameList{i} ' started!'])
    load([datapath '\' filenameList{i}],'s');
    sessionIdx = find(configTable.Day == (day) & ~strcmp(configTable.BehavType,'Baseline')); 
    tPeakMeanTemp = []; fPeakMeanTemp = [];
    tPeakTrialTemp = []; fPeakTrialTemp = [];
    tPeakTemp = []; fPeakTemp= [];
    for j = 1:length(sessionIdx)
        tempSession = sessionIdx(j);
        tempS = [nan(size(s{tempSession},1),1) s{tempSession}]; tempS = tempS(selectNeuron,:);      
        behavMatrix = importdata([behavpath '\' configTable.BehavFile{tempSession}]);
        [tempTPM, tempFPM,tempT,tempF,tempTPT,tempFPT] = getPeakAct(behavMatrix,tempS,peakDetection,peakFrames);
        tPeakMeanTemp = cat(2,tPeakMeanTemp,tempTPM);fPeakMeanTemp = cat(2,fPeakMeanTemp,tempFPM) ;
        tPeakTrialTemp = cat(2,tPeakTrialTemp,tempTPT);fPeakTrialTemp = cat(2,fPeakTrialTemp,tempFPT);
        tPeakTemp = cat(3,tPeakTemp,tempT);fPeakTemp = cat(3,fPeakTemp,tempF);
    end
    tPeakMean = cat(2,tPeakMean,mean(tPeakMeanTemp,2));fPeakMean = cat(2,fPeakMean,mean(fPeakMeanTemp,2)) ;
    tPeakTrial = cat(3,tPeakTrial,tPeakTrialTemp);fPeakTrial = cat(3,fPeakTrial,fPeakTrialTemp);
    tPeak = cat(4,tPeak,tPeakTemp);fPeak = cat(4,fPeak,fPeakTemp);
end


%% plot 1 - the averaged activity over days in tone-evoked period
figure; hold on;
dffFlat = [tPeakMeanDff(:)' fPeakMeanDff(:)'];
for i = 1:length(filenameList); spikeFlat(:,i) = [reshape(tPeakMean(:,:,i),1,[]) reshape(fPeakMean(:,:,i),1,[])];end
[N,edges,bin] = histcounts(dffFlat,linspace(0,prctile(dffFlat,98),21));
for j = 1:20 
    tempX(j) = mean(edges(j:j+1));
    for i = 1:length(filenameList)
        tempMean(j,i) = mean(spikeFlat(bin==j,i)); tempSEM(j,i) = std(spikeFlat(bin==j,i))/sqrt(N(j));     
    end
end
for i = 1:length(filenameList); plot(tempX,tempMean(:,i),'Color',matlabColors(i)); end
tempLegend = [];
for i = 1:length(filenameList)
    f = fillErrorbarPlot(tempX,tempMean(:,i)', tempSEM(:,i)',matlabColors(i),'LineStyle','none');
    f.FaceAlpha = 0.1;
    scatter(dffFlat,spikeFlat(:,i),5,matlabColors(i), 'filled' , 'o', 'MarkerFaceAlpha', 0.1); 
end
ylabel('spike'); xlabel('dff'); 
xlim([0 prctile(dffFlat,98)]);
ylim([0 prctile(spikeFlat(:),99)]);
legend(shortFilename{:},'Location','Best');
title('Tone-evoked activity')
%% plot 1.1 - the averaged activity over days in tone-evoked period
figure; 
for j = 1:length(filenameList); subplot_tight(2,3,j,[0.15,0.06]);hold on; title(shortFilename{j})
for i= 1:7
    scatter([tPeakMeanDff(:,i)' fPeakMeanDff(:,i)'],...
        [reshape(tPeakMean(:,i,j),1,[]) reshape(fPeakMean(:,i,j),1,[])],...
        5, 'filled' , 'o','MarkerFaceAlpha', 0.4);
end
ylabel('spike'); xlabel('dff'); 
xlim([0 prctile([tPeakMeanDff(:)' fPeakMeanDff(:)'],99)]);
ylim([0 prctile([reshape(tPeakMean(:,:,j),1,[]) reshape(fPeakMean(:,:,j),1,[])],99)]);
end

%% plot 2.1 - PSTH of different groups of neurons by dff
binEdges = 0:0.01:0.15; binSelection = [2 4 6 8 10];
tPeakMeanDffAvg = mean(tPeakMeanDff,2); [N, edges, tbin] = histcounts(tPeakMeanDffAvg,binEdges);
for i = 1:length(filenameList)
    tPSTH(:,:,i) = mean(reshape(tPeak(:,:,:,:,i),size(tPeak,1), size(tPeak,2),[]),3);
    fPSTH(:,:,i) = mean(reshape(fPeak(:,:,:,:,i),size(fPeak,1), size(fPeak,2),[]),3);
end
tDffPSTH = mean(reshape(tPeakDff,size(tPeakDff,1), size(tPeakDff,2),[]),3);
fDffPSTH = mean(reshape(fPeakDff,size(fPeakDff,1), size(fPeakDff,2),[]),3);
figure; for i = 1:length(binSelection)
    plot(smoothdata(mean(tPSTH(:,:,i),1) - mean(mean(tPSTH(:,1:5,i),1),2),'gaussian',3)); hold on; end
figure;
for i = 1:length(filenameList)
    subplot_tight(3,3,i);
    for j = 1:length(binSelection)
        tempAct = mean(tPSTH(tbin == binSelection(j) ,:,i),1) - mean(mean(tPSTH(tbin == binSelection(j),1:3,i),1),2);
        plot(tempAct); hold on; 
    end
    if i == 1; legend('dff 0.01-0.02','dff 0.03-0.04','dff 0.05-0.06','dff 0.07-0.08','dff 0.09-0.1');end
    title(shortFilename{i})
end
subplot_tight(3,3,length(filenameList)+1); title('dff')
for j = 1:length(binSelection)
    tempAct = mean(tDffPSTH(tbin == binSelection(j) ,:),1) - mean(mean(tDffPSTH(tbin == binSelection(j),1:3),1),2);
    plot(tempAct); hold on; 
end

%% plot 2.2 - PSTH of different groups of neurons by percentile
binEdges = 0:0.01:0.15; binSelection = floor([.05 .2 .4 .6 .8 .95] * length(tPeakMeanDffAvg));
tPeakMeanDffAvg = mean(tPeakMeanDff,2); [~, tSortIdx] = sort(tPeakMeanDffAvg,'ascend');
for i = 1:length(filenameList)
    tPSTH(:,:,i) = mean(reshape(tPeak(:,:,:,:,i),size(tPeak,1), size(tPeak,2),[]),3);
    fPSTH(:,:,i) = mean(reshape(fPeak(:,:,:,:,i),size(fPeak,1), size(fPeak,2),[]),3);
end
tDffPSTH = mean(reshape(tPeakDff,size(tPeakDff,1), size(tPeakDff,2),[]),3);
fDffPSTH = mean(reshape(fPeakDff,size(fPeakDff,1), size(fPeakDff,2),[]),3);
figure;
for i = 1:length(filenameList)
    subplot_tight(3,3,i);
    for j = 1:length(binSelection)-1
        tempIdx= tSortIdx(binSelection(j): binSelection(j+1)-1);
        tempAct = mean(tPSTH( tempIdx,:,i),1) - mean(mean(tPSTH(tempIdx,1:3,i),1),2);
        plot(tempAct); hold on; 
    end
    if i == 1; legend('5-20 perc','20-40 perc','40-60 perc','60-80 perc','80-95 perc');end
    title(shortFilename{i})
end
subplot_tight(3,3,length(filenameList)+1); title('dff')
for j = 1:length(binSelection)-1
    tempIdx= tSortIdx(binSelection(j): binSelection(j+1)-1);
    tempAct = mean(tDffPSTH(tempIdx ,:),1) - mean(mean(tDffPSTH(tempIdx,1:3),1),2);
    plot(tempAct); hold on; 
end
%% plot 2.3 - Correlation of peak activity 

actMatT = cat(2,reshape(tPeakMean,size(tPeakMean,1)*size(tPeakMean,2),size(tPeakMean,3)), tPeakMeanDff(:));
actMatF = cat(2,reshape(fPeakMean,size(fPeakMean,1)*size(fPeakMean,2),size(fPeakMean,3)), fPeakMeanDff(:));
actMat = cat(1,actMatT,actMatF);
figure; subplot(1,2,1);imagesc(corr(actMat)); colorbar; caxis([0.8 1]); title('Correlation of Mean Activity')

actMatT = cat(2,reshape(tPeakTrial,size(tPeakTrial,1)*size(tPeakTrial,2)*size(tPeakTrial,3),[]), tPeakTrialDff(:));
actMatF = cat(2,reshape(fPeakTrial,size(fPeakTrial,1)*size(fPeakTrial,2)*size(fPeakTrial,3),[]), fPeakTrialDff(:));
actMat = cat(1,actMatT,actMatF);
subplot(1,2,2); imagesc(corr(actMat)); colorbar; caxis([0.5 1]); title('Correlation of Trial-by-Trial Activity ')
%% plot 3 - get number of responsive neurons 
tAct = reshape(tPeak,size(tPeak,1),size(tPeak,2),size(tPeak,3)*size(tPeak,4),size(tPeak,5));
fAct = reshape(fPeak,size(fPeak,1),size(fPeak,2),size(fPeak,3)*size(fPeak,4),size(fPeak,5));
peakDetection = false;peakFrames = 1:6;smoothWindow = 0;
for i = 1:length(filenameList)
    for j = 1:size(tAct,1)
        tempAct = squeeze(tAct(j,:,:,i)); [h, tempAuc, tempIdx] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
        ttestT(j,i) = h;rocAucT(j,i) = tempAuc; peakIdxT(j,i) = tempIdx;
        tempAct = squeeze(fAct(j,:,:,i)); [h, tempAuc] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
        ttestF(j,i) = h;rocAucF(j,i) = tempAuc;
    end
end
tActDff = reshape(tPeakDff,size(tPeakDff,1),size(tPeakDff,2),[]);fActDff = reshape(fPeakDff,size(fPeakDff,1),size(fPeakDff,2),[]);
peakDetection = true;peakFrames = 1:12;smoothWindow = 3;
for j = 1:size(tActDff,1)
    tempAct = squeeze(tActDff(j,:,:)); [h, tempAuc, tempIdx] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
    ttestTDff(j) = h; rocAucTDff(j) = tempAuc; peakIdxT(j,i) = tempIdx;
    tempAct = squeeze(fActDff(j,:,:)); [h, tempAuc] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
    ttestFDff(j) = h; rocAucFDff(j) = tempAuc;
end

%% plot 3.1 - 
figure;
for i = 1:length(filenameList)
    h = subplot(2,3,i);
    bothFlag = ((ttestF(:,i)==1)' & ttestFDff==1);
    dffFlag = ((ttestF(:,i)==0)' & ttestFDff==1);
    spkFlag = ((ttestF(:,i)==1)' & ttestFDff==0);
    fn_plotPieMultiPanel({{bothFlag,dffFlag,spkFlag}}, 'legend', {{'both','dff only','spike only'}},...
        'title',shortFilename(i),'axes',h);
end
figure; hold on;
for i = 1:length(filenameList)
    bins = 0.4:0.005:1.0; binsPlot = bins(1:end-1) + 0.005/2;
    %[N, ~,~] = histcounts(rocAucT(:,i), 0.4:0.005:0.8,'Normalization','probability');
    cdfplot([rocAucT(:,i);rocAucF(:,i)]); xlim([0.5 1.0])
end
h = cdfplot([rocAucTDff rocAucFDff] ); xlim([0.5 1.0])
set( h, 'Color', [0 0 0], 'LineWidth', 1.5); legend({shortFilename{:},'dff'},'Location','Best')
title('Distribution of ROC value')

figure;  hold on;
respFlag = ttestTDff == 1 | ttestFDff == 1;
h = cdfplot([rocAucT(respFlag,3);rocAucF(respFlag,3)]); set( h, 'Color', matlabColors(1), 'LineWidth', 1.5);
h = cdfplot([rocAucT(~respFlag,3);rocAucF(~respFlag,3)]);set( h, 'Color', matlabColors(1),'LineStyle','--', 'LineWidth', 1.5);
h = cdfplot([rocAucTDff(respFlag) rocAucFDff(respFlag)]); set( h, 'Color', [0 0 0], 'LineWidth', 1.5);
h = cdfplot([rocAucTDff(~respFlag) rocAucFDff(~respFlag)]); set( h, 'Color', [0 0 0], 'LineStyle','--', 'LineWidth', 1.5);
title('ROC value distribution'); legend('spk-resp','spk-not resp','dff-resp','dff-not resp','Location','Best')
xlabel('AUC value of ROC (measurement of signal-to-noise)');
%% plot 4 - calculate the SI of all neurons
spikeSI = abs(tPeakMean- fPeakMean) ./ (abs(tPeakMean) +  abs(fPeakMean));
spikeSI = reshape(spikeSI,size(spikeSI,1)*size(spikeSI,2),[]);
spikeSI(isnan(spikeSI)) = 0;
spikeSIDff = abs(tPeakMeanDff- fPeakMeanDff) ./ (abs(tPeakMeanDff) +  abs(fPeakMeanDff));
spikeSIDff = spikeSIDff(:);
figure; subplot(1,3,1); hold on;
for i = 1:length(filenameList); cdfplot(spikeSI(:,i));end; h = cdfplot(spikeSIDff);
set( h, 'Color', [0 0 0], 'LineWidth', 1.5); legend({shortFilename{:},'dff'},'Location','Best');
title('Distribution of SI')
subplot(1,3,2); hold on;
for i = 1:length(filenameList); cdfplot(spikeSI(:,i)-spikeSIDff);end
legend(shortFilename,'Location','Best');title('Distribution of deltaSI')
subplot(1,3,3); 
tempMat = cat(2,spikeSI,spikeSIDff);imagesc(corr(tempMat));
title('Correlation of SI'); colorbar; caxis([0.4 1])
%% plot 5 - 



%% functions 
function [h, tempAuc, peakIdx] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow)
    tempBaseline = tempAct(1:3,:); tempBaselineMean = mean(tempBaseline,1);
    if smoothWindow ~= 0; tempMeanAct = smoothdata(mean(tempAct,2),'gaussian',smoothWindow);
    else; tempMeanAct = mean(tempAct,2); end
    toneFrame = 5; % number of pretone frames
    selectToneFrame = toneFrame + peakFrames;
    [~,tempPeakIdx] = max(tempMeanAct(selectToneFrame));
    peakIdx = tempPeakIdx + selectToneFrame(1)-1; 
    if peakDetection; tempPeak = tempAct(peakIdx,:); 
    else; tempPeak = mean(tempAct(selectToneFrame,:),1); end
    [h,p] = ttest(tempPeak,tempBaselineMean,'tail','right');
    % do an roc for each tone 
    rocAct = [tempBaseline(:)' tempPeak]; 
    rocAct(rocAct<1e-4) = randn(1,sum(rocAct<1e-4))*1e-4;% introduce small noise to spks to avoid domination of 0s
    rocLabel = [zeros(1,numel(tempBaseline)) ones(1,length(tempPeak))];
    [tpr, fpr, threshold] = roc(rocLabel, rocAct);
    tempAuc = trapz([0 fpr 1],[0 tpr 1]);
end

function [tPeakMean, fPeakMean, tAct, fAct,tPeakTrial,fPeakTrial] = getPeakAct(behavMatrix,act,peakDetection,peakFrames)
    tFrame = floor(behavMatrix(behavMatrix(:,4)==1 | behavMatrix(:,4)==2,12)/2);
    fFrame = floor(behavMatrix(behavMatrix(:,4)==3 | behavMatrix(:,4)==4,12)/2);
    selectFrame = -5 : 30; toneFrame = abs(selectFrame(1)); selectToneFrame = toneFrame + peakFrames;
    tAct = []; fAct = [];
    for k = 1:length(tFrame)
        tAct(:,:,k) = act(:,tFrame(k) + selectFrame);
        fAct(:,:,k) = act(:,fFrame(k) + selectFrame);
    end
    if peakDetection
        tActMean = mean(tAct,3); tPeakMean = max(smoothdata(tActMean(:,selectToneFrame),2,'gaussian',3),[],2);
        fActMean = mean(fAct,3); fPeakMean = max(smoothdata(fActMean(:,selectToneFrame),2,'gaussian',3),[],2);
        tPeakTrial = squeeze(max(smoothdata(tAct(:,selectToneFrame,:),2,'gaussian',3),[],2));
        fPeakTrial = squeeze(max(smoothdata(fAct(:,selectToneFrame,:),2,'gaussian',3),[],2));
        
    else
        tPeakTrial = squeeze(mean(tAct(:,selectToneFrame,:),2));tPeakMean = mean(tPeakTrial,2);
        fPeakTrial = squeeze(mean(fAct(:,selectToneFrame,:),2));fPeakMean = mean(fPeakTrial,2);
    end
end

function f = fillErrorbarPlot(xdata,yMean, ySEM,varargin)
    x = [xdata, fliplr(xdata)];
    y = [yMean+ySEM, fliplr(yMean-ySEM)];
    f = fill(x,y,varargin{:});
end