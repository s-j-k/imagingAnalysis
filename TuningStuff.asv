tmpNorm=TCpretone_reorder-median(TCpretone_reorder(1:10,:,:));
tmpNormFull=tmpNormFull';stdTmpNorm=std(tmpNormFull);
tmpNormFull=tmpNormFull';tmpNormFull=reshape(tmpNorm,[30000,4067]);
for ii=1:size(tmpNormFull,1)
meanTmpNorm(ii)=mean(tmpNormFull(ii,1:size(tmpNormFull,2)));
if max(tmpNormFull(ii,1:size(tmpNormFull,2)))>stdTmpNorm(ii)*20
iscell(ii)=1;
else
iscell(ii)=0;
end
end
find(iscell==0);
iscell=logical(iscell);
imagesc(squeeze(mean(TCpretone_reorder(:,:,iscell),1))');
tmpReshaped=reshape(TCpretone_reorder,[100,15,20,size(TCpretone,1)]);
tmpTrialsMeanReshaped=squeeze(mean(tmpReshaped,2));
tmpAllCellsReshape=squeeze(mean(tmpTrialsMeanReshaped,3))';

%%
% % need to add the green channel manually first
greenPath='W:\su\CODE\imagingAnalysis\process2p-master\examples\config\20230704_sk132_tuning_00001_Tuning\population\';
cd(greenPath)
load([greenPath 'tuning.mat']);
greenSuite2pPath='E:\sk132\tuning\suite2p\plane0';
cd(greenSuite2pPath);
load('Fall.mat');
%suite2p data
greenF=F;greenFneu=Fneu;greenIsCell=iscell;greenRedCell=redcell;greenOps=ops;greenSpks=spks;greenStat=stat;
% for the extracted data
greenNeuronEachPlane=neuronEachPlane;greenPeakFrames=peakFrames;greenTuning=tuning;
greenTCpretone_reorder=TCpretone_reorder(:,:,:,greenTuning.responsiveCellFlag);
greenToneLabel=toneLabel;
sum(greenTuning.responsiveCellFlag(:)==1) % how many responsive axon segments 

cd('O:\sjk\sk132\tuning_red\suite2p\plane0')
% cd('O:\sjk\sk126\tuning\sk126_0708_tuning2\suite2p_red\plane0');
load('Fall.mat')
load('O:\sjk\sk132\tuning_red\20230704_tuning_00001_Tuning1\population\tuning.mat')
% load('O:\sjk\sk126\tuning\sk126_0708_tuning2\20230708_sk126_00002_Tuning\population\tuning.mat');

%%
% do some preprocessing

tone = [64001 64000 53817 4000 9514 ...
16000 6727 64002 19027 26909 32000 ...
64003 11314 38055 45255 8000 13454 4757 5657 22627];
% 64001 is white noise, 4264 upsweep, 6424 downsweep
toneLabel = strsplit(int2str(round(tone)));

[order,toneindex] = sort(tone);

pureTone = [64000 53817 4000 9514 ...
16000 6727 19027 26909 32000 ...
11314 38055 45255 8000 13454 4757 5657 22627];
% 64001 is white noise, 4264 upsweep, 6424 downsweep
pureToneLabel = strsplit(int2str(round(pureTone)));
[pureToneOrder,pureToneIndex] = sort(pureTone);pureToneOrder=pureToneOrder';
nTones=length(pureTone);

%%
% remove the 2 sweeps and whitenoise
numcells=size(TCpretone_reorder,4);
redTmp=reshape(TCpretone_reorder,[25,20*10,numcells]);
redTmpNorm=redTmp-median(redTmp(1:2,:,:)); % subtract baseline
redTmpNormRe=reshape(redTmpNorm,[25,20,10,numcells]);
truncatedReorderTC=redTmpNormRe(:,1:17,:,:); % this is the whole TC for all reps for tones 1:17
%redo the baseline calculation
gnumcells=size(greenTCpretone_reorder,4);
tmp = reshape(greenTCpretone_reorder,[100,20*15,gnumcells]);
tmpnorm = tmp-median(tmp(1:10,:,:)); % subtract baseline
tmpnormre = reshape(tmpnorm,[100,20,15,gnumcells]); % reshape
% plot to see what it looks like
figure; plot(squeeze(mean(median(tmpnormre(:,:,:,:),3),2)))
gTruncatedReorderTC=tmpnormre(:,1:17,:,:); 
%truncatedReorderTC 100 frames, 17 tones, 15 repetitions, 77 cells 
medRepetition=squeeze(median(truncatedReorderTC,3));
gMedRepetition=squeeze(median(gTruncatedReorderTC,3));
%what does it look like?
figure;plot(squeeze(mean(gMedRepetition,3)),'Color',[0.8 0.8 0.8]);title('sk132 Green MGB Mean tone evoked response for each Tone');
hold on; plot(mean(mean(gMedRepetition,3),2),'LineWidth',2);xline(10);ylim([-0.5 3]);
figure;plot(squeeze(mean(medRepetition,3)),'Color',[0.8 0.8 0.8]);title('sk132 Red AC Mean tone evoked responses');
hold on;plot(mean(mean(medRepetition,3),2),'LineWidth',2);xline(8);

figure;
imagesc(squeeze(mean(medRepetition(10:40,:,:)))'); % mean 
title('sk 132 Mean response Pure Tones Red AC Cells');ylabel('Cell Number');xlabel('Tone Label');
xticks([1:17]);xticklabels(pureToneLabel(pureToneIndex));xlabel('Frequency(kHz)')
figure;
imagesc(squeeze(mean(gMedRepetition(10:50,:,:)))');
ylabel('Cell Number');xlabel('Tone Label');title('sk132 Mean response Pure Tones Green MGB Boutons');
xticks([1:17]);xticklabels(pureToneLabel(pureToneIndex));xlabel('Frequency(kHz)')

%%
% get the bf of all of the cell 
[bfValue,bfPureTones]=max(squeeze((mean(medRepetition(10:25,:,:)))));
% [gBfValue,gBfPureTones]=max(squeeze((mean(gMedRepetition(20:50,:,:)))));
%get all of the neurons with the same BF and get their index
% get their tuning curve based on the index and mean them together
% plot the meaned traces
% normalize by peak of response? so everything is 0 to 1 
% to do this, divide by the max amplitude per cell 
for pp = 1:(nTones-3)
   tempBFIdx=find(bfPureTones==pp);
%    gTempBFIdx=find(gBfPureTones==pp);
   meanBFTC(pp,:)=max(squeeze(mean(medRepetition(10:25,:,tempBFIdx),3))); 
%    gMeanBFTC(pp,:)=max(squeeze(mean(gMedRepetition(20:50,:,gTempBFIdx),3)));
%    meanBFTC(pp,:)=mean(squeeze((mean(medRepetition(1:30,:,tempBFIdx)))),2); 
end
normMeanBFTC=(meanBFTC-min(meanBFTC')')./(max(meanBFTC')-min(meanBFTC'))'; % this uses the min and the max
% gNormMeanBFTC=(gMeanBFTC-min(gMeanBFTC')')./(max(gMeanBFTC')-min(gMeanBFTC'))';
% normalizes 0 to 1 bc jenni is v smart
%% plot it
C = colormap('jet'); 
colormapIndex = round(linspace(1,size(C,1),nTones-3));
tempNTones = 17; jjet = jet; 
colormapIndex = round(linspace(1,size(jjet,1),tempNTones)); 
toneColor = jjet(colormapIndex,:);

figure; hold on;qq=1;
tonesNoNan=find(sum(isnan(normMeanBFTC),2)==0);
for qq = 1:length(tonesNoNan)%1:(nTones-3)
plot(meanBFTC(tonesNoNan(qq),:)','LineWidth',1.5,'color',toneColor(tonesNoNan(qq),:))
      %     plot(smooth(normMeanBFTC(qq,:),12)','LineWidth',1.5,'color',toneColor(qq,:))
%     shadedErrorBeditar(smooth(normMeanBFTC(qq,:),7)','LineWidth',1.5,'color',toneColor(qq,:))
end
ylim([0 15]);
xlim([1 17]);
xticks([1:17]);
xticklabels(pureToneLabel(pureToneIndex));xlabel('Frequency(kHz)')
ylabel('dF/F');title('SK132 (Red, AC) Mean Tuning Curves Grouped by BF')
legend(num2str(pureToneOrder(tonesNoNan)));
% now do it again but normalized
qq=1;figure; hold on;
for qq = 1:length(tonesNoNan)%1:(nTones-3)
    plot(normMeanBFTC(tonesNoNan(qq),:)','LineWidth',1.5,'color',toneColor(tonesNoNan(qq),:))
end
ylim([0 1.2]);
xlim([1 17]);
xticks([1:17]);
xticklabels(pureToneLabel(pureToneIndex));xlabel('Frequency(kHz)')
ylabel('dF/F');title('SK132 (Red, AC) Normalized Mean Tuning Curves Grouped by BF')
legend(num2str(pureToneOrder(tonesNoNan)));

%%
%for MGB boutons

figure; hold on;qq=1;
gTonesNoNan=find(sum(isnan(gNormMeanBFTC),2)==0);
for qq = 1:length(gTonesNoNan)%1:(nTones-3)
plot(gMeanBFTC(gTonesNoNan(qq),:)','LineWidth',1.5,'color',toneColor(gTonesNoNan(qq),:))
      %     plot(smooth(normMeanBFTC(qq,:),12)','LineWidth',1.5,'color',toneColor(qq,:))
%     shadedErrorBeditar(smooth(normMeanBFTC(qq,:),7)','LineWidth',1.5,'color',toneColor(qq,:))
end
ylim([0 15]);
xlim([1 17]);
xticks([1:17]);
xticklabels(pureToneLabel(pureToneIndex));xlabel('Frequency(kHz)')
ylabel('dF/F');title('SK132 (green, MGB) Mean Tuning Curves Grouped by BF')
legend(num2str(pureToneOrder(gTonesNoNan)));
% now do it again but normalized
qq=1;figure; hold on;
for qq = 1:length(gTonesNoNan)%1:(nTones-3)
    plot(gNormMeanBFTC(gTonesNoNan(qq),:)','LineWidth',1.5,'color',toneColor(gTonesNoNan(qq),:))
end
ylim([0 1.2]);
xlim([1 17]);
xticks([1:17]);
xticklabels(pureToneLabel(pureToneIndex));xlabel('Frequency(kHz)')
ylabel('dF/F');title('SK132 (Green, MGB) Normalized Mean Tuning Curves Grouped by BF')
legend(num2str(pureToneOrder(gTonesNoNan)));
%%
nFramesPerTrial=1700;
nFramesPerTrial=25*17;
trialMedian = medRepetition; nTones=17;
% trialMedianTrialTC = reshape(trialMedian,[nFramesPerTrial,nNeuron]);
% frameRate=17;
nNeuron=size(TCpretone_reorder,4);
frameRate=7;
trialMedianTrialTC = reshape(trialMedian,[nFramesPerTrial,nNeuron]);

nPlanes=1;
figure; hold on
for ii=1:numcells
%     tuningFig = figure('visible','off');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);% Enlarge figure to full screen.
        timeAxis = (0:(nTones*100)-1) * nPlanes / frameRate;
%         timeAxis = (0:(nTones*25)-1) * nPlanes / frameRate; % by the modulus value
        tempNTones = 17; jjet = jet; colormapIndex = round(linspace(1,size(jjet,1),20)); 
        toneColor= jjet(colormapIndex,:);
        for k = 1:tempNTones
%             tempIdx =((k-1)*100+1):(k*100);
            tempIdx =((k-1)*25+1):(k*25);
            p1 = plot(timeAxis(tempIdx), trialMedianTrialTC(tempIdx,ii),'LineWidth',1.5,'color',toneColor(k,:));
        end
end


% figure; hold on;qq=1;
% for qq = 1:8
%       plot(normMeanBFTC(qq,:)','LineWidth',1.5,'color',toneColor(qq,:))
% end
% ylim([-0 1.2]);xlim([1 17]);
% labelsIndex=pureToneLabel(pureToneIndex);
% xticklabels(labelsIndex(1:8));xlabel('Frequency(kHz)')
% ylabel('dF/F');title('SK132 (Red, AC) Mean Activity Across Repetitions')
% legend(char(pureToneLabel(pureToneIndex)));

%%
% hold off;
% cellIndex=25;
% tempTraceTC=mean(TC_trial(1:tempSampleRate:nFramesPerTrial,1:15,cellIndex),2);
% plot(timeAxis,tempTraceTC,'color',[0.4 0.4 0.4], 'LineWidth',0.5);

%% trying to make the ROI plot and failing (: 

load('F:\sk134\10_10\suite2p\tuning_workspace.mat')
%%
iscellFlag = iscell(:,1);
tempRoi=stat(logical(iscellFlag))';
j=1; % planes 
roisCoord = cell(1,j); % only 1 functional channel for suite2p
for k = 1:length(tempRoi) % neuron in each plane 
    bound = boundary(double(tempRoi{k}.xpix)', double(tempRoi{k}.ypix)',1); % restricted bound
    tempCoord = [tempRoi{k}.xpix(bound)' tempRoi{k}.ypix(bound)'];
    roisCoord{k} = tempCoord;
end
roisBound=roisCoord';
refImg=ops.meanImg;
rotatedrefImg=imrotate(refImg,-180);  
rotatedrefImg=flip(rotatedrefImg,1);
rotatedrefImg=flip(rotatedrefImg,2);

%%
figure;
imagesc(rotatedrefImg);colormap gray;hold on;
ylim([0 size(rotatedrefImg,1)]);xlim([0 size(rotatedrefImg,2)]);
    for j = 1:length(tempRoi)
        if ~isempty(roisBound{j})
            x = roisBound{j}(:,1); %freehand rois have the outlines in x-y coordinates
            y = roisBound{j}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
            if tuning.responsiveCellFlag(j)
                patch(x,y,C(colormapIndex(bfPureTones(j)),:),'EdgeColor','none');
                % flip horizontally & rotate 90 degrees to the right
            else
                patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
            end
        end
    end
    
%% let's ignore that for now...
% make heat map with normalized activity of all of the neurons
% first normalize the activity for each cell (fuck)
[~,sortedBFTones]=sort(bfPureTones);
[~,gSortedBFTones]=sort(gBfPureTones);
evokedResponseTuning=squeeze(mean(medRepetition(20:50,:,sortedBFTones))); 
normResponseTuning=(evokedResponseTuning-min(evokedResponseTuning))./(max(evokedResponseTuning)-min(evokedResponseTuning)); 

gEvokedResponseTuning=squeeze(mean(gMedRepetition(20:50,:,gSortedBFTones)));
gNormResponseTuning=(gEvokedResponseTuning-min(gEvokedResponseTuning))./(max(gEvokedResponseTuning)-min(gEvokedResponseTuning));
% now do the sort & plot part
figure;hold on; 
% medRepetition=squeeze(median(truncatedReorderTC,3));
% imagesc(squeeze(mean(normMedRepetition(20:50,:,sortedBFTones)))'); % mean 
imagesc(flipud(normResponseTuning')); % mean 
ylabel('ROI'); caxis([0,1]); hold on; colormap(flipud(gray));title('sk132 FOV2 Red, AC Normalized Mean Evoked Response to Pure Tones')
xticks([1:17]);xticklabels(pureToneLabel(pureToneIndex));xlabel('Frequency(kHz)');
axis tight 
colorbar

figure;hold on; 
% medRepetition=squeeze(median(truncatedReorderTC,3));
% imagesc(squeeze(mean(normMedRepetition(20:50,:,sortedBFTones)))'); % mean 
imagesc(flipud(gNormResponseTuning')); % mean 
ylabel('ROI'); caxis([0,1]); hold on; colormap(flipud(gray));title('sk132 Green, MGB Boutons Normalized Mean Evoked Response to Pure Tones')
xticks([1:17]);xticklabels(pureToneLabel(pureToneIndex));xlabel('Frequency(kHz)');
axis tight 
colorbar
%%
%now we make the correlation... across all cells
%
figure;
coeffMatrix=corrcoef(normResponseTuning);
imagesc(coeffMatrix); title('sk132 Red, AC Normalized Tuning Correlation matrix');

figure;gCoeffMatrix=corrcoef(gNormResponseTuning);
imagesc(gCoeffMatrix); title('sk132 Green, MGB Normalized Tuning Correlation matrix');

% nwo get the distance between pairs of cells

%stat variable has the ROIs

isCellROI=stat(iscell==1);
% isCellROI=stat;
% toneColorMap=colormap(jet(17));C = colormap('jet');
toneColorMap=colormap(jet(20));C = colormap('jet');
% jjet = jet; colormapIndex = round(linspace(1,size(jjet,1),17)); 
jjet = jet; colormapIndex = round(linspace(1,size(jjet,1),20)); 
toneColor = jjet(colormapIndex,:);
colormapIndex = round(linspace(1,size(C,1),nTones));colormap gray;
figure;imagesc(ops.meanImg);hold on; title('sk132 Red, AC Cell ROI Centroids');
for jj=1:size(normResponseTuning,2)
% for jj=1:length(isCellROI)
    centerYROI(jj)=mean(isCellROI{jj}.ypix);
    centerXROI(jj)=mean(isCellROI{jj}.xpix);
    plot(centerXROI(jj),centerYROI(jj),'.','MarkerSize',20,'Color',C(colormapIndex(bfPureTones(jj)),:));
%     plot(centerXROI(jj),centerYROI(jj),'.','MarkerSize',15,'Color',[0.7 0.7 0.7]);
end
colormap gray;

figure;imagesc(greenOps.meanImg);hold on; title('sk132 Green, MGB Bouton ROI Centroids');
% gIsCellROI=greenStat(greenIsCell==1); % for all cells
gIsCellROI=gIsCellROI(greenTuning.responsiveCellFlag);
for jj=1:size(normResponseTuning,2)
% for jj=1:size(gIsCellROI,2) % for all cells
    gCenterYROI(jj)=mean(gIsCellROI{jj}.ypix);
    gCenterXROI(jj)=mean(gIsCellROI{jj}.xpix);
    plot(gCenterXROI(jj),gCenterYROI(jj),'.','MarkerSize',10,'Color',C(colormapIndex(gBfPureTones(jj)),:));
%     plot(gCenterXROI(jj),gCenterYROI(jj),'.','MarkerSize',10,'Color',[.7 .7 .7]); % for all cells
end


%now find the Euclidean Distance between each cell
dumbind = [];distE=[];
otherdumbind =1;
for ii = 1:length(centerXROI)
    dumbind = 1:length(centerXROI);
%     dumbind(ii)=[]; % because using triangle friend tril() function
    for jj=1:length(dumbind)
        distance=sqrt((centerXROI([dumbind(jj)])-centerXROI([ii]))^2+(centerYROI([dumbind(jj)])-centerYROI([ii]))^2);
        distE(otherdumbind) = distance;
        otherdumbind =otherdumbind+1;
        %ii,jj)=distance;
    end
end

% 
gind = [];gdistE=[];ii=1;otherdumbind =1;
for ii = 1:length(gCenterXROI)
    gind = 1:length(gCenterXROI);
%     gind(ii)=[];
    for jj=1:length(gind)
        gDistance=sqrt((gCenterXROI([gind(jj)])-gCenterXROI([ii]))^2+(gCenterYROI([gind(jj)])-gCenterYROI([ii]))^2);
        gDistE(otherdumbind)=gDistance;
        otherdumbind =otherdumbind+1;
    end
end


figure;imagesc(distE); % should be zero along the inflection point...
xlabel('Cell Number');
ylabel('Cell Number');
title('SK132 (Red, AC) Pairwise Cell Distance (px) Correlation');colorbar;

figure;imagesc(gDistE)
xlabel('Cell Number');
ylabel('Cell Number');
title('SK132 (Green, MGB) Pairwise Cell Distance (px) Correlation');colorbar;

%take out the entries that are comparing the same cell to itself
plotDistE=distE;
gPlotDistE=gDistE;

% plotDistE(tril(true(size(plotDistE)),-1))=[];
% plotCoefMatrix=coeffMatrix;
% plotCoefMatrix(tril(true(size(plotCoefMatrix)),-1))=[];
% plotCoefMatrix=plotCoefMatrix(plotCoefMatrix~=1);

% plotDistE=distE(distE~=0);

% gPlotDistE=gDistE(gDistE~=0);gPlotCoefMatrix=gCoeffMatrix(gCoeffMatrix~=1);
% 
% [v, w] = unique( plotDistE, 'stable' );
% duplicate_indices = setdiff( 1:numel(plotDistE), w );
% plotDistE(duplicate_indices)=[];
% [v, w] = unique( plotCoefMatrix, 'stable' );
% duplicate_indices = setdiff( 1:numel(plotCoefMatrix), w );
% plotCoefMatrix(duplicate_indices)=[];
% distCorrModel=fitlm(plotDistE,plotCoefMatrix);

% %maybe change this to also use triangle friend
% [gv, gw] = unique(gPlotDistE, 'stable' );
% gDuplicate_indices = setdiff( 1:numel(gPlotDistE), gw );
% gPlotDistE(gDuplicate_indices)=[];
% [gv, gw] = unique(gPlotCoefMatrix, 'stable' );
% gDuplicate_indices = setdiff( 1:numel(gPlotCoefMatrix), gw );
% gPlotCoefMatrix(gDuplicate_indices)=[];
% gDistCorrModel=fitlm(gPlotDistE,gPlotCoefMatrix);

% figure;scatter(plotDistE,plotCoefMatrix);
% xlabel('Distance in px');ylabel('Tuning Curve Correlation');
% title('SK132 (Red, AC) Pairwise Cell Distance by Tuning Curve Correlation');
% 
% figure;scatter(gPlotDistE,gPlotCoefMatrix);
% xlabel('Distance in px');ylabel('Tuning Curve Correlation');
% title('SK132 (Green, MGB) Pairwise Cell Distance by Tuning Curve Correlation');
% 
% % plot it
% figure;distCorrPlot=plot(distCorrModel,'Marker','o','MarkerSize',5,'MarkerFaceColor',[122 197 205]./255);
% xlabel('Distance in px');
% ylabel('Tuning Curve Correlation');
% title('SK132 (Red, AC) Pairwise Cell Distance by Tuning Curve Correlation');
% 
% figure;gDistCorrPlot=plot(gDistCorrModel,'Marker','o','MarkerSize',5,'MarkerFaceColor',[122 197 205]./255);
% xlabel('Distance in px');
% ylabel('Tuning Curve Correlation');
% title('SK132 (Green, MGB) Pairwise Cell Distance by Tuning Curve Correlation');
%%
% correlation plot is not really informative 
% first compute the difference in BF in frequency
bfPureTonesFreq=pureToneOrder(bfPureTones);ee=1;
otherdumbind =1;
for ee = 1:length(bfPureTonesFreq)
    smartind = 1:length(bfPureTonesFreq);
    %     smartind(ee)=[];
    for rr=1:length(smartind)
        freqDelta=abs(log2(bfPureTonesFreq([smartind(rr)])/bfPureTonesFreq([ee])));
        allFreqDelta(otherdumbind) = freqDelta;%(ee,smartind(rr))=freqDelta;
        otherdumbind =otherdumbind+1;
    end
end
plotFreqDelta=allFreqDelta;
% there are zeroes that repeat throughout because some cells have the same BF
% plotFreqDelta(tril(true(size(plotFreqDelta)),-1))=[]; % this should work
% reshape(plotFreqDelta,[
% but gives weird problems in the order
% for tt=1:size(allFreqDelta,1)
%     plotFreqDelta(tt,:)=allFreqDelta(tt,;

% use identity matrix
% plotFreqDelta(logical(eye(size(plotFreqDelta))))=[];

gBfPureTonesFreq=pureToneOrder(gBfPureTones);ee=1;
otherdumbind =1;rr=1;
for ee = 1:length(gBfPureTonesFreq)
    smartind = 1:length(gBfPureTonesFreq);
    %     smartind(ee)=[];
    for rr=1:length(smartind)
        gFreqDelta=abs(log2(gBfPureTonesFreq([smartind(rr)])/gBfPureTonesFreq([ee])));
        gAllFreqDelta(otherdumbind) = gFreqDelta;%(ee,smartind(rr))=freqDelta;
        otherdumbind =otherdumbind+1;
    end
end
gPlotFreqDelta=gAllFreqDelta;
%%
figure;distFreqCorrModel=fitlm(plotDistE,plotFreqDelta);
plot(distFreqCorrModel);ylabel('Difference in BF in Octave');xlabel('Distance between pairs of ROIs');
title('sk126 Red, AC Difference in BF by Distance between pairs of ROIs');

figure;gDistFreqCorrModel=fitlm(gPlotDistE,gPlotFreqDelta);
plot(gDistFreqCorrModel);
ylabel('Difference in BF in Octave');xlabel('Distance between pairs of ROIs');
title('sk132 Green, MGB Difference in BF by Distance between pairs of ROIs');
%%
%%now for cells and axon ROIs


%now find the Euclidean Distance between each cell
dumbind = [];bothDistE=[];
otherdumbind =1;
for ii = 1:length(centerXROI)
    dumbind = 1:length(centerXROI);
%     dumbind(ii)=[]; % because using triangle friend tril() function
    for jj=1:length(dumbind)
        bothDistance=sqrt((centerXROI([dumbind(jj)])-gCenterXROI([ii]))^2+(centerYROI([dumbind(jj)])-gCenterYROI([ii]))^2);
        bothDistE(otherdumbind) = bothDistance;
        otherdumbind =otherdumbind+1;

    end
end


for ee = 1:length(gBfPureTonesFreq)
    smartind = 1:length(gBfPureTonesFreq);
    %     smartind(ee)=[];
    for rr=1:length(smartind)
        bothFreqDelta=abs(log2(bfPureTonesFreq([smartind(rr)])/gBfPureTonesFreq([ee])));
        bothAllFreqDelta(otherdumbind) = bothFreqDelta;%(ee,smartind(rr))=freqDelta;
        otherdumbind =otherdumbind+1;
    end
end
bothPlotFreqDelta=bothAllFreqDelta;

%% troubleshooting 
% shortDist=find(plotDistE<100);
% shortDistOctave=plotFreqDelta(shortDist);
% ans=mean(shortDistOctave)
% medDist=find(plotDistE>1000);medDistOctave=plotFreqDelta(medDist);
% ans=mean(medDistOctave)
% 
% [centerXROI(2) centerYROI(2)]
% [centerXROI(3) centerYROI(3)]
% 
% bfPureTonesFreq
