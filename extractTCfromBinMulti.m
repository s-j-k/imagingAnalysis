function TC = extractTCfromBinMulti(mouse,channel)

global info;

if strcmp(mouse,'sk78')
%     suite2ppath = ['V:\sjk\sk78\suite2p_sk78_010\'];
% %     sbxpath = ['C:\Users\sjkim1\Desktop\scratch\' mouse  '\'];
%     sbxpath='D:\sk78\SfnPoster\day7\';

    suite2ppath = ['V:\sjk\sk78\027\suite2p\'];
    sbxpath='V:\sjk\sk78\027\';
end

if strcmp(mouse,'sk83')
%     suite2ppath = ['C:\Users\sjkim1\Desktop\scratch\' mouse '\FOV2\suite2p\'];
    suite2ppath=['V:\sjk\ACMGB_2\sk83\tonotopy\FOV1\suite2p\'];
    sbxpath='V:\sjk\ACMGB_2\sk83\tonotopy\FOV1\';
%     sbxpath=['C:\Users\sjkim1\Desktop\scratch\' mouse '\FOV2\'];
%     sbxpath = ['C:\Users\sjkim1\Desktop\scratch\' mouse '\'];
end



if strcmp(mouse,'sk70')
%     suite2ppath = ['C:\Users\sjkim1\Desktop\scratch\' mouse '\suite2p\'];
%     sbxpath = ['C:\Users\sjkim1\Desktop\scratch\' mouse '\'];
    suite2ppath='X:\sjk\ACMGB_1\sk70\tonotopy\sk70_037_002\suite2p_green\';
    sbxpath='X:\sjk\ACMGB_1\sk70\tonotopy\sk70_037_002\';
end

% if strcmp(mouse,'sk82')
% %     suite2ppath = ['C:\Users\sjkim1\Desktop\scratch\' mouse '\FOV2\suite2p\'];
%     suite2ppath=['C:\Users\sjkim1\Desktop\scratch\sk82\FOV1\suite2p\'];
%     sbxpath=['C:\Users\sjkim1\Desktop\scratch\sk82\FOV1\'];
% %     sbxpath=['C:\Users\sjkim1\Desktop\scratch\' mouse '\FOV2\'];
% %     sbxpath = ['C:\Users\sjkim1\Desktop\scratch\' mouse '\'];
% end
% 
% if strcmp(mouse,'sk84')
%     suite2ppath = ['C:\Users\sjkim1\Desktop\scratch\' mouse '\FOV2\suite2p\'];
%     sbxpath=['C:\Users\sjkim1\Desktop\scratch\' mouse '\FOV2\'];
% %     sbxpath = ['C:\Users\sjkim1\Desktop\scratch\' mouse '\'];
% end

if strcmp(mouse,'sk76')
    suite2ppath = ['X:\sjk\ACMGB_1\sk76\suite2p\'];
    sbxpath = ['X:\sjk\ACMGB_1\sk76\'];
end



if strcmp(channel,'red')
    readChannel = 'data.bin';
elseif strcmp(channel,'green')
    readChannel='data_chan2.bin';
else
    error(['channel should be a string that is either red (ch2) or green (ch1)']);
end

cd(sbxpath);
files = dir('*.sbx');
nFiles = length(files);
names = cell(nFiles,1);
for i=1:nFiles, 
    names{i} = files(i).name(1:end-4); 
end
% names=names(1);
% nFiles=1;
% names = names(2); nFiles=1 ;% take only the second one

cd(sbxpath);
nFrames = nan(nFiles,1);
for ii=1:nFiles
    sbxread(names{ii},1,1);
    nFrames(ii) = info.max_idx;       
end

nPlanes = info.otparam(3);
if isnan(nPlanes), nPlanes = 1; end
if nPlanes>1 % Here, the nb of frames / plane is NOT cumulative
    nFrames_oneplane = nan(nFiles,nPlanes);
    nFrames_oneplane(logical(mod(nFrames,2)),:) = [round(nFrames(logical(mod(nFrames,2)))/nPlanes) round(nFrames(logical(mod(nFrames,2)))/nPlanes)-1];
    nFrames_oneplane(~mod(nFrames,2),:) = [nFrames(~mod(nFrames,2))/nPlanes nFrames(~mod(nFrames,2))/nPlanes];
    nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];
else
    nFrames_oneplane = [zeros(1,nPlanes);nFrames];
end

cd(suite2ppath);
TC = cell(nPlanes,1);
for i=1:nPlanes
    cd([suite2ppath '\plane' num2str(i-1)]);    
    data = load('Fall.mat');
    
    % check if nb of frames per plane computed is true in suite2p files
     nFrameThisPlane = size(data.Fneu,2);
    if nFrameThisPlane==sum(nFrames_oneplane(:,i))
%     if nFrameThisPlane==sum(nFrames_oneplane(:,i))
        disp(['Correct nb of frame for plane ' num2str(i) '. Good to go!']);
    else
        error(['Bad nb of frame computed for plane ' num2str(i) ', check it out! check .sbx is the right one']);
    end
    
    ly = data.ops.Ly;
    lx = data.ops.Lx;
    fileID = fopen(readChannel,'r'); % open binary file, 
    
    
%     % pixel based analysis - todo 0802\
%     if grid==0
%         % extract TC from all pixels
%         win = [1 1]; % size spatial window, currently 1 x 1 px
%     %     kall = 0;
%         for j=1:nFiles        
%             
%             % first, split the stack half 
%             nimg = nFrames_oneplane(j+1,i);
%             k=0;
%             blksize = 2000;%7000; % nb of frames loaded at a time (depend on RAM)
%             to_read = min(blksize,nimg-k);  
%             while to_read>0
%                     a = fread(fileID,ly*lx*to_read,'*int16');
%                     A = reshape(a,lx,ly,[]); 
%         %             avgA = mean(A,3); % figure;subplot(2,2,1);imagesc(avgA);colormap gray;
%                     fun = @(block_struct) nanmean(block_struct.data(:));
%                     b = blockproc(A(:,:,1),win,fun);            
%                     tempTC = nan(size(b,1)*size(b,2),size(A,3));
%                     tempTC(:,1) = b(:);
%                     for kk=2:size(A,3)
%                         if kk<(size(A,3)/2)
%                             b = blockproc(A(:,:,kk),win,fun);
%                             tempTC(:,kk) = b(:);
%                         
%                         
%                         else kk>=(size(A,3)/2)
%                             b=blockproc(A(:,:,kk),win,fun);
%                             tempTC(:,kk)=b(:);
%                         end
%                         save([mouse 'tempTC_' kk '.mat'])
%                     end
%                     TC{i} = [TC{i} tempTC]; % figure;plot(mean(TC{1})');
%                     k = k+to_read;
%                     to_read = min(blksize,nimg-k);
%                 end
%                 
%             end
%             disp(['File ' num2str(j) '/' num2str(nFiles) ' DONE! ' num2str(size(TC{i},1)) ' frames so far.']);
%         end
%         fclose all;
    
    % for every pixel you hav ethe TC, then you can go back to the px data
    % and avg across for a given window
    % this way you can average for a given window that you specify after
    % the fact 
    
    
%     if grid==1   
        % extract TC in ROI grid
        win = [20 20]; % size spatial window, currently 30 x 30 px
    %     kall = 0;
        for j=1:nFiles        
            k=0;
            nimg = nFrames_oneplane(j+1,i);
            blksize = 2000;%7000; % nb of frames loaded at a time (depend on RAM)
            to_read = min(blksize,nimg-k);   
            while to_read>0
                a = fread(fileID,ly*lx*to_read,'*int16');
                A = reshape(a,lx,ly,[]); 
    %             avgA = mean(A,3); % figure;subplot(2,2,1);imagesc(avgA);colormap gray;
                
                fun = @(block_struct) nanmean(block_struct.data(:));
                b = blockproc(A(:,:,1),win,fun);            
                tempTC = nan(size(b,1)*size(b,2),size(A,3));
                tempTC(:,1) = b(:);
                for kk=2:size(A,3)
                    b = blockproc(A(:,:,kk),win,fun);
                    tempTC(:,kk) = b(:);
                end
                TC{i} = [TC{i} tempTC]; % figure;plot(mean(TC{1})');
                k = k+to_read;
                to_read = min(blksize,nimg-k);
                for ee=1:size(A,3)
                    imwrite2tif(A(:,:,ee),info,[channel '_Stack' num2str(ee+k-2000) '.tif'],'int16');
                end
            end
            disp(['File ' num2str(j) '/' num2str(nFiles) ' DONE! ' num2str(size(TC{i},1)) ' frames so far.']);
        end
        fclose all;
    end
    
    
    
      
%% save TC if wanted
 % now save the TC for the pixel based analysis 
 % 
fileTC=[ mouse '_tc_' channel '_grid_' num2str(grid) '.mat'];
save(fileTC)
%% get tuning

singleneuronanalysis = true;
savefig = true;
fileName = [mouse '_' channel '_tonotopy'];
cd(sbxpath)
if savefig
    tuningFolderName = ['TuningCurve2_' fileName];
    if exist([sbxpath '/' tuningFolderName],'dir') ~= 7
        mkdir(tuningFolderName);
        mkdir([tuningFolderName '/singleNeuron']);
        mkdir([tuningFolderName '/population']);
    end
end

% Default values
nTones = 17;
nFramesPerTone = 100/nPlanes; % 25
nFramesPerTrial = nFramesPerTone * nTones; % 425
%     startTrial = 3; % get rid of first 2 trials where overall fluor is high
% nTrials = round((nFrames+1)/100/17);
nTrials = round((nFrames+1)/100/17)-1;
startTrial = 1; 
nFrames = nFramesPerTrial*nTrials; % 4250       

tc = [];
%%
% stuck here! 
for p=1:nPlanes
    tc = nancat(1,[TC{1}],[TC{2}]);
end
nROIs = size(tc,1);

TC = tc';
completeTC = nan(sum(nFrames)+100,nROIs); % just so that if TC size is a lot larger, it can still work
completeTC(2:nFrames,:) = TC(1:nFrames-1,:); % shift one frame so that all the tone onset is on 101, 201, 301...
TC2 = completeTC(1:nFrames,:);
%%
allNeuronMean = reshape(TC2,nFramesPerTone,nTones,nTrials,nROIs);
allNeuronMean = allNeuronMean(:,:,startTrial:end,:);
allNeuronMean = squeeze(nanmean(nanmean(nanmean(allNeuronMean,2),3),4));
[~,peakIndex] = max(smoothdata(allNeuronMean,'movmean',4));

% startFrame = peakIndex-ceil(6/nPlanes);
% endFrame = peakIndex+ceil(6/nPlanes);
% peakFrames  = startFrame : endFrame;
% pretoneFrames = length(peakFrames);
pretoneFrames = 10;
completeTCpretone = nan(nFrames+100,nROIs); % just so that if TC size is a lot larger, it can still work
% completeTCpretone(2+pretoneFrames:size(TC2,1)+1+pretoneFrames,:) = TC2; % shift one frame so that all the tone onset is on 101, 201, 301...
completeTCpretone(pretoneFrames:size(TC2,1)+pretoneFrames-1,:) = TC2; % shift one frame so that all the tone onset is on 101, 201, 301...
TCpretoneTemp = completeTCpretone(1:nFrames,:);
    
toneorder = [45255 8000 13454 4757 5657,...
    22627 64000 53817 4000 9514,...
    16000 6727 19027 26909 32000,...
    11314 38055];    
toneindex = [9;4;5;12;2;10;16;3;11;13;6;14;15;17;1;8;7];

TCreshape=reshape(TC2,nFramesPerTone,nTones,nTrials,nROIs); %reshape(tc,25,17,9) 25 = frames per tone, 17 = number of tones, 9 = trials
TCpretoneReshape = reshape(TCpretoneTemp,nFramesPerTone,nTones,nTrials,nROIs);

TCreorder=zeros(size(TCreshape));
TCpretoneReorder = zeros(size(TCpretoneReshape));
for x=1:nTones
    index=toneindex(x);%+1; %the "+1" is because the first 20 frames are 38 kHz that wraps around
    TCreorder(:,x,:,:) = TCreshape(:,index,:,:);
    TCpretoneReorder(:,x,:,:) = TCpretoneReshape(:,index,:,:);
end
TCtone = TCreorder;
TCpretone = TCpretoneReorder;
TCreorder=permute(reshape(TCreorder,nFramesPerTrial,nTrials,nROIs),[2 1 3]); % now TC reorder is the size of trials*frame a trial * neuron

% use the same basline for all trials, different 
% have sem in the tuning curve data
TCreorder = TCreorder(startTrial:end,:,:); % [nTrials-1,(nFramesPerTone*nTones),nROIs]
baseline = prctile(reshape(TCreorder, [(nTrials-startTrial+1)*nFramesPerTone*nTones,nROIs]),25);
TCreorder = TCreorder ./ repmat(reshape(baseline,[1 1 nROIs]),nTrials-startTrial+1,nFramesPerTone*nTones,1);

TCtrialMean = squeeze(nanmean(TCreorder));
TCtrialMedian = squeeze(nanmedian(TCreorder));
TCtrialSEM = squeeze(nanstd(TCreorder)./sqrt(nTrials));

TCtoneMean = reshape(TCtrialMean,[nFramesPerTone,nTones,nROIs]);
TCtoneMedian = reshape(TCtrialMedian,[nFramesPerTone,nTones,nROIs]);

peakFrames = 23:36;
tuningData = reshape(TCreorder, [nTrials-startTrial+1, nFramesPerTone,nTones,nROIs]);
tuningData = reshape(tuningData(:,peakFrames,:,:),[(nTrials-startTrial+1)*length(peakFrames),nTones,nROIs]);
tuningMean = squeeze(nanmean(tuningData));
tuningMedian = squeeze(nanmedian(tuningData));
tuningSEM = squeeze(nanstd(tuningData)/sqrt(size(tuningData,1)));
tones=toneorder(toneindex);


figure;
subplot(2,2,1);hold on;
PlotColorCurves(mean(TCtoneMedian,3)');
set(gca, 'YTick', 1:17);
set(gca, 'YTickLabel', tones);
colormap jet;
title(mouse);
subplot(2,2,3);hold on;
[~,bfs] = max(tuningMean);
spatial_bfs = reshape(bfs',size(b,1),size(b,2));
C = colormap(jet);
indexC = linspace(1,size(C,1),17);
PlotColorMap(spatial_bfs);colormap(jet);
clim([1 17]);
subplot(2,2,2);hold on;
AvTCtoneMean = mean(TCtoneMean,3);
plot(AvTCtoneMean(:));
xlim([1 1700]);
set(gca, 'XTick', 1:100:1700);
set(gca, 'XTickLabel', round(tones/1000));
subplot(2,2,4);hold on;
PlotColorCurves(1:17,[1 17]);colormap(jet);
clim([1 17]);
set(gca, 'XTick', 1:17);
set(gca, 'XTickLabel', round(tones/1000));