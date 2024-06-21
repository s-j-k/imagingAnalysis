global nPlanes

% -------------- Label of Sounds in the Protocol ---------------------
% ONLY toneLabel and toneindex matters here in this function
% Toneindex sort the tones according to whatever way you want
% E.g here, the 9th tone is placed at position 1 in the new matrix
%     and the 4th tone is placed at position 2, etc.
% ToneLabel is the label that is displayed on graphs.
% Tone does not matter, here I use the tones to create the labels
%     of toneLabel of a cell array of strings
tone = [22627 64001 64000 53817 4000 9514 ...
16000 6727 64002 19027 26909 32000 ...
64003 11314 38055 45255 8000 13454 4757 5657];
% 64001 is white noise, 4264 upsweep, 6424 downsweep
toneLabel = strsplit(int2str(round(tone)));
[order,toneindex] = sort(tone);

% ----------------- Tone Presentation Protocol ---------------------
nTones = length(tone);
nTrials = 10;
nFramesPerTone = 150/nPlanes; % 50
nFramesPerTrial = nFramesPerTone * nTones; % 850
startTrial = 1; % the first tone is on frame 0
nFrames = nFramesPerTrial*nTrials; % 4250
frameRate = 42.77/nPlanes; % in the future, do not hard code this.
% ----------------- Tone Presentation Protocol ---------------------
pretoneFrames = 10;
baselineFrames = 5;

toneOnset = 20/nPlanes;
peakFrameBin = ceil(0.66 * frameRate); %how you select peak activity
% within 660ms after tone onset

% --------------- Method of Calculating Dff and Smoothing Data-----------------------
% If smoothArg NOT DECLARED, TC will NOT be smoothed (e.g. DECONNVOLVED SPIKE DATA)
smoothWindow = 5; smoothArg = {'gaussian',smoothWindow};
% If dffArg is NOT DECLARED, Dff will NOT be calculated (e.g. DECONNVOLVED SPIKE DATA)
% IF dffArg = {}, then the default method of dff is used (rolling median of 1000)
dffArg = {'method', 'movMean','dffWindow',2000,'baselineCorrectionPostDff',...
    true,'baselineCorrectionWindow',2000};

% --------------- Method for Determing Tone-evoked Act-----------------------
toneActSel  = 'peak';
toneActSumBin = 1:ceil(0.33 * frameRate);

% -------------------------- Selection of Analysis-----------------------
% Note that the FIRST item in significant test will be used to determine
% reponsiveness of the neuron
significanceTestList = {'signrank','ttest','anova'};
significanceTestAlpha = {0.01,0.01,0.05};
snrAnalysisList = {}; %snrAnalysisList = {'roc'};
otherAnalysisList = {'suppression'};
% ---------------------------- Selection of Figure-----------------------
saveSingleNeuronFlag = true;