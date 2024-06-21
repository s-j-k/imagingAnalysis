% Generate example neuronal traces
% Example of CD036.

mouse = 'CD036';
path = 'M:\Celine\exci_variables\';

cellselection = load([path mouse '-cellselection.mat']);
cellselection = cellselection.cellselection;
takenCells = cellselection{1,1};
findgooddays = cellselection{1,2};
nKeptCells = sum(takenCells);   

nframes_psth = pretone*round(acq) + posttone*round(acq);       
tonepsth = load([path mouse '-tonepsth.mat']);
tonepsth = tonepsth.tonepsth;       
tonepsth = tonepsth(:,1);

trial_activity = tonepsth{1};
figure;plot(squeeze(trial_activity(:,1,:))+ repelem((1:nKeptCells)/10,75,1),'k');
ylim([0 (nKeptCells+1)/10]);

figure;plot(squeeze(trial_activity(:,1,:))+ repelem((1:nKeptCells)/10,75,1),'k');
ylim([0 (nKeptCells+1)/10]);

results = load([path mouse '-results_nosignals.mat']);
results = results.results;
nFrames_oneplane = results{1};
matrix = results{3};
ishere = results{4};
ks = results{5};
ks2 = results{6};
ctx = results{7};
acq = results{8};
TONEF = results{9};
REWF = results{10};
dayss_behav = results{12};
startfrom = results{13};
nDays_behav = length(dayss_behav);
if strcmp(mouse,'cd017')
    nop = matrix(:,BLOC)==4; % remove bad trials
else
    nop = ~ismember(matrix(:,BLOC),matrix(:,BLOC)); % take everything
end
ok = ~nop & ctx(:,1);

global_here = [];
for p=1:nPlanes % loop through planes
    here = ishere{p}(:,dayss_behav+1); 
    here(isnan(here) | here~=1) = 0;
    nCellsHere = sum(here); 
    global_here = [global_here logical(here)']; % days x neurons
end  
[nmaxcell,daywithmaxcells] = max(sum(global_here,2));
% trace from day with maximum cell number

trial_activity = tonepsth{daywithmaxcells}(:,:,logical(global_here(daywithmaxcells,:)));
figure;plot(squeeze(trial_activity(:,1,:))+ repelem(1:nmaxcell,75,1));



sigs = {signals{1}(:,:,3),signals{2}(:,:,3)};
taken = round(5*15.63);
start = 2000;
s = [sigs{1}(:,start:start+taken);sigs{2}(:,start:start+taken)];
figure;plot((s(1:200,:) + repelem((1:200)'/10,1,taken+1))');
ylim([0 201/10]);
