%% List of mice and experiments
mice = {'cd017','cd036','cd019','cd042','cd041','cd044','cd037','zz033'};
exp = [1;1;2;1;2;1;1;2]; % 1=GNG, 2= PASSIVE

%%
m=6;
mouse = mice{m};
path = 'M:\Celine\exci_variables\';
    
% load mouse signal and behavioral info
results = load([path mouse '-results_nosignals.mat']);
results = results.results;
nFrames_oneplane = results{1};
signals = results{2};
matrix = results{3};
ishere = results{4};
ks = results{5};
ks2 = results{6};
ctx = results{7};
acq = results{8};
TONEF = results{9};
REWF = results{10};
acq = results{11};
dayss_behav = results{12};
startfrom = results{13};
nDays_behav = length(dayss_behav);
if strcmp(mouse,'cd017')
    nop = matrix(:,BLOC)==4; % remove bad trials
else
    nop = ~ismember(matrix(:,BLOC),matrix(:,BLOC)); % take everything
end

ok = ~nop & ctx(:,1);
    
% load tone PSTH (DF/F)
tonepsth = load([path mouse '-tonepsth.mat']);
tonepsth = tonepsth.tonepsth;
tonepsth_reinf = tonepsth(:,1);

wrapper = @(x) permute(x,[2 1 3]);
perm_tonepsth = cellfun(wrapper,tonepsth_reinf(:,1),'uniformoutput',false);
m_tonepsth = cellfun(@nanmean,perm_tonepsth,'uniformoutput',false);
    
nCells = size(tonepsth_reinf{1},3);
nDays = size(tonepsth_reinf,1);

fig=figure;
for day=1:nDays 
    subplot(1,nDays,day);hold on;
    PlotColorCurves(squeeze(m_tonepsth{day})');
    clim([0 0.5]);
end


fig=figure;
for day=1:nDays 
    subplot(1,nDays,day);hold on;
    plot(squeeze(m_tonepsth{day}(:,:,13)),'k');
    ylim([-0.02 0.14]);    
end



figure;plot(tonepsth_reinf{day}(:,4,13),'k');

fa_color = [255 159 25]/255; % orange
cr_color = [255 207 25]/255; % yellow
m_color = [47 122 182]/255;
h_color = [47 196 182]/255;
outcome_colors = {h_color,m_color,fa_color,cr_color};
day = 1;
ylimm = [-0.2 1.5];
figure;
a=1;
c=13;
for trial = 55:59
    subplot(2,5,a);
    trialtype = matrix(ok & matrix(:,DAY)==day+1,RESP);
    outcome = trialtype(trial);
    plot(tonepsth_reinf{day}(:,trial,c),'color',outcome_colors{outcome});
    ylim(ylimm);hold on;PlotHVLines(16,'v','k');
    title(['Cell ' num2str(c) ', trial ' num2str(trial)]);
    a=a+1;
end

c=526;
for trial = 55:59
    subplot(2,5,a);
    trialtype = matrix(ok & matrix(:,DAY)==day+1,RESP);
    outcome = trialtype(trial);
    plot(tonepsth_reinf{day}(:,trial,c),'color',outcome_colors{outcome});
    ylim(ylimm);hold on;PlotHVLines(16,'v','k');
    title(['Cell ' num2str(c) ', trial ' num2str(trial)]);
    a=a+1;
end

day = 15;
ylimm = [-0.2 1.5];
figure;
a=1;
c=13;
for trial = 55:59
    subplot(2,5,a);
    trialtype = matrix(ok & matrix(:,DAY)==day+1,RESP);
    outcome = trialtype(trial);
    plot(tonepsth_reinf{day}(:,trial,c),'color',outcome_colors{2});
    ylim(ylimm);hold on;PlotHVLines(16,'v','k');
    title(['Cell ' num2str(c) ', trial ' num2str(trial)]);
    a=a+1;
end

c=526;
for trial = 55:59
    subplot(2,5,a);
    trialtype = matrix(ok & matrix(:,DAY)==day+1,RESP);
    outcome = trialtype(trial);
    plot(tonepsth_reinf{day}(:,trial,c),'color',outcome_colors{outcome});
    ylim(ylimm);hold on;PlotHVLines(16,'v','k');
    title(['Cell ' num2str(c) ', trial ' num2str(trial)]);
    a=a+1;
end



figure;plot(tonepsth_reinf{day}(:,55,c),'k');
ylim(ylimm);hold on;PlotHVLines(16,'v','k');
title(['Cell ' num2str(c) ', trial ' num2str(55)]);

figure;plot(tonepsth_reinf{day}(:,56,c),'k');
ylim(ylimm);hold on;PlotHVLines(16,'v','k');
title(['Cell ' num2str(c) ', trial ' num2str(56)]);

figure;plot(tonepsth_reinf{day}(:,57,c),'k');
ylim(ylimm);hold on;PlotHVLines(16,'v','k');
title(['Cell ' num2str(c) ', trial ' num2str(57)]);

figure;plot(tonepsth_reinf{day}(:,58,c),'k');
ylim(ylimm);hold on;PlotHVLines(16,'v','k');
title(['Cell ' num2str(c) ', trial ' num2str(58)]);

figure;plot(tonepsth_reinf{day}(:,59,c),'k');
ylim(ylimm);hold on;PlotHVLines(16,'v','k');
title(['Cell ' num2str(c) ', trial ' num2str(59)]);

for day=1:nDays       
    for c=13:nCells        
        fig=figure;
%         plot(tonepsth_reinf{day}(:,1:5,c),'k');
        PlotColorCurves(tonepsth_reinf{day}(:,1:100,c)');
        clim([0 1]);
        PlotHVLines(15,'v','w');
        title(['Cell ' num2str(c)]);
        pause();
        close(fig);
    end
end