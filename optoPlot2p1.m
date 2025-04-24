load('O:\sjk\DATA\imagingData\2p-opto\sk241\003_005_lightoff_lighton\suite2p\plane0\Fall.mat')

manualCell=[9,3,22,40,17,38];

% avgTrace=mean(F(manualCell,1:2000));
figure;subplot(1,2,1);hold on;title('light off cell 4');

% avgTraceNorm=avgTrace-median(avgTrace); 
% plot(smooth(avgTraceNorm,20));
normTrace=F(4,1:500)-median(F(4,800:1000));
plot(smooth(normTrace,20));ylim([-10 600]);
xline(20,'color', [.75 .75 .75]); xline(120,'color', [.75 .75 .75]);...
    xline(220,'color', [.75 .75 .75]); xline(320,'color', [.75 .75 .75]);...
    xline(420,'color', [.75 .75 .75]);
subplot(1,2,2);
avgTraceOn=F(4,2000:2500);
avgTraceOnNorm=avgTraceOn-median(avgTraceOn);
plot(smooth(avgTraceOnNorm,20)); title('light on cell 4');
xlim([0 500]);ylim([-10 600]); box off;
xline(20,'color', [.75 .75 .75]); xline(120,'color', [.75 .75 .75]);...
    xline(220,'color', [.75 .75 .75]); xline(320,'color', [.75 .75 .75]);...
    xline(420,'color', [.75 .75 .75]);
%%
load('O:\sjk\DATA\imagingData\2p-opto\sk241\003_004_lgihtoff_lighton\suite2p\plane0\Fall.mat')

% manualCell=[47,10,18,11,23,28,36,42];
figure;subplot(1,2,1);hold on;title('light off cell 4');
 cellnum=5;
normTrace=F(cellnum,1:500)-median(F(cellnum,800:1000));
plot(smooth(normTrace,20));ylim([-10 600]);
xline(20,'color', [.75 .75 .75]); xline(120,'color', [.75 .75 .75]);...
    xline(220,'color', [.75 .75 .75]); xline(320,'color', [.75 .75 .75]);...
    xline(420,'color', [.75 .75 .75]);
subplot(1,2,2);
avgTraceOn=F(cellnum,2000:2500);
avgTraceOnNorm=avgTraceOn-median(F(cellnum,800:1000));
plot(smooth(avgTraceOnNorm,20)); title('light on cell 4');
xlim([0 500]);ylim([-10 600]); box off;
xline(20,'color', [.75 .75 .75]); xline(120,'color', [.75 .75 .75]);...
    xline(220,'color', [.75 .75 .75]); xline(320,'color', [.75 .75 .75]);...
    xline(420,'color', [.75 .75 .75]);
%%
load('O:\sjk\DATA\imagingData\2p-opto\sk241\003_006_lightoff_lighton\suite2p\plane0\Fall.mat')

% manualCell=[47,10,18,11,23,28,36,42];
figure;subplot(1,2,1);hold on;title('light off cell 4');
 cellnum=2;
normTrace=F(cellnum,1:500)-median(F(cellnum,800:1000));
plot(smooth(normTrace,20));ylim([-10 600]);
xline(20,'color', [.75 .75 .75]); xline(120,'color', [.75 .75 .75]);...
    xline(220,'color', [.75 .75 .75]); xline(320,'color', [.75 .75 .75]);...
    xline(420,'color', [.75 .75 .75]);
subplot(1,2,2);
avgTraceOn=F(cellnum,2000:2500);
avgTraceOnNorm=avgTraceOn-median(F(cellnum,2000:4000));
plot(smooth(avgTraceOnNorm,20)); title('light on cell 4');
xlim([0 500]);ylim([-10 600]); box off;
xline(20,'color', [.75 .75 .75]); xline(120,'color', [.75 .75 .75]);...
    xline(220,'color', [.75 .75 .75]); xline(320,'color', [.75 .75 .75]);...
    xline(420,'color', [.75 .75 .75]);

%%

cellnum=23;cellnum=cellnum+1;
Fbatch=reshape(F(cellnum,:),[],100);
figure;subplot(1,2,1); hold on;
title(['Cell ' num2str(cellnum) ', Average Response']);
traceTemp=Fbatch(1:20,:);
% normTrace=traceTemp-median(traceTemp(1:2,:,:));
% plot(smooth(mean(normTrace(20:2000),5)));