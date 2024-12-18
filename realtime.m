clamppath='Z:\su\DATA\imagingData\202207\sk70\20220715\';
cd(clamppath);
clampex=abfload('Z:\su\DATA\imagingData\202207\sk70\20220715\2022_07_15_0003.abf');
% 
%     {'IN 0'     }
%     {'Photogate'}
%     {'Water'    }
%     {'Sound'    }
%     {'Frame'    }
%     {'Dig_licks'}
sound=clampex(:,4);
frames=clampex(:,5);
%%
% now load the data
load('C:\Users\sjkim1\Desktop\do\sk70_001_001_realtime.mat');
figure;plot(([1:length(rtdata)]/30),smooth(rtdata(:,7),20));
%%

% find the end of the imaging session
lastFrame=7573000;
% we know the length of the session from the imaging data
ans=length(rtdata);
sampleDiff=lastFrame/ans;

xsound=((1:length(sound))/sampleDiff)/1800;
% looks like clampex is sampling at 324 samples/sec?
%%
figure;ax(1)=subplot(2,1,1);
% plot(([1:length(rtdata)]/1800),smooth(rtdata(:,7),20));
plot(([1:length(rtdata)]/1800),rtdata(:,7));
ax(2)=subplot(2,1,2);
plot(xsound,sound);

linkaxes(ax,'x');


%%

windowSize = 1; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

nrtdata=filter(b,a,rtdata);

figure;ax(1)=subplot(2,1,1);
% plot(([1:length(rtdata)]/1800),smooth(rtdata(:,7),20));
plot(([1:length(nrtdata)]/1800),nrtdata(:,6));
ax(2)=subplot(2,1,2);
plot(xsound,sound);

linkaxes(ax,'x');

%%

% now find indices of sound starting and stopping
% there are 300 sounds, each 1 sec apart, 120ms 
sStart=0.04430*1800*sampleDiff;
sEnd=0.04492*1800*sampleDiff;
