% [fname_all,pathname]=uigetfile('*.abf','Select ABF file','N:\jenni\clamp\','Multiselect','On'); % change path of directory if necessary

% Create a folder where the data is going to be saved
% mkdir(savedir);
animal = 'RoofBuddy2';
savedir = ['N:\' animal '\Microphone\']; % change path of directory if necessary

datatable = readtable(['N:\Jenni\' animal '\ProcessingProgress.xlsx'],'Format','auto');
daytraining = table2array(datatable(:,1));
framerate = table2array(datatable(:,8));
protocol = table2array(datatable(:,7));
path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,10));
pclampfile = table2array(datatable(:,11));
uniquedaytraining = unique(daytraining);
site = table2array(datatable(:,5));
uniquesite = unique(site);
% for ii = 1:length(uniquesite)
%     analysefilelistind(ii) = find(site==uniquesite(ii),1);
% end

pclamppath = ['K:\Jenni\' animal '\pclamp\' ];
% load corresponding behavior file
filenb = 144; % filenb behavior and pclamp Angie file 4 and 28, Mel 5 and
load(['K:\Jenni\' animal '\behavior\' behaviorfile{filenb}])

% % Check nb of files
% if ischar(fname_all)==1
%     numfiles = 1;
% else
%     numfiles = length(fname_all);
% end
% disp(['There will be ',num2str(numfiles),' abf files to process ...']);

% Loop through file to process
% for filectr = 1:numfiles

% Load the ABF file
% if numfiles == 1
%     fname = fname_all;
% else
%     fname = char(fname_all(filectr));
% end
% disp(['Loading ' fname]);
fname = pclampfile{filenb};
[abfdata,samplingrate,header] = abfload([pclamppath,[fname, '.abf']]);
savename = [savedir,strrep(fname,'.abf','_processed')];
Fs = 1/(samplingrate*(10^-6)); % capture the sampling rate in Hz

% Identify channel with all potential channel types
photogate = find(strcmp(header.recChNames, 'Photogate'));
sound_channel = find(strcmp(header.recChNames, 'Sound'));
framechannel = find(strcmp(header.recChNames, 'Frame'));

% Identify the samples for each frame, frame data is always channel 1
frametransitions = diff (abfdata(:,framechannel));%abs(diff(abfdata(:,framechannel))); %calculate derivative of frame signal to identify transitions
ind = find(frametransitions>0.4); %identify when there is a positive change in the signal (frame start)
ind2 = diff(ind); %the frame signal transition sometimes takes more than one sample, so need to identify the big transition points
ind3 = find(ind2>80); %length of ind3 should now be the total number of frames acquired
numsamples_perframe = round(mean(ind2(ind3))); %average number of samples per frame, e.g. 25600 for 100KhZ at 3.91 frames per second

ind_end = ind(ind3); %now just identify the exact samples where frame transition occurs [[this will be "last sample in frame n"]]
ind_end = [ind_end;ind_end(end)+numsamples_perframe];
ind_start = [ind_end(1:(length(ind_end)-1))]; %create index of frame start

numframes = length(ind_end); % tot nb of frames +1

filectr = 1;
frame_recording = abfdata(:,framechannel);
%     lickingchannel = find(strcmp(header.recChNames, 'Dig_licks'));
%     rewardchannel = find(strcmp(header.recChNames, 'Water'));
mic_recording{filectr} = abfdata(:,photogate);
sound_recording{filectr} = abfdata(:,sound_channel);
% load corresponding behavrio file and get when the stim where played +
% length stim duration

stimframeind = batimagingdata.frametone;
lengthstimincalamp = (unique(batimagingdata.durationseq)/1000+0.2)*Fs;% add 200 ms of buffer 100ms before and 100 ms after
correspondingframe = ind_start(stimframeind);

if strfind(protocol{filenb},'Mel')==1
    stimwavname = dir(['N:\Jenni\wavfiles\Mel\' '*.wav']);
    for si = 1:length({stimwavname.name})
        wav_tmp{si} = audioread(['N:\Jenni\wavfiles\Mel\' stimwavname(si).name]);%MelStim' num2str(si) '.wav']);
    end
    correctedstimseq = repmat([2:10,1],1,10);

elseif strfind(protocol{filenb},'Angie')==1
    stimwavname = dir(['N:\Jenni\wavfiles\Angie\' '*.wav']);
    for si = 1:length({stimwavname.name})
        wav_tmp{si} = audioread(['N:\Jenni\wavfiles\Angie\' stimwavname(si).name]);
    end
    correctedstimseq = repmat([2:17,1],1,10);

elseif strfind(protocol{filenb},'Pulse')==1
    stimwavname = dir(['N:\Jenni\wavfiles\PE\' '*.wav']);
    for si = 1:length({stimwavname.name})
        wav_tmp{si} = audioread(['N:\Jenni\wavfiles\PE\' stimwavname(si).name]);
    end
    correctedstimseq = repmat([2:10,1],1,10);
elseif  strfind(protocol{filenb},'Tonotopy')==1
end

clear Matmicstim Matsoundstim
for si = 1:length(stimframeind)
    if strfind(protocol{filenb},'PE')==1
        Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+0.1*Fs));
        Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+0.1*Fs));
    elseif strfind(protocol{filenb},'tonotopy')==1
        Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+0.05*Fs));
        Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+0.05*Fs));
    else
        Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):((correspondingframe(si))+(lengthstimincalamp)));
        Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+(lengthstimincalamp)));
    end
end

if strfind(protocol{filenb},'tonotopy')==1
    h=figure(); hold on;
    tonid = unique(batimagingdata.toneseq);
         count = 0;
    for tid = [1:2:10]
        clear toneindunder10 toneind10 toneind20
        toneindunder10 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq<10));
        toneind10 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq==10));
        toneind20 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq==20));
        subplot(5,3,1+count);spectrogram(mean(Matsoundstim(1:(0.050*Fs),toneindunder10),2),400,300,[],Fs,'yaxis'); ylim([0 100]);
        subplot(5,3,2+count);spectrogram(mean(Matsoundstim(1:(0.050*Fs),toneind10),2),400,300,[],Fs,'yaxis'); ylim([0 100]);
        subplot(5,3,3+count);spectrogram(mean(Matsoundstim(1:(0.050*Fs),toneind20),2),400,300,[],Fs,'yaxis'); ylim([0 100]);
        count = count+3;
    end
            saveas(h,['N:\Jenni\' animal '\figure\examplestim_' behaviorfile{filenb} '.pdf'])
elseif strfind(protocol{filenb},'Bat linear sweep and white noise')==1
    h=figure(); hold on;
    tonid = unique(batimagingdata.toneseq);
    count = 0;
    for tid = [1:3]
        clear toneindunder10 toneind10 toneind20
        toneindunder10 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq<10));
        toneind10 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq==10));
        toneind20 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq==20));
        toneind50 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq==50));
        toneind100 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq==100));
        toneind200 = intersect(find(batimagingdata.toneseq==tonid(tid)),find(batimagingdata.durationseq==200));

        subplot(3,3,1+count);spectrogram(mean(Matmicstim(1:(0.050*250000),toneindunder10),2),400,300,[],250000,'yaxis'); ylim([0 100]);
        subplot(3,3,2+count);spectrogram(mean(Matmicstim(1:(0.050*500000),toneind100),2),400,300,[],250000,'yaxis'); ylim([0 100]);
        subplot(3,3,3+count);spectrogram(mean(Matmicstim(1:(0.050*500000),toneind200),2),400,300,[],250000,'yaxis'); ylim([0 100]);
%         subplot(3,3,4+count);spectrogram(mean(Matsoundstim(1:(0.070*250000),toneind20),2),400,300,[],250000,'yaxis'); ylim([0 100]);
%         subplot(3,3,5+count);spectrogram(mean(Matsoundstim(1:(0.100*250000),toneind20),2),400,300,[],250000,'yaxis'); ylim([0 100]);
%         subplot(3,3,6+count);spectrogram(mean(Matsoundstim(1:(0.200*250000),toneind20),2),400,300,[],250000,'yaxis'); ylim([0 100]);

        count = count+3;
    end
            saveas(h,['N:\Jenni\' animal '\figure\examplestim_' behaviorfile{filenb} '.pdf'])
else
    h =figure(1); hold on;
    for si = 6:length({stimwavname.name})
        if strfind(protocol{filenb},'Pulse')==1
%             subplot(1,2,1);spectrogram(wav_tmp{si},40,30,[],200000,'yaxis'); % plot theorithecal wav file
            subplot(1,5,si-5);spectrogram(mean(Matmicstim(1:0.05*Fs,find(correctedstimseq==si)),2),40,30,[],250000,'yaxis'); ylim([0 100]); %(0.1*250000):(0.15*250000),si)),400,300,[],250000,'yaxis'); ylim([0 100]); %ylim([0,100000]); % plot sound input
       title(batimagingdata.wavnames{si})
        else
            subplot(1,2,1);spectrogram(wav_tmp{si},400,300,[],200000,'yaxis'); % plot theorithecal wav file
            subplot(1,2,2);spectrogram(mean((Matmicstim(1:round((length(wav_tmp{si})/200000)*250000),find(correctedstimseq==si))),2),400,300,[],250000,'yaxis'); ylim([0 100]); %(0.1*250000):(0.15*250000),si)),400,300,[],250000,'yaxis'); ylim([0 100]); %ylim([0,100000]); % plot sound input
        end
        
        %     subplot(1,3,3);spectrogram((Matmicstim(1:(1.4*250000),si)),400,300,[],250000,'yaxis');
        %     % plot mic recording
        %      figure(1); hold on; subplot(2,5,si);spectrogram(mean(Matmicstim(:,[(1:10:100)]+(si-1)),2),400,300,[],250000,'yaxis');
    end
            print(['N:\Jenni\' animal '\figure\examplestim_' behaviorfile{filenb} '_2'] ,'-dpdf','-fillpage')
end
% end
disp('done')