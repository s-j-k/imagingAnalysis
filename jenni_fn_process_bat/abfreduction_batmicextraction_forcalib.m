% [fname_all,pathname]=uigetfile('*.abf','Select ABF file','N:\jenni\clamp\','Multiselect','On'); % change path of directory if necessary

% Create a folder where the data is going to be saved
savedir = 'N:\RoofBuddy2\Microphone\'; % change path of directory if necessary
% mkdir(savedir);
animal = 'RoofBuddy2';
datatable = readtable('C:\Users\kklab\Desktop\TPM_code_software_2018\TDT_NLW\BatSim\CalibfileRight.xlsx');
% daytraining = table2array(datatable(:,1));
framerate = table2array(datatable(:,4));
protocol = table2array(datatable(:,3));
% path = table2array(datatable(:,9));
behaviorfile = table2array(datatable(:,6));
pclampfile = table2array(datatable(:,7));
% uniquedaytraining = unique(daytraining);
% site = table2array(datatable(:,5));
% uniquesite = unique(site);
% for ii = 1:length(uniquesite)
%     analysefilelistind(ii) = find(site==uniquesite(ii),1);
% end

pclamppath = ['D:\pClamp\Data\' ];
% load corresponding behavior file
filenb = 10;%6; % filenb behavior and pclamp Angie file 4 and 28, Mel 5 and
load(['D:\batdata\' behaviorfile{filenb}])

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
frametransitions = diff(abfdata(:,framechannel));%abs(diff(abfdata(:,framechannel))); %calculate derivative of frame signal to identify transitions
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
lengthstimincalamp = (unique(batimagingdata.durationseq)/1000+0.2)*Fs;% add 200 ms of buffe 100ms before and 100 ms after
correspondingframe = ind_start(stimframeind);

if strfind(protocol{filenb},'Mel')==1
    %     stimwavname = dir(['N:\Jenni\wavfiles\Mel\' '*.wav']);
    %     for si = 1:length({stimwavname.name})
    %         wav_tmp{si} = audioread(['N:\Jenni\wavfiles\Mel\' stimwavname(si).name]);%MelStim' num2str(si) '.wav']);
    %     end
    %     correctedstimseq = repmat([2:10,1],1,10);
    tonelist = unique(batimagingdata.toneseq);
    for si = 1:length(tonelist)
        ordertone(si) = find(batimagingdata.toneseq==tonelist(si),1);
    end
    
elseif strfind(protocol{filenb},'Angie')==1
    %     stimwavname = dir(['N:\Jenni\wavfiles\Angie\' '*.wav']);
    %     for si = 1:length({stimwavname.name})
    %         wav_tmp{si} = audioread(['N:\Jenni\wavfiles\Angie\' stimwavname(si).name]);
    %     end
    %     correctedstimseq = repmat([2:17,1],1,10);
    tonelist = unique(batimagingdata.toneseq);
    for si = 1:length(tonelist)
        ordertone(si) = find(batimagingdata.toneseq==tonelist(si),1);
    end
    
elseif strfind(protocol{filenb},'Pulse')==1
    %     stimwavname = dir(['N:\Jenni\wavfiles\PE\' '*.wav']);
    %     for si = 1:length({stimwavname.name})
    %         wav_tmp{si} = audioread(['N:\Jenni\wavfiles\PE\' stimwavname(si).name]);
    %     end
    tonelist = unique(batimagingdata.toneseq);
    for si = 1:length(tonelist)
        ordertone(si) = find(batimagingdata.toneseq==tonelist(si),1);
    end
elseif strfind(protocol{filenb},'Tonotopy')==1
    tonelist = unique(batimagingdata.toneseq);
    for si = 1:length(tonelist)
        ordertone(si) = find(batimagingdata.toneseq==tonelist(si),1);
    end
elseif strfind(protocol{filenb},'Linear')==1
    tonelist = unique(batimagingdata.toneseq);
    for si = 1:length(tonelist)
        ordertone(si) = find(batimagingdata.toneseq==tonelist(si),1);
    end
end

clear Matmicstim Matsoundstim
for si = 1:length(stimframeind)
    if strfind(protocol{filenb},'Pulse')==1
        Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+0.05*Fs));
        Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+0.05*Fs));
    elseif strfind(protocol{filenb},'Mel')==1
        Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):((correspondingframe(si))+(lengthstimincalamp)));
        Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+(lengthstimincalamp)));
    elseif  strfind(protocol{filenb},'Angie')==1
        Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):((correspondingframe(si))+(lengthstimincalamp)));
        Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+(lengthstimincalamp)));
    elseif strfind(protocol{filenb},'Tonotopy')==1
        Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):((correspondingframe(si))+max(lengthstimincalamp)));
        Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+max(lengthstimincalamp)));
    elseif strfind(protocol{filenb},'Linear')==1
        Matmicstim(:,si) = mic_recording{filectr}((correspondingframe(si)):((correspondingframe(si))+max(lengthstimincalamp)));
        Matsoundstim(:,si) = sound_recording{filectr}((correspondingframe(si)):(correspondingframe(si)+max(lengthstimincalamp)));
    end
end

if strfind(protocol{filenb},'Tonotopy')==1
    h =figure(1); hold on;
    for si = 1:length(tonelist)
        subplot(2,5,si);spectrogram((Matmicstim(1:55001,ordertone(si))),400,300,[],250000,'yaxis'); ylim([0 100]);
        maxdb(si)=max(20*log10(abs(Matmicstim(1:55001,ordertone(si))-0.4)/20e-5));% the channel for the mic is now ins mPa, there's a 20dB gain, so 0.0004 pa, so minus 0.4%20*log10((abs(mean(Matmicstim(1:55001,si:10:100),2))*(0.9/20))/20e-5)
    end
    saveas(h,['C:\Users\kklab\Desktop\TPM_code_software_2018\TDT_NLW\BatSim\Tones' num2str(date) '.pdf'])
elseif strfind(protocol{filenb},'Linear')==1
    h =figure(1); hold on;
    for si = 1:length(tonelist)
        subplot(2,5,si);spectrogram((Matmicstim(1:55001,ordertone(si))),400,300,[],250000,'yaxis'); ylim([0 100]);
        maxdb(si)=max(20*log10(abs(Matmicstim(1:55001,ordertone(si))-0.4)/20e-5));% the channel for the mic is now ins mPa, there's a 20dB gain, so 0.0004 pa, so minus 0.4%20*log10((abs(mean(Matmicstim(1:55001,si:10:100),2))*(0.9/20))/20e-5)
    end
    saveas(h,['C:\Users\kklab\Desktop\TPM_code_software_2018\TDT_NLW\BatSim\Linear' num2str(date) '.pdf'])
elseif strfind(protocol{filenb},'Pulse')==1
    h =figure(1); hold on;
    for si = 1:length(tonelist)
        subplot(2,5,si);spectrogram((Matmicstim(1:12501,ordertone(si))),400,300,[],250000,'yaxis'); ylim([0 100]);
        maxdb(si)=max(20*log10(abs(Matmicstim(1:12501,ordertone(si))-0.4)/20e-5));% the channel for the mic is now ins mPa, there's a 20dB gain, so 0.0004 pa, so minus 0.4%20*log10((abs(mean(Matmicstim(1:55001,si:10:100),2))*(0.9/20))/20e-5)
    end
    saveas(h,['C:\Users\kklab\Desktop\TPM_code_software_2018\TDT_NLW\BatSim\Pulse' num2str(date) '.pdf'])
elseif strfind(protocol{filenb},'Angie')==1
    h =figure(1); hold on;
    for si = 1:length(tonelist)
        subplot(2,10,si);spectrogram((Matmicstim(1:61751,ordertone(si))),400,300,[],250000,'yaxis'); ylim([0 100]);
        maxdb(si)=max(20*log10(abs(Matmicstim(1:61751,ordertone(si))-0.4)/20e-5));% the channel for the mic is now ins mPa, there's a 20dB gain, so 0.0004 pa, so minus 0.4%20*log10((abs(mean(Matmicstim(1:55001,si:10:100),2))*(0.9/20))/20e-5)
    end
    saveas(h,['C:\Users\kklab\Desktop\TPM_code_software_2018\TDT_NLW\BatSim\Angie' num2str(date) '.pdf'])
    elseif strfind(protocol{filenb},'Mel')==1
    h =figure(1); hold on;
    for si = 1:length(tonelist)
        subplot(2,10,si);spectrogram((Matmicstim(1:400000,ordertone(si))),400,300,[],250000,'yaxis'); ylim([0 100]);
        maxdb(si)=max(20*log10(abs(Matmicstim(1:400000,ordertone(si))-0.4)/20e-5));% the channel for the mic is now ins mPa, there's a 20dB gain, so 0.0004 pa, so minus 0.4%20*log10((abs(mean(Matmicstim(1:55001,si:10:100),2))*(0.9/20))/20e-5)
    end
    saveas(h,['C:\Users\kklab\Desktop\TPM_code_software_2018\TDT_NLW\BatSim\Mel' num2str(date) '.pdf'])
    
end


% for si = 1:length({stimwavname.name})
%    clear h
%     h =figure(si); hold on;
%     if strfind(protocol{filenb},'Pulse')==1
%         subplot(1,2,1);spectrogram(wav_tmp{si},40,30,[],200000,'yaxis'); % plot theorithecal wav file
%         subplot(1,2,2);spectrogram(mean(Matsoundstim(1:round((length(wav_tmp{si})/200000)*250000)+500,find(correctedstimseq==si)),2),40,30,[],250000,'yaxis'); ylim([0 100]); %(0.1*250000):(0.15*250000),si)),400,300,[],250000,'yaxis'); ylim([0 100]); %ylim([0,100000]); % plot sound input
%     else
%         subplot(1,2,1);spectrogram(wav_tmp{si},400,300,[],200000,'yaxis'); % plot theorithecal wav file
%         subplot(1,2,2);spectrogram(mean((Matsoundstim(1:round((length(wav_tmp{si})/200000)*250000),find(correctedstimseq==si))),2),400,300,[],250000,'yaxis'); ylim([0 100]); %(0.1*250000):(0.15*250000),si)),400,300,[],250000,'yaxis'); ylim([0 100]); %ylim([0,100000]); % plot sound input
%     end
%         saveas(h,['N:\Jenni\' animal '\figure\examplestim_' stimwavname(si).name '.pdf'])
%
%     %     subplot(1,3,3);spectrogram((Matmicstim(1:(1.4*250000),si)),400,300,[],250000,'yaxis');
% %     % plot mic recording
%     %      figure(1); hold on; subplot(2,5,si);spectrogram(mean(Matmicstim(:,[(1:10:100)]+(si-1)),2),400,300,[],250000,'yaxis');
% end
% end
disp('done')