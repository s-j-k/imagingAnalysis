pclamppath = 'O:\sjk\DATA\imagingData\2p-opto\sk272\';
fname = '2025_07_15_0028'; 
savedir =pclamppath;
[abfdata,samplingrate,header] = abfload([pclamppath,[fname, '.abf']]);
savename = [savedir,strrep(fname,'.abf','_processed')];
Fs = 1/(samplingrate*(10^-6)); % capture the sampling rate in Hz

% Identify channel with all potential channel types
optoTrig = find(strcmp(header.recChNames, 'Opto_trig'));
soundChannel = find(strcmp(header.recChNames, 'Sound'));
frameChannel = find(strcmp(header.recChNames, 'Frame'));

% Identify the samples for each frame, frame data is always channel 1
frametransitions = diff (abfdata(:,frameChannel));%abs(diff(abfdata(:,framechannel))); %calculate derivative of frame signal to identify transition
ind = find(frametransitions>0.4); %identify when there is a positive change in the signal (frame start)
ind2 = diff(ind); %the frame signal transition sometimes takes more than one sample, so need to identify the big transition points
ind3 = find(ind2>80); %length of ind3 should now be the total number of frames acquired
numsamples_perframe = round(mean(ind2(ind3))); %average number of samples per frame, e.g. 25600 for 100KhZ at 3.91 frames per second

soundtransitions= diff(abfdata(:,soundChannel));
soundInd=find(soundtransitions>0.1);
soundInd2=diff(soundInd);
soundInd3=find(soundInd2>100);

ind_end = ind(ind3); %now just identify the exact samples where frame transition occurs [[this will be "last sample in frame n"]]
ind_end = [ind_end;ind_end(end)+numsamples_perframe];
ind_start = [ind_end(1:(length(ind_end)-1))]; %create index of frame start

numframes = length(ind_end); % tot nb of frames +1

filectr = 1;
frame_recording = abfdata(:,frameChannel);
%     lickingchannel = find(strcmp(header.recChNames, 'Dig_licks'));
%     rewardchannel = find(strcmp(header.recChNames, 'Water'));
opto_recording = abfdata(:,optoTrig);
sound_recording = abfdata(:,soundChannel);
% load corresponding behavrio file and get when the stim where played +
% length stim duration


figure; ax1=subplot(2,1,1);plot(sound_recording);hold on;
ax2=subplot(2,1,2);plot(frame_recording);
linkaxes([ax1 ax2],'xy')