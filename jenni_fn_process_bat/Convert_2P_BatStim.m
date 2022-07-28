function Convert_2P_BatStim(owner)

switch owner
    case 'Angie'
        cd C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\Angie_stim_IC
        load('Stim2PechoFMB.mat')
        filename = fields([Stim2P]);
        nbfiles = length(fields([Stim2P]));;
        SF = 250000;
        p = 4;
        q = 5;
        for fn = 1:nbfiles
            clear bla
            %     if mod(length(bli.sound),5)>0
            downsamplesound{fn} = resample(Stim2P.(filename{fn}),p,q);
            %     elseif fn == nbfiles
            %         downsamplesound{fn} = resample(bli.sound(1:length(bli.sound)-34),p,q);
            %     else
            %         downsamplesound{fn} = resample(bli.sound,p,q);
            %     end
            figure();spectrogram(downsamplesound{fn},400,300,[],200000,'yaxis'); title(filename{fn})
            audiowrite([filename{fn} '.wav'],(int16(downsamplesound{fn}*100)),SF*(p/q),'BitsPerSample',16)
            bla = audioread([filename{fn} '.wav']);
            h = figure(); hold on; spectrogram(bla,400,300,[],200000,'yaxis'); title(filename{fn})
            axis tight; %saveas(h,'C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\wavhypPEspectrum2ms.pdf')
%             
        end
    case 'Bernhard'
        cd C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\Benhard_pup_cals
        load('Vocalizations.mat')
        SF = 250000;
        p = 4;
        q = 5;
        [trial,ind]= find(vertcat(Vocs.SNR)>5);
        filename = trial;
        nbfiles = length(trial);
        for fn = 1:nbfiles
            clear bla 
            %     if mod(length(bli.sound),5)>0
            downsamplesound{fn} = resample(Vocs(trial(fn)).Sound{ind(fn)}(1:end-1),p,q);
            %     elseif fn == nbfiles
            %         downsamplesound{fn} = resample(bli.sound(1:length(bli.sound)-34),p,q);
            %     else
            %         downsamplesound{fn} = resample(bli.sound,p,q);
            %     end
            figure();spectrogram(downsamplesound{fn},400,300,[],200000,'yaxis'); title(['voc trial ' trial(fn) ' ind ' ind(fn)])
%             audiowrite(['voc trial ' num2str(trial(fn)) ' ind ' num2str(ind(fn)) '.wav'],int16(downsamplesound{fn}*100),SF*(p/q),'BitsPerSample',16)
%             bla = audioread(['voc trial ' num2str(trial(fn)) ' ind ' num2str(ind(fn)) '.wav']);
%             h = figure(); hold on; spectrogram(bla,400,300,[],200000,'yaxis'); 
%             axis tight; %saveas(h,'C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\wavhypPEspectrum2ms.pdf')
        end
    case 'Mel'
        cd C:\Users\jlawlor3\Dropbox\aav_histo/2P_stimuli/Mel_stim_SC/ %cd /home/jennifer/Dropbox/aav_histo/2P_stimuli/Mel_stim_SC/
        load('stimuli_for_2P_bats.mat')
        nbfiles = length(auditory_stimuli_zero_pad);
        SF = 500000;
        p = 2;
        q = 5;
        for fn = 1:nbfiles
            %     if mod(length(bli.sound),5)>0
            downsamplesound{fn} = resample(auditory_stimuli_zero_pad{fn},p,q);
            figure();spectrogram(downsamplesound{fn},400,300,[],200000,'yaxis'); %title(['voc trial ' trial(fn) ' ind ' ind(fn)])
%             if fn == 5
%                 scrambled = [downsamplesound{fn}(64001:end),downsamplesound{fn}(1:32000),downsamplesound{fn}(32001:64000)];
%                 audiowrite(['MelStim' num2str(fn) 'scrambled.wav'],downsamplesound{fn},SF*(p/q),'BitsPerSample',16)
%             end
            audiowrite(['MelStim' num2str(fn) '.wav'],int16(downsamplesound{fn}),SF*(p/q),'BitsPerSample',16)
            bla = audioread(['MelStim' num2str(fn) '.wav']);
            h = figure(); hold on; spectrogram(bla,400,300,[],200000,'yaxis'); title(['MelStim' num2str(fn) '.wav'])
        end
        
    case 'Sam'
        cd C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\SamStims\final_stims165
        p = 9;
        q = 8;
        [xt,fs]=audioread('stim368_thunder.wav');
        downsamplesound = resample(xt,p,q);
        %     elseif fn == nbfiles
        %         downsamplesound{fn} = resample(bli.sound(1:length(bli.sound)-34),p,q);
        %     else
        %         downsamplesound{fn} = resample(bli.sound,p,q);
        %     end
        figure();spectrogram(downsamplesound,400,300,[],50000,'yaxis'); title('stim368_thunder.wav')
        audiowrite(['stim368_thunder_resample.wav'],downsamplesound,50000,'BitsPerSample',16)
        
    case 'HypSweep'
        %     alpha	= 15*n*pi/1024;
        %     beta    = 5*n*pi/1024;
        %     t  	= (1.001:1:n+.001)./n;
        %     f1      = zeros(1,n);
        %     f2      = zeros(1,n);
        %     f1  	= sin(alpha./(.8-t)).*(0.1<t).*(t<0.68);
        %     f2  l= sin(beta./(.8-t)).*(0.1<t).*(t<0.75);
        %     M  	= round(0.65*n);
        %     P 	= floor(M/4);
        %     enveloppe = ones(1,M); % the rising cutoff function
        %     enveloppe(1:P) = (1+sin(-pi/2+((1:P)-ones(1,P))./(P-1)*pi))/2;
        %     enveloppe(M-P+1:M) = reverse(enveloppe(1:P));
        %     env 	= zeros(1,n);
        %     env(ceil(n/10):M+ceil(n/10)-1) = enveloppe(1:M);
        %     sig     = (f1+f2).*env;
        %
        %
        %     t = (0:1/1e3:2)/1000;
        %     f0 = 450;
        %     f1 = 10000;
        %     T = 1e-3;
        %
        %     w = cosh([flip(chirp(t,f0,T,f1,'linear'))]);
        %     pspectrum(w',1e6,'spectrogram','TimeResolution',0.0001,'OverlapPercent',99,'Leakage',0.85);
        %
        %     Fs = 100000;
        %     t = 0:1/Fs:0.002;
        %     f = 1000;
        %     fi = [];
        % %     temp =0
        % %     for i = 1:length(t)
        % %         fi = [fi temp];
        % %         temp = temp + (0.000025*f)
        % %     end
        %     fi = cosh(linspace(0,1.32,201));
        %     y = sin(2*pi*t.*fi);
        %     plot(flip(y))
        %     figure(); hold on;
        %     subplot(1,2,1); plot(fi)
        %     subplot(1,2,2);plot(flip(y));
        %     pspectrum(y',2e5,'spectrogram','TimeResolution',0.0001,'OverlapPercent',99,'Leakage',0.85);
        
        
        %     sound(flip(y),Fs)
        
        %     Fs=200000; % sample rate
        %     tf=0.002; % 2 seconds
        %     t=0:1/Fs:tf-1/Fs;
        %     f1=25000;
        %     f2=90000; % start @ 100 Hz, go up to 400Hz
        %     semi_t=0:1/Fs:(tf/2-1/Fs);
        %     sl=2*(f2-f1/2);
        %     f1=f1*semi_t+(sl.*semi_t/2);
        %     f2=f1(end)+f2*semi_t-sl.*semi_t/2;
        %     f=[f1 f2];
        %     y=1.33*cos(2*pi*f.*t);
        
        hypchirp = load('C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\hypchirptest.txt'); % sample rate for creation 2000000
        %     h = figure(1); hold on; subplot(1,2,1); pspectrum(logchirp,2e6,'spectrogram','TimeResolution',0.0005,'OverlapPercent',99,'Leakage',0.85);
        hypchirpresample = resample(hypchirp(1:4000),1,10); % extra 0.1ms to avoid aliasing at the end of the stim
        % add 0.1 ramp
        linrampup = 0+1/20:1/20:1;
        linrampdown = flip(linrampup);
        hypchirpresample(1:20) = (hypchirpresample(1:20)).*linrampup';
        hypchirpresample(end-19:end) = (hypchirpresample(end-19:end)).*linrampdown';
        
        %     subplot(1,2,2); hold on;
        %     pspectrum(logchirpresample,2e5,'spectrogram','TimeResolution',0.0005,'OverlapPercent',99,'Leakage',0.85);
        %     axis tight
        %     saveas(h,'K:\Jenni\logchirptest.pdf')
            h = figure(1); hold on;
        for ii = 1:10
            PE{ii} = [hypchirpresample; zeros(1,400*ii)';(hypchirpresample./(20*log10(7/5)))];% factore 10 between 70 dB SPL and 50 in Pa %;%(5/7))];% zeros(1,200)';; zeros(1,200)'
            reversePE{ii} = flipud(PE{ii}); 
            %                     subplot(2,5,ii); pspectrum(PE{ii},2e5,'spectrogram','TimeResolution',0.0005,'OverlapPercent',99,'Leakage',0.85);
            audiowrite(['C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\PE' num2str(2*ii) 'ms.wav'],int16(PE{ii}),2e5);
            audiowrite(['C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\reversePE' num2str(2*ii) 'ms.wav'],int16(reversePE{ii}),2e5);
        end
        simple_PE = [hypchirpresample];
        simple_Echo = (hypchirpresample./(20*log10(7/5)));
        audiowrite(['C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\PE_simple.wav'],int16(simple_PE),2e5);
        audiowrite(['C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\PE_simpleecho.wav'],int16(simple_PE),2e5);

        %     saveas(h,'C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\hypPEspectrum.pdf')
        bla = audioread(['C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\reversePE2ms.wav']);
        h = figure(); hold on; pspectrum(bla,2e5,'spectrogram','TimeResolution',0.0005,'OverlapPercent',99,'Leakage',0.85);
        axis tight; saveas(h,'C:\Users\jlawlor3\Dropbox\aav_histo\2P_stimuli\PEpairs\wavhypPEspectrum2ms.pdf')
        
        
        
end
end