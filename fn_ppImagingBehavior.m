function fn_ppImagingBehavior(animal)
data_folder='D:\';
animal_list ={animal};
plotcellsindividually=1;
% plotcellsprepostprobe=0;
distancetoplaqueanalysis=0;

for fl=length(animal_list)

    Bfolder='New folder\';
    behaviour_folder=[data_folder cell2mat(animal_list(fl))  filesep Bfolder filesep 'behavior\']; 
    suite2ppathB=[data_folder cell2mat(animal_list(fl))  filesep Bfolder filesep 'suite2p\']; 
    distancetoplaquepath =[data_folder cell2mat(animal_list(fl))  filesep Bfolder filesep ]; 

    cd(behaviour_folder);
    IBfiles = dir('*.mat');
    name1IBf = IBfiles(1).name(1:end-3);
    infosbxa = load([name1IBf 'mat']); 
    infosbx = infosbxa.info;
    
    if isempty(infosbx.otparam)
        nPlanes = 1;
    else
        nPlanes = infosbx.otparam(3);
    end
    
    if infosbx.scanmode == 0
        factor = 2; 
    else
        factor = 1;
    end
    framerate = factor*(infosbx.resfreq/infosbx.config.lines);
    
 
    for i=1:nPlanes
        cd([suite2ppathB 'plane' num2str(i-1)]);
        TCBdata = load('Fall.mat'); 
        TCB(i).plane=TCBdata;
        if distancetoplaqueanalysis==1
            cd(distancetoplaquepath); 
            distancetoplaquedata = load('neurondist.mat'); 
            NeuronDist(i).plane=distancetoplaquedata.neurondist{i};
        end
    end
    
    if nPlanes>2
        frameswindowtouse=40;
        frameswindowtouse1=20;
        framesbaseline=8;
        framesbaseline1=1;
    else
        frameswindowtouse=80;
        frameswindowtouse1=40;
        framesbaseline=16;
        framesbaseline1=1;
    end
    
    BehaviorImagingData(fl).mouse=cell2mat(animal_list(fl));
    BehaviorImagingData(fl).behaviorTrace=TCB;
    
    % Column names
    SESSION = 1;
    TRIAL = 2;
    TONE = 3;
    RESPONSE = 4;
    H = 1; M = 2; FA = 3; CR = 4;
    NOLICK_PERIOD = 5;
    RESP_TIME = 6;
    DELAY_AFTER_RESP = 7;
    TOTAL_TRIAL_DUR_MINUS_RESP_TIME = 8;
    LICKFrame = 9;%frame when they lick
    TOTAL_TRIAL_DUR = 10; % toc
    TONEF = 11;
    FRAME = 12;
    CONTEXT=13;
  
    framerateplane=framerate/nPlanes;
    ToneDuration = 0.1; % 100 ms

    
    cd(behaviour_folder);
    Bfiles = struct2table(dir('*.txt'));
    nBfiles = size(Bfiles,1);

    Bfilesnames = cellfun(@(x)(strsplit(x,'.')),Bfiles.name,'UniformOutput',false);
    row = cellfun(@(x) size(x,2),Bfilesnames);
    iwant = Bfilesnames(row==2);
    firstframe= 1;
    lastframe= 0;
    for bfilnb=1:nBfiles
        behaviordatatmp = ([behaviour_folder filesep cell2mat(iwant{bfilnb})]);
        Bdataraw=importdata([behaviordatatmp(1:end-3) '.txt']);
        Bdata(bfilnb).session=Bdataraw;
        NofTrials = size(Bdataraw,1);
        framespersession=BehaviorImagingData(fl).behaviorTrace(1).plane.ops.nframes_per_folder(bfilnb);
        nchannels = BehaviorImagingData(fl).behaviorTrace(1).plane.ops.nchannels;
        framespersessionperplane=framespersession/nchannels;
        
        ToneId{bfilnb} = Bdataraw(:,TONE);
        Hit{bfilnb} =find(Bdataraw(:,RESPONSE)== H); 
        Falarm{bfilnb} = find(Bdataraw(:,RESPONSE)== FA); 
        Crejection{bfilnb} = find(Bdataraw(:,RESPONSE)== CR); 
        Miss{bfilnb} = find(Bdataraw(:,RESPONSE) == M);
        Contextvect{bfilnb} = Bdataraw(:,CONTEXT);
        Targetone{bfilnb} = Bdataraw(1,TONE);
        pairoftones{bfilnb} = unique(Bdataraw(:,TONE),'stable');
        Foiltonepos{bfilnb} = find(pairoftones{bfilnb}~=Targetone{bfilnb});
        Foiltone{bfilnb} =pairoftones{bfilnb}(Foiltonepos{bfilnb});
        Targetind{bfilnb} = find(ToneId{bfilnb}==Targetone{bfilnb});
        Foilind{bfilnb} = find(ToneId{bfilnb}~=Targetone{bfilnb});
        Probetrialsind{bfilnb}= find(Bdataraw(:,CONTEXT)== 0);
        Reinforcedtrialsind{bfilnb}= find(Bdataraw(:,CONTEXT)== 1);
        
        Probeind{bfilnb} = find(Contextvect{bfilnb}==0); 
  
        lastframe = framespersessionperplane+lastframe;

        Framenbpre = floor(Bdataraw(:,FRAME)/(nPlanes));
        Framenb{bfilnb}= Framenbpre -framesbaseline;
               
        WindowActivity = frameswindowtouse; 
        WindowActivity2 = frameswindowtouse1;
        WindowActivity1= 9;
        lickframenb{bfilnb}= floor(Bdataraw(:,LICKFrame)/nPlanes)-framesbaseline;
 
        
        ToneEvokedFrames{bfilnb} = cellfun(@(y) (y:y+WindowActivity),mat2cell(Framenb{bfilnb},ones(1,length(Framenb{bfilnb})),1),'UniformOutput',false);
        trialnb{bfilnb}=size(ToneEvokedFrames{bfilnb},1);
        
        
        for nbplane=1:nPlanes

%            smoothdatatouse{nbplane}= smoothdata(BehaviorImagingData(fl).behaviorTrace(nbplane).plane.F,2,'gaussian',4);
            smoothdatatouse{nbplane}= BehaviorImagingData(fl).behaviorTrace(nbplane).plane.F;
            traceall{nbplane}= smoothdatatouse{nbplane};%./mean(smoothdatatouse{nbplane},2);%
            cellsroiind{nbplane} = find(BehaviorImagingData(fl).behaviorTrace(nbplane).plane.iscell(:,1)==1);
            trace{bfilnb,nbplane} = traceall{nbplane}(cellsroiind{nbplane},firstframe:lastframe);%only iscells
            nbcellsperplane{nbplane} =size(trace{bfilnb,nbplane},1);
  
            ToneEvokedMat{bfilnb,nbplane} = zeros(trialnb{bfilnb},WindowActivity+1,nbcellsperplane{nbplane});
            RewardEvokedMat{bfilnb,nbplane} = zeros(trialnb{bfilnb},WindowActivity2+1,nbcellsperplane{nbplane});
            ToneEvokedMatraw{bfilnb,nbplane} = zeros(trialnb{bfilnb},WindowActivity+1,nbcellsperplane{nbplane});
            PreToneMat{bfilnb,nbplane} = zeros(trialnb{bfilnb},WindowActivity1+1,nbcellsperplane{nbplane});

            for cell=1:nbcellsperplane{nbplane}
                tracetemp{cell}=trace{bfilnb,nbplane}(cell,:);
                for trial=1:trialnb{bfilnb}
                    clear evockedtrace meanfirstframesET fstart fend
                    fstart=Framenb{bfilnb}(trial,:); 
                    fstartreward=lickframenb{bfilnb}(trial,:); 
                    fend=fstart+WindowActivity; 
                    if fstartreward>size(tracetemp{cell},2)
                        fstartreward=fstart;
                    end 
                    fendreward=fstartreward+WindowActivity2;
                    evockedtrace=tracetemp{cell}(1,fstart:fend);
                    meanfirstframesET=mean(evockedtrace(framesbaseline1:framesbaseline));
                    ToneEvokedMat{bfilnb,nbplane}(trial,:,cell) = (evockedtrace-meanfirstframesET)/meanfirstframesET;
                    ToneEvokedMatraw{bfilnb,nbplane}(trial,:,cell) = evockedtrace;
                    evockedtracereward=tracetemp{cell}(1,fstartreward:fendreward);
                    RewardEvokedMat{bfilnb,nbplane}(trial,:,cell) = evockedtracereward;
                end

            end

            lastvaluewindow=frameswindowtouse-framesbaseline;
            
            WindowTest = framesbaseline:lastvaluewindow; 
            WindowTest1 = 1:framesbaseline;% 1:10

%             for iii = 1:size(ToneEvokedMat{bfilnb,nbplane},3)
%                 testresults{bfilnb,nbplane}(iii,:) = ttest(mean(ToneEvokedMat{bfilnb,nbplane}(:,WindowTest,iii),2),mean(ToneEvokedMat{bfilnb,nbplane}(:,WindowTest1,iii),2));
%             end
                
        end
        firstframe= framespersessionperplane+firstframe;
        
    end

       
    % Pull all files together per session
    cumsumtrials = cumsum([0 cellfun(@(x) size(x,1),ToneId)]);
    BlockId = cell2mat(cellfun(@(x,y) ones(1,size(x,1))*y,ToneId,mat2cell([1:bfilnb],1,ones(1,bfilnb)),'UniformOutput',false));
    cumsumtrials(end) = [];
    TargetIndConca = cell2mat(cellfun(@(x,y) (x+y)',Targetind,mat2cell(cumsumtrials,1,ones(1,length(cumsumtrials))),'UniformOutput',0));
    FoilindConca = cell2mat(cellfun(@(x,y) (x+y)',Foilind,mat2cell(cumsumtrials,1,ones(1,length(cumsumtrials))),'UniformOutput',0));
    ProbeConca = cell2mat(cellfun(@(x,y) (x+y)',Probeind,mat2cell(cumsumtrials,1,ones(1,length(cumsumtrials))),'UniformOutput',0));

    HitConca = cell2mat(cellfun(@(x,y) (x+y)',Hit,mat2cell(cumsumtrials,1,ones(1,length(cumsumtrials))),'UniformOutput',0));
    FAConca = cell2mat(cellfun(@(x,y) (x+y)',Falarm,mat2cell(cumsumtrials,1,ones(1,length(cumsumtrials))),'UniformOutput',0));
    CRConca = cell2mat(cellfun(@(x,y) (x+y)',Crejection,mat2cell(cumsumtrials,1,ones(1,length(cumsumtrials))),'UniformOutput',0));
    MissConca = cell2mat(cellfun(@(x,y) (x+y)',Miss,mat2cell(cumsumtrials,1,ones(1,length(cumsumtrials))),'UniformOutput',0));
    
    if ~isempty(ProbeConca)
        HitProbe = HitConca(ismember(HitConca,ProbeConca));
        FAProbe = FAConca(ismember(FAConca,ProbeConca));
        MissProbe = MissConca(ismember(MissConca,ProbeConca));
        CRProbe = CRConca(ismember(CRConca,ProbeConca));
        TargetProbe = TargetIndConca(ismember(TargetIndConca,ProbeConca));
        FoilProbe = FoilindConca(ismember(FoilindConca,ProbeConca));
        
        TargetIndConca(ismember(TargetIndConca,ProbeConca)) = [];
        FoilindConca(ismember(FoilindConca,ProbeConca)) = [];
        HitConca(ismember(HitConca,ProbeConca)) = [];
        FAConca(ismember(FAConca,ProbeConca)) = [];
        CRConca(ismember(CRConca,ProbeConca)) = [];
        MissConca(ismember(MissConca,ProbeConca)) = [];
        
        probestartsat=min(ProbeConca);
        probefinishat=max(ProbeConca);
        Hitpreprobe=HitConca(HitConca < probestartsat);
        Hitpostprobe=HitConca(HitConca > probefinishat);
        FApreprobe=FAConca(FAConca < probestartsat);
        FApostprobe=FAConca(FAConca > probefinishat);
        CRpreprobe=CRConca(CRConca < probestartsat);
        CRpostprobe=CRConca(CRConca > probefinishat);
        Misspreprobe=MissConca(MissConca < probestartsat);
        Misspostprobe=MissConca(MissConca > probefinishat);
    end
    

    for nbplane = 1:nPlanes
%     clear ROINorm
        ToneEvokedMatConca{nbplane} = vertcat(ToneEvokedMat{:,nbplane});%for tone evoked responses
        ToneEvokedMatrawConca{nbplane} = vertcat(ToneEvokedMatraw{:,nbplane});%for tone evoked responses
        RewardEvokedMat{nbplane} = vertcat(RewardEvokedMat{:,nbplane});%for tone evoked responses
        LickFrameConca{nbplane} = vertcat(lickframenb{1,:});
        SelectedCells= intersect(find(testresults{1,nbplane}),find(testresults{2,nbplane}));%for tone evoked responses
        SelectedCells_save{nbplane}=SelectedCells;
%         ROINorm{nbplane}= zeros(length(BlockId),WindowActivity+1,length(SelectedCells));
        ROINorm{nbplane}= zeros(length(BlockId),WindowActivity+1,size(cellsroiind{nbplane},1));
  
        
%         define plot 
        colorsHFA=[46,196,182;255,159,28]./255; %hit and false alarm same as behaviour plots
        colorsHFAl=[203,243,240;255,235,210]./255;% same but light colors for individual plots hit and false alarm same as behaviour plots
        colorsMCR=[14,62,127;196,46,60]./255;% miss and correct rejects average
        colorsMCRl=[158,177,203;231,171,177]./255;% miss and correct rejects light for individual plots
        
        
        lastvaluewindow=frameswindowtouse-framesbaseline;
        xvalue = -framesbaseline:lastvaluewindow;% sets 0 at tone onset before ita was 53 63
        mylim=[-0.02 0.06]; % y lim for average plots
        

% %         for ii = 1:length(SelectedCells)%1:100%[1:100]%%size(ToneEvokedMat{1},3)
        for ii = 1:size(cellsroiind{nbplane},1)%1:100%[1:100]%%size(ToneEvokedMat{1},3)
            
            
            ROINorm{nbplane}(:,:,ii) = (ToneEvokedMatConca{nbplane}(:,:,ii));
            if plotcellsindividually==1
                fig1=figure; hold on; 
                subplot(2,8,1); hold on; title('Hit Reinf.'); plot(xvalue,(ROINorm{nbplane}(HitConca,:,ii).'),'Color',colorsHFAl(1,:),'LineWidth',0.5); plot(xvalue,median(ROINorm{nbplane}(HitConca,:,ii),1),'Color',colorsHFA(1,:),'LineWidth',2);
                xlabel('Frame nb'); ylabel('f/f0');          
                xlim([xvalue(1),xvalue(end)]);ylim([min(min(ROINorm{nbplane}(:,:,ii))) max(max(ROINorm{nbplane}(:,:,ii)))]);
%                 licktime=median(LickFrameConca{nbplane}(HitConca,:));
%                 xline(licktime,'-'); %,{'reward'}
%                 xline(0,'-'); %,{'tone'}

                if ~isempty(MissConca)
                    subplot(2,8,2); hold on; title('Miss Reinf.'); plot(xvalue,(ROINorm{nbplane}(MissConca,:,ii).'),'Color',colorsMCRl(1,:),'LineWidth',0.5); plot(xvalue,median(ROINorm{nbplane}(MissConca,:,ii),1),'Color',colorsMCR(1,:),'LineWidth',2);
                    xlim([xvalue(1),xvalue(end)]);ylim([min(min(ROINorm{nbplane}(:,:,ii))) max(max(ROINorm{nbplane}(:,:,ii)))]);
                    xlabel('Frame nb');
                else
                    subplot(2,8,2); hold on; title('Miss Reinf.');

                end
                if ~isempty(ProbeConca)
                    if ~isempty(HitProbe)
                        subplot(2,8,3); hold on; title('Hit Probe'); plot(xvalue,(ROINorm{nbplane}(HitProbe,:,ii).'),'Color',colorsHFAl(1,:),'LineWidth',0.5); plot(xvalue,median(ROINorm{nbplane}(HitProbe,:,ii),1),'Color',colorsHFA(1,:),'LineWidth',2);
                        xlim([xvalue(1),xvalue(end)]);ylim([min(min(ROINorm{nbplane}(:,:,ii))) max(max(ROINorm{nbplane}(:,:,ii)))]);
                        xlabel('Frame nb');
                    end
                    if ~isempty(MissProbe)

                        subplot(2,8,4); hold on; title('Miss Probe'); plot(xvalue,(ROINorm{nbplane}(MissProbe,:,ii).'),'Color',colorsMCRl(1,:),'LineWidth',0.5); plot(xvalue,median(ROINorm{nbplane}(MissProbe,:,ii),1),'Color',colorsMCR(1,:),'LineWidth',2);
                        xlim([xvalue(1),xvalue(end)]);ylim([min(min(ROINorm{nbplane}(:,:,ii))) max(max(ROINorm{nbplane}(:,:,ii)))]);
                        xlabel('Frame nb');
                    end

                end
                subplot(2,8,10); hold on; title('FA Reinf.'); plot(xvalue,(ROINorm{nbplane}(FAConca,:,ii).'),'Color',colorsHFAl(2,:),'LineWidth',0.5); plot(xvalue,median(ROINorm{nbplane}(FAConca,:,ii),1),'Color',colorsHFA(2,:),'LineWidth',2);
                xlim([xvalue(1),xvalue(end)]);ylim([min(min(ROINorm{nbplane}(:,:,ii))) max(max(ROINorm{nbplane}(:,:,ii)))]);
                xlabel('Frame nb');
                
                subplot(2,8,9); hold on; title('CR Reinf.'); plot(xvalue,(ROINorm{nbplane}(CRConca,:,ii).'),'Color',colorsMCRl(2,:),'LineWidth',0.5); plot(xvalue,median(ROINorm{nbplane}(CRConca,:,ii),1),'Color',colorsMCR(2,:),'LineWidth',2);
                ylabel('f/f0');  
                xlim([xvalue(1),xvalue(end)]);ylim([min(min(ROINorm{nbplane}(:,:,ii))) max(max(ROINorm{nbplane}(:,:,ii)))]);
                xlabel('Frame nb');

                if ~isempty(ProbeConca)
                    if isempty(FAProbe)
                        subplot(2,8,12); hold on;

                    else
                        subplot(2,8,12); hold on; title('FA Probe'); plot(xvalue,(ROINorm{nbplane}(FAProbe,:,ii).'),'Color',colorsHFAl(2,:),'LineWidth',0.5); plot(xvalue,median(ROINorm{nbplane}(FAProbe,:,ii),1),'Color',colorsHFA(2,:),'LineWidth',2);
  
                    end
                    xlim([xvalue(1),xvalue(end)]);ylim([min(min(ROINorm{nbplane}(:,:,ii))) max(max(ROINorm{nbplane}(:,:,ii)))]);
                    xlabel('Frame nb');

                    subplot(2,8,11); hold on; title('CR Probe'); plot(xvalue,(ROINorm{nbplane}(CRProbe,:,ii).'),'Color',colorsMCRl(2,:),'LineWidth',0.5); plot(xvalue,median(ROINorm{nbplane}(CRProbe,:,ii),1),'Color',colorsMCR(2,:),'LineWidth',2);
                    xlim([xvalue(1),xvalue(end)]);ylim([min(min(ROINorm{nbplane}(:,:,ii))) max(max(ROINorm{nbplane}(:,:,ii)))]);
                    xlabel('Frame nb');
                end


                subplot(2,8,5); hold on; title('Hit Reinf.'); imagesc(ROINorm{nbplane}(HitConca,:,ii));  axis tight; xlabel('Frame nb'); ylabel('Trial nb'); xline(15,'-','color','w', 'linewidth', 2);

                subplot(2,8,6); hold on;  title('Miss Reinf.');imagesc(ROINorm{nbplane}(MissConca,:,ii)); axis tight; xlabel('Frame nb');xline(framesbaseline,'-','color','w', 'linewidth', 2);

                if ~isempty(ProbeConca)
                    subplot(2,8,7); hold on;  title('Hit Probe'); imagesc(ROINorm{nbplane}(HitProbe,:,ii)); axis tight; xlabel('Frame nb');xline(framesbaseline,'-','color','w', 'linewidth', 2);

                    subplot(2,8,8); hold on;  title('Miss Probe'); imagesc(ROINorm{nbplane}(MissProbe,:,ii)); axis tight; xlabel('Frame nb');xline(framesbaseline,'-','color','w', 'linewidth', 2);
  
                end
                subplot(2,8,14); hold on; title('FA Reinf.');imagesc(ROINorm{nbplane}(FAConca,:,ii)); axis tight; xlabel('Frame nb');xline(framesbaseline,'-','color','w', 'linewidth', 2);

                subplot(2,8,13); hold on;  title('CR Reinf.');imagesc(ROINorm{nbplane}(CRConca,:,ii)); axis tight; xlabel('Frame nb');xline(framesbaseline,'-','color','w', 'linewidth', 2);ylabel('Trial nb');

                if ~isempty(ProbeConca)
                    subplot(2,8,16); hold on; title('FA Probe');imagesc(ROINorm{nbplane}(FAProbe,:,ii)); axis tight;xlabel('Frame nb');xline(framesbaseline,'-','color','w', 'linewidth', 2);

                    subplot(2,8,15); hold on;  title('CR Probe');imagesc(ROINorm{nbplane}(CRProbe,:,ii)); axis tight;   xlabel('Frame nb');xline(framesbaseline,'-','color','w', 'linewidth', 2);
                end

                suptitle(['Plane:' num2str(nbplane) ' Cell:' num2str(ii)])
                saveas(fig1,fullfile([behaviour_folder filesep 'evokedresponsebycell\'],['trial average per cell plane' num2str(nbplane) 'cell' num2str(ii) '.jpg']));
                close(fig1); 

            end
        end

            fig3=figure; hold on; %
            toneOnset= 0;
            
            Hcells=squeeze(mean(ROINorm{nbplane}(HitConca,:,:),1));
            Mcells=squeeze(mean(ROINorm{nbplane}(MissConca,:,:),1));
            FAcells=squeeze(mean(ROINorm{nbplane}(FAConca,:,:),1));
            CRcells=squeeze(mean(ROINorm{nbplane}(CRConca,:,:),1));
            Tcells=squeeze(mean(ROINorm{nbplane}(TargetIndConca,:,:),1));
            Fcells=squeeze(mean(ROINorm{nbplane}(FoilindConca,:,:),1));
            if ~isempty(ProbeConca)
                HcellsP=squeeze(mean(ROINorm{nbplane}(HitProbe,:,:),1));
                McellsP=squeeze(mean(ROINorm{nbplane}(MissProbe,:,:),1));
                FAcellsP=squeeze(mean(ROINorm{nbplane}(FAProbe,:,:),1));
                CRcellsP=squeeze(mean(ROINorm{nbplane}(CRProbe,:,:),1));
                TcellsP=squeeze(mean(ROINorm{nbplane}(TargetProbe,:,:),1));
                FcellsP=squeeze(mean(ROINorm{nbplane}(FoilProbe,:,:),1));
            end
            Hcells_tosave{nbplane}=Hcells;
            Mcells_tosave{nbplane}=Mcells;
            FAcells_tosave{nbplane}=FAcells;
            CRcells_tosave{nbplane}=CRcells;
            Tcells_tosave{nbplane}=Tcells;
            Fcells_tosave{nbplane}=Fcells;
            if ~isempty(ProbeConca)
                HcellsP_tosave{nbplane}=HcellsP;
                McellsP_tosave{nbplane}=McellsP;
                FAcellsP_tosave{nbplane}=FAcellsP;
                CRcellsP_tosave{nbplane}=CRcellsP;
                TcellsP_tosave{nbplane}=TcellsP;
                FcellsP_tosave{nbplane}=FcellsP;
            end
            
            nNeurosperplane{nbplane}=size(Hcells,2);
                        SEM1 = nanstd(Hcells'/sqrt(size(Hcells,2)));
            subplot(2,4,1); hold on; title('Hit Reinf.'); shadedErrorBar(xvalue,mean(Hcells,2),SEM1,'lineProps',{'Color',colorsHFA(1,:),'LineWidth',2});
          xlim([xvalue(1),xvalue(end)]);ylim(mylim); 
            xlabel('Frame nb'); ylabel('f/f0');
           
            if ~isempty(MissConca)
%                 subplot(2,4,2); hold on; title('Miss Reinf.');  plot(xvalue,mean(Mcells,2),'Color',colorsMCR(1,:),'LineWidth',2);
                SEM2 = nanstd(Mcells'/sqrt(size(Mcells,2)));
                subplot(2,4,2); hold on; title('Miss Reinf.'); shadedErrorBar(xvalue,mean(Mcells,2),SEM2,'lineProps',{'Color',colorsMCR(1,:),'LineWidth',2});
%                 plot(xvalue,Mcells','Color',colorsMCRl(1,:));
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
            else
                subplot(2,4,2); hold on; title('Miss Reinf.');
            end
            if ~isempty(ProbeConca)
%                 subplot(2,4,3); hold on; title('Hit Probe'); plot(xvalue,mean(HcellsP,2),'Color',colorsHFA(1,:),'LineWidth',2);
                SEM3 = nanstd(HcellsP'/sqrt(size(HcellsP,2)));
                subplot(2,4,3); hold on; title('Hit Probe'); shadedErrorBar(xvalue,mean(HcellsP,2),SEM3,'lineProps',{'Color',colorsHFA(1,:),'LineWidth',2});
%                 plot(xvalue,HcellsP','Color',colorsHFAl(1,:), 'LineWidth',0.5);
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
%                 subplot(2,4,4); hold on; title('Miss Probe'); plot(xvalue,mean(McellsP,2),'Color',colorsMCR(1,:),'LineWidth',2);
                SEM4 = nanstd(McellsP'/sqrt(size(McellsP,2)));
                subplot(2,4,4); hold on; title('Miss Probe'); shadedErrorBar(xvalue,mean(McellsP,2),SEM4,'lineProps',{'Color',colorsMCR(1,:),'LineWidth',2});
%                 plot(xvalue,McellsP','Color',colorsMCRl(1,:),'LineWidth',0.5);
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
            end
%             subplot(2,4,6); hold on; title('FA Reinf.');  plot(xvalue,mean(FAcells,2),'Color',colorsHFA(2,:),'LineWidth',2);
            SEM6 = nanstd(FAcells'/sqrt(size(FAcells,2)));
            subplot(2,4,6); hold on; title('FA Reinf.');  shadedErrorBar(xvalue,mean(FAcells,2),SEM6,'lineProps',{'Color',colorsHFA(2,:),'LineWidth',2});
%             plot(xvalue,FAcells','Color',colorsHFAl(2,:),'LineWidth',0.5);
            xlim([xvalue(1),xvalue(end)]);ylim(mylim);
            xlabel('Frame nb');
%             subplot(2,4,5); hold on; title('CR Reinf.');  plot(xvalue,mean(CRcells,2),'Color',colorsMCR(2,:),'LineWidth',2);
            SEM5 = nanstd(CRcells'/sqrt(size(CRcells,2)));
            subplot(2,4,5); hold on; title('CR Reinf.');  shadedErrorBar(xvalue,mean(CRcells,2), SEM5,'lineProps',{'Color',colorsMCR(2,:),'LineWidth',2});
%             plot(xvalue,CRcells','Color',colorsMCRl(2,:),'LineWidth',0.5);
            xlim([xvalue(1),xvalue(end)]);ylim(mylim);
            xlabel('Frame nb');ylabel('f/f0');
            if ~isempty(ProbeConca)
                if isempty(FAProbe)
                    subplot(2,4,8); hold on;
                else
%                     subplot(2,4,8); hold on; title('FA Probe');  plot(xvalue,mean(FAcellsP,2),'Color',colorsHFA(2,:),'LineWidth',2);
                    SEM8 = nanstd(FAcellsP'/sqrt(size(FAcellsP,2)));
                    subplot(2,4,8); hold on; title('FA Probe');  shadedErrorBar(xvalue,mean(FAcellsP,2),SEM8,'lineProps',{'Color',colorsHFA(2,:),'LineWidth',2});
%                     plot(xvalue,FAcellsP','Color',colorsHFAl(2,:),'LineWidth',0.5);
                end
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
%                 subplot(2,4,7); hold on; title('CR Probe');  plot(xvalue,mean(CRcellsP,2),'Color',colorsMCR(2,:),'LineWidth',2);
                SEM7 = nanstd(CRcellsP'/sqrt(size(CRcellsP,2)));
                subplot(2,4,7); hold on; title('CR Probe');  shadedErrorBar(xvalue,mean(CRcellsP,2), SEM7,'lineProps',{'Color',colorsMCR(2,:),'LineWidth',2});
%                 plot(xvalue,CRcellsP','Color',colorsMCRl(2,:),'LineWidth',0.5);
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
            end
            
            for nbplot=1:8
                subplot(2,4,nbplot);xline(toneOnset,'-'); %,{'Tone presentation'}
            end
                
        
            fig5=figure; hold on; %


            % Plot activity per outcome
%             xvalue = -15:53;
            toneOnset= 0;
%             mylim=[-0.02 0.06];
                 
            SEM1 = nanstd(Hcells'/sqrt(size(Hcells,2)));
            subplot(2,8,1); hold on; title('Hit Reinf.'); shadedErrorBar(xvalue,mean(Hcells,2),SEM1,'lineProps',{'Color',colorsHFA(1,:),'LineWidth',2});
            xlim([xvalue(1),xvalue(end)]);ylim(mylim); 
            xlabel('Frame nb'); ylabel('f/f0');
  
            
            if ~isempty(MissConca)
                SEM2 = nanstd(Mcells'/sqrt(size(Mcells,2)));
                subplot(2,8,2); hold on; title('Miss Reinf.'); shadedErrorBar(xvalue,mean(Mcells,2),SEM2,'lineProps',{'Color',colorsMCR(1,:),'LineWidth',2});
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
            else
                subplot(2,8,2); hold on; title('Miss Reinf.');
            end
            if ~isempty(ProbeConca)
                SEM3 = nanstd(HcellsP'/sqrt(size(HcellsP,2)));
                subplot(2,8,3); hold on; title('Hit Probe'); shadedErrorBar(xvalue,mean(HcellsP,2),SEM3,'lineProps',{'Color',colorsHFA(1,:),'LineWidth',2});
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
                SEM4 = nanstd(McellsP'/sqrt(size(McellsP,2)));
                subplot(2,8,4); hold on; title('Miss Probe'); shadedErrorBar(xvalue,mean(McellsP,2),SEM4,'lineProps',{'Color',colorsMCR(1,:),'LineWidth',2});
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
            end
            SEM6 = nanstd(FAcells'/sqrt(size(FAcells,2)));
            subplot(2,8,10); hold on; title('FA Reinf.');  shadedErrorBar(xvalue,mean(FAcells,2),SEM6,'lineProps',{'Color',colorsHFA(2,:),'LineWidth',2});
            xlim([xvalue(1),xvalue(end)]);ylim(mylim);
            xlabel('Frame nb');
            SEM5 = nanstd(CRcells'/sqrt(size(CRcells,2)));
            subplot(2,8,9); hold on; title('CR Reinf.');  shadedErrorBar(xvalue,mean(CRcells,2), SEM5,'lineProps',{'Color',colorsMCR(2,:),'LineWidth',2});
            xlim([xvalue(1),xvalue(end)]);ylim(mylim);
            xlabel('Frame nb');ylabel('f/f0');
            if ~isempty(ProbeConca)
                if isempty(FAProbe)
                    subplot(2,8,12); hold on;
                else
                    SEM8 = nanstd(FAcellsP'/sqrt(size(FAcellsP,2)));
                    subplot(2,8,12); hold on; title('FA Probe');  shadedErrorBar(xvalue,mean(FAcellsP,2),SEM8,'lineProps',{'Color',colorsHFA(2,:),'LineWidth',2});
                end
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
                SEM7 = nanstd(CRcellsP'/sqrt(size(CRcellsP,2)));
                subplot(2,8,11); hold on; title('CR Probe');  shadedErrorBar(xvalue,mean(CRcellsP,2), SEM7,'lineProps',{'Color',colorsMCR(2,:),'LineWidth',2});
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');
            end
            
            for nbplot=1:16
                subplot(2,8,nbplot);xline(toneOnset,'-'); %,{'Tone presentation'}
            end
            
            
            subplot(2,8,5); hold on; title('Hit');[~,sortIndex] = sort(mean(Hcells(:,:)),'descend'); imagesc(Hcells(:,sortIndex)');  
            axis tight; xlabel('Frame nb'); ylabel('Neuron'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
            if ~isnan(Hcells)
                caxis([prctile(Hcells(:),5) prctile(Hcells(:),95)]);
            end
            subplot(2,8,6); hold on;  title('Miss');[~,sortIndex] = sort(mean(Mcells(:,:)),'descend'); imagesc(Mcells(:,sortIndex)');  
            axis tight; xlabel('Frame nb'); ylabel('Neuron'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
            if ~isnan(Mcells)
                caxis([prctile(Mcells(:),5) prctile(Mcells(:),95)]);
            end
            if ~isempty(ProbeConca)
            subplot(2,8,7); hold on;  title('Hit probe');[~,sortIndex] = sort(mean(HcellsP(:,:)),'descend'); imagesc(HcellsP(:,sortIndex)');  
            axis tight; xlabel('Frame nb'); ylabel('Neuron'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
            
            if ~isnan(HcellsP)
                caxis([prctile(HcellsP(:),5) prctile(HcellsP(:),95)]);
            end
            
            subplot(2,8,8); hold on;  title('Miss probe'); [~,sortIndex] = sort(mean(McellsP(:,:)),'descend');imagesc(McellsP(:,sortIndex)');  
            axis tight; xlabel('Frame nb'); ylabel('Neuron'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
            if ~isnan(McellsP)
                caxis([prctile(McellsP(:),5) prctile(McellsP(:),95)]);
            end
            end
            subplot(2,8,14); hold on; title('FA');[~,sortIndex] = sort(mean(FAcells(:,:)),'descend');imagesc(FAcells(:,sortIndex)');  
            axis tight; xlabel('Frame nb'); ylabel('Neuron'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
            if ~isnan(FAcells)
                caxis([prctile(FAcells(:),5) prctile(FAcells(:),95)]);
            end
            subplot(2,8,13); hold on;  title('CR');[~,sortIndex] = sort(mean(CRcells(:,:)),'descend');imagesc(CRcells(:,sortIndex)');  
            axis tight; xlabel('Frame nb'); ylabel('Neuron'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
            if ~isnan(CRcells)
                caxis([prctile(CRcells(:),5) prctile(CRcells(:),95)]);
            end
            if ~isempty(ProbeConca)
            subplot(2,8,16); hold on; title('FA probe');[~,sortIndex] = sort(mean(FAcellsP(:,:)),'descend');imagesc(FAcellsP(:,sortIndex)');  
            axis tight; xlabel('Frame nb'); ylabel('Neuron'); xline(framesbaseline,'-','color','w', 'linewidth', 2); 
            if ~isnan(FAcellsP)
                caxis([prctile(FAcellsP(:),5) prctile(FAcellsP(:),95)]);
            end
            
            subplot(2,8,15); hold on;  title('CR probe');[~,sortIndex] = sort(mean(CRcellsP(:,:)),'descend');imagesc(CRcellsP(:,sortIndex)');  
            axis tight; xlabel('Frame nb'); ylabel('Neuron'); xline(framesbaseline,'-','color','w', 'linewidth', 2); 
            if ~isnan(CRcellsP)
                caxis([prctile(CRcellsP(:),5) prctile(CRcellsP(:),95)]);
            end
            end
            suptitle(['Mouse:' animal_list{fl} ' Plane:' num2str(nbplane) ' Neurons:' num2str(size(Hcells,2))])
%             saveas(fig5,fullfile([behaviour_folder filesep 'evokedresponsebycell\'],['mean response per plane raster' num2str(nbplane) '.jpg']));
%             close(fig5); 
                        
           if distancetoplaqueanalysis==1 
               fig6=figure; hold on; %


                % Plot activity per outcome
    %             xvalue = -15:53;
                toneOnset= 0;
    %             mylim=[-0.02 0.06];

%               less100=find(NeuronDist(nbplane).plane<100);
               less100=find(NeuronDist(nbplane).plane<60);
               less200=find(NeuronDist(nbplane).plane>=60 & NeuronDist(nbplane).plane<200);
               more200=find(NeuronDist(nbplane).plane>=200);
               
               less60{nbplane}= less100;
               morethan60{nbplane}=less200;

               Hcells100=Hcells(:,less100);
               Hcells200=Hcells(:,less200);
               Hcellsmore=Hcells(:,more200);

               Mcells100=Mcells(:,less100);
               Mcells200=Mcells(:,less200);
               Mcellsmore=Mcells(:,more200);            

               FAcells100=FAcells(:,less100);
               FAcells200=FAcells(:,less200);
               FAcellsmore=FAcells(:,more200);

               CRcells100=CRcells(:,less100);
               CRcells200=CRcells(:,less200);
               CRcellsmore=CRcells(:,more200);
               if ~isempty(ProbeConca)
               HPcells100=HcellsP(:,less100);
               HPcells200=HcellsP(:,less200);
               HPcellsmore=HcellsP(:,more200);

               MPcells100=McellsP(:,less100);
               MPcells200=McellsP(:,less200);
               MPcellsmore=McellsP(:,more200);            

               FAPcells100=FAcellsP(:,less100);
               FAPcells200=FAcellsP(:,less200);
               FAPcellsmore=FAcellsP(:,more200);

               CRPcells100=CRcellsP(:,less100);
               CRPcells200=CRcellsP(:,less200);
               CRPcellsmore=CRcellsP(:,more200);
               end




                SEM1_100H = nanstd(Hcells100'/sqrt(size(Hcells100,2)));
                SEM1_200H = nanstd(Hcells200'/sqrt(size(Hcells200,2)));
                SEM1_moreH = nanstd(Hcellsmore'/sqrt(size(Hcellsmore,2)));
                subplot(2,8,1); hold on; title('Hit Reinf.'); 
                shadedErrorBar(xvalue,mean(Hcells100,2)',SEM1_100H,'lineProps',{'Color','k','LineWidth',2});
                shadedErrorBar(xvalue,mean(Hcells200,2)',SEM1_200H,'lineProps',{'Color',colorsHFA(1,:),'LineWidth',2});
    %             shadedErrorBar(xvalue,mean(Hcellsmore,2)',SEM1_100,'lineProps',{'Color',colorsHFAl(1,:),'LineWidth',2});
                xlim([xvalue(1),xvalue(end)]);ylim(mylim); 
                xlabel('Frame nb'); ylabel('f/f0');


                if ~isempty(MissConca)
                    SEM2_100M = nanstd(Mcells100'/sqrt(size(Mcells100,2)));
                    SEM2_200M = nanstd(Mcells200'/sqrt(size(Mcells200,2)));
                    subplot(2,8,2); hold on; title('Miss Reinf.'); 
                    shadedErrorBar(xvalue,mean(Mcells100,2),SEM2_100M,'lineProps',{'Color','k','LineWidth',2});
                    shadedErrorBar(xvalue,mean(Mcells200,2),SEM2_200M,'lineProps',{'Color',colorsMCR(1,:),'LineWidth',2});
                    xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                    xlabel('Frame nb');
                else
                    subplot(2,8,2); hold on; title('Miss Reinf.');
                end
                if ~isempty(ProbeConca)
                    SEM3_100HP = nanstd(HPcells100'/sqrt(size(HPcells100,2)));
                    SEM3_200HP = nanstd(HPcells200'/sqrt(size(HPcells200,2)));
                    subplot(2,8,3); hold on; title('Hit Probe'); 
                    shadedErrorBar(xvalue,mean(HPcells100,2),SEM3_100HP,'lineProps',{'Color','k','LineWidth',2});
                    shadedErrorBar(xvalue,mean(HPcells200,2),SEM3_200HP,'lineProps',{'Color',colorsHFA(1,:),'LineWidth',2});
                    xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                    xlabel('Frame nb');

                    SEM4_100MP = nanstd(MPcells100'/sqrt(size(MPcells100,2)));
                    SEM4_200MP = nanstd(MPcells200'/sqrt(size(MPcells200,2)));
                    subplot(2,8,4); hold on; title('Miss Probe'); 
                    shadedErrorBar(xvalue,mean(MPcells100,2),SEM4_100MP,'lineProps',{'Color','k','LineWidth',2});
                    shadedErrorBar(xvalue,mean(MPcells200,2),SEM4_200MP,'lineProps',{'Color',colorsMCR(1,:),'LineWidth',2});
                    xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                    xlabel('Frame nb');
                end
                SEM6_100FA = nanstd(FAcells100'/sqrt(size(FAcells100,2)));
                SEM6_200FA = nanstd(FAcells200'/sqrt(size(FAcells200,2)));
                subplot(2,8,10); hold on; title('FA Reinf.');  
                shadedErrorBar(xvalue,mean(FAcells100,2),SEM6_100FA,'lineProps',{'Color','k','LineWidth',2});
                shadedErrorBar(xvalue,mean(FAcells200,2),SEM6_200FA,'lineProps',{'Color',colorsHFA(2,:),'LineWidth',2});
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');

                SEM5_100CR = nanstd(CRcells100'/sqrt(size(CRcells100,2)));
                SEM5_200CR = nanstd(CRcells200'/sqrt(size(CRcells200,2)));
                subplot(2,8,9); hold on; title('CR Reinf.');  
                shadedErrorBar(xvalue,mean(CRcells100,2), SEM5_100CR,'lineProps',{'Color','k','LineWidth',2});
                shadedErrorBar(xvalue,mean(CRcells200,2), SEM5_200CR,'lineProps',{'Color',colorsMCR(2,:),'LineWidth',2});         
                xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                xlabel('Frame nb');ylabel('f/f0');
                if ~isempty(ProbeConca)
                    if isempty(FAProbe)
                        subplot(2,8,12); hold on;
                    else
                        SEM8_100CR = nanstd(FAPcells100'/sqrt(size(FAPcells100,2)));
                        SEM8_200CR = nanstd(FAPcells200'/sqrt(size(FAPcells200,2)));
                        subplot(2,8,12); hold on; title('FA Probe');  
                        shadedErrorBar(xvalue,mean(FAPcells100,2),SEM8_100CR,'lineProps',{'Color','k','LineWidth',2});
                        shadedErrorBar(xvalue,mean(FAPcells200,2),SEM8_200CR,'lineProps',{'Color',colorsHFA(2,:),'LineWidth',2});
                    end
                    xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                    xlabel('Frame nb');
                    SEM7_100CRP = nanstd(CRPcells100'/sqrt(size(CRPcells100,2)));
                    SEM7_200CRP = nanstd(CRPcells200'/sqrt(size(CRPcells200,2)));
                    subplot(2,8,11); hold on; title('CR Probe');  
                    shadedErrorBar(xvalue,mean(CRPcells100,2), SEM7_100CRP,'lineProps',{'Color','k','LineWidth',2});
                    shadedErrorBar(xvalue,mean(CRPcells200,2), SEM7_200CRP,'lineProps',{'Color',colorsMCR(2,:),'LineWidth',2});
                    xlim([xvalue(1),xvalue(end)]);ylim(mylim);
                    xlabel('Frame nb');
                end

                for nbplot=1:16
                    subplot(2,8,nbplot);xline(toneOnset,'-'); %,{'Tone presentation'}
                end


                subplot(2,8,5); hold on; title('Hit');[~,sortIndex] = sort(mean(Hcells(:,:)),'descend'); imagesc(Hcells(:,sortIndex)');  
                axis tight; xlabel('Frame nb'); ylabel('Neuron (distance to plaque)'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
                if ~isnan(Hcells)
                caxis([prctile(Hcells(:),5) prctile(Hcells(:),95)]);
                end
                subplot(2,8,6); hold on;  title('Miss');[~,sortIndex] = sort(mean(Mcells(:,:)),'descend'); imagesc(Mcells(:,sortIndex)');  
                axis tight; xlabel('Frame nb'); ylabel('Neuron (distance to plaque)'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
                if ~isnan(Mcells)
                    caxis([prctile(Mcells(:),5) prctile(Mcells(:),95)]);
                end
                if ~isempty(ProbeConca)
                subplot(2,8,7); hold on;  title('Hit probe');[~,sortIndex] = sort(mean(HcellsP(:,:)),'descend'); imagesc(HcellsP(:,sortIndex)');  
                axis tight; xlabel('Frame nb'); ylabel('Neuron (distance to plaque)'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
                if ~isempty(HcellsP)
                    caxis([prctile(HcellsP(:),5) prctile(HcellsP(:),95)]);
                end
                
                subplot(2,8,8); hold on;  title('Miss probe'); [~,sortIndex] = sort(mean(McellsP(:,:)),'descend');imagesc(McellsP(:,sortIndex)');  
                axis tight; xlabel('Frame nb'); ylabel('Neuron (distance to plaque)'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
                if ~isnan(McellsP)
                    caxis([prctile(McellsP(:),5) prctile(McellsP(:),95)]);
                end
                end
                subplot(2,8,14); hold on; title('FA');[~,sortIndex] = sort(mean(FAcells(:,:)),'descend');imagesc(FAcells(:,sortIndex)');  
                axis tight; xlabel('Frame nb'); ylabel('Neuron (distance to plaque)'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
                if ~isempty(FAcells)
                    caxis([prctile(FAcells(:),5) prctile(FAcells(:),95)]);
                end
                subplot(2,8,13); hold on;  title('CR');[~,sortIndex] = sort(mean(CRcells(:,:)),'descend');imagesc(CRcells(:,sortIndex)');  
                axis tight; xlabel('Frame nb'); ylabel('Neuron (distance to plaque)'); xline(framesbaseline,'-','color','w', 'linewidth', 2);
                if ~isempty(CRcells)
                    caxis([prctile(CRcells(:),5) prctile(CRcells(:),95)]);
                end
                if ~isempty(ProbeConca)
                subplot(2,8,16); hold on; title('FA probe');[~,sortIndex] = sort(mean(FAcellsP(:,:)),'descend');imagesc(FAcellsP(:,sortIndex)');  
                axis tight; xlabel('Frame nb'); ylabel('Neuron (distance to plaque)'); xline(framesbaseline,'-','color','w', 'linewidth', 2); 
                if ~isnan(FAcellsP)
                    caxis([prctile(FAcellsP(:),5) prctile(FAcellsP(:),95)]);
                end
                subplot(2,8,15); hold on;  title('CR probe');[~,sortIndex] = sort(mean(CRcellsP(:,:)),'descend');imagesc(CRcellsP(:,sortIndex)');  
                axis tight; xlabel('Frame nb'); ylabel('Neuron (distance to plaque)'); xline(framesbaseline,'-','color','w', 'linewidth', 2); 
                if ~isempty(CRcellsP)
                    caxis([prctile(CRcellsP(:),5) prctile(CRcellsP(:),95)]);
                end
                end
                suptitle(['Mouse:' animal_list{fl} ' Plane:' num2str(nbplane) ' Neurons:' num2str(size(Hcells,2))])
    %             saveas(fig6,fullfile([behaviour_folder filesep 'evokedresponsebycell\'],['mean response per plane distance to plaque' num2str(nbplane) '.jpg']));
    %             close(fig6); 
           end
           

    end 
    
    cumsumneurons = cumsum([0 cellfun(@(x) size(x,3),ToneEvokedMatConca)]);
    totalneurons=cumsumneurons(end);       
% disp(nbplane);
    if distancetoplaqueanalysis==1   
        neurondistindexunder60=less60;
        neurondistindex60orhigher=morethan60;

        cumsumneurons(end)=[];
        NeuronIndless60 = cell2mat(cellfun(@(x,y) (x+y)',neurondistindexunder60,mat2cell(cumsumneurons,1,ones(1,length(cumsumneurons))),'UniformOutput',0));
        NeuronIndmore60 = cell2mat(cellfun(@(x,y) (x+y)',neurondistindex60orhigher,mat2cell(cumsumneurons,1,ones(1,length(cumsumneurons))),'UniformOutput',0));
    else
        NeuronIndless60 = [];
        NeuronIndmore60 = [];
    end
    DataImagingE(fl).SessionNb = 1;
    DataImagingE(fl).Mouse = cell2mat(animal_list(fl));
%     DataImagingE(fl).ImagingFileName = behaviorfilesession;
    DataImagingE(fl).BehaviorFileName = iwant;
    DataImagingE(fl).ToneEvoked = ToneEvokedMat;
    DataImagingE(fl).Ttest = testresults;
    DataImagingE(fl).Hit = Hit;
    DataImagingE(fl).Miss = Miss;
    DataImagingE(fl).Falarm = Falarm;
    DataImagingE(fl).Crejection = Crejection;
    DataImagingE(fl).Probeind = Probeind;
    DataImagingE(fl).Targetind = Targetind;
    DataImagingE(fl).Foilind = Foilind;
    DataImagingE(fl).ToneEvokedMatConca = ToneEvokedMatConca;
    DataImagingE(fl).ToneEvokedMatConcaRAW = ToneEvokedMatrawConca;
    DataImagingE(fl).RewardEvokedMat = RewardEvokedMat;
    DataImagingE(fl).SelectedCells = SelectedCells_save; %DataImagingE(snb);
    DataImagingE(fl).TargetIndConca = TargetIndConca;
    DataImagingE(fl).FoilIndConca = FoilindConca;
    DataImagingE(fl).ProbeIndConca = ProbeConca;
    if ~isempty(ProbeConca)
    DataImagingE(fl).TargetProbe =TargetProbe;
    DataImagingE(fl).FoilProbe=FoilProbe;
    end
    DataImagingE(fl).BlockId = BlockId;
    DataImagingE(fl).HitConca = HitConca;
    DataImagingE(fl).MissConca = MissConca;
    DataImagingE(fl).FAConca = FAConca;
    DataImagingE(fl).CRConca = CRConca;
%     DataImagingE(fl).RawTraceMat = TraceMat;
    if ~isempty(ProbeConca)
        DataImagingE(fl).HitProbe = HitProbe;
        DataImagingE(fl).MissProbe = MissProbe;
        DataImagingE(fl).FAProbe = FAProbe;
        DataImagingE(fl).CRProbe = CRProbe;
    else
        DataImagingE(fl).HitProbe = [];
        DataImagingE(fl).MissProbe = [];
        DataImagingE(fl).FAProbe = [];
        DataImagingE(fl).CRProbe = [];
    end
    DataImagingE(fl).LickFrameConca = LickFrameConca;
    DataImagingE(fl).ROINorm = ROINorm;
    if distancetoplaqueanalysis==1
        DataImagingE(fl).NeuronDist = NeuronDist;
    else
        DataImagingE(fl).NeuronDist = [];
    end
    
    DataImagingE(fl).Hcells = Hcells_tosave;
    DataImagingE(fl).Mcells = Mcells_tosave;
    DataImagingE(fl).FAcells = FAcells_tosave;
    DataImagingE(fl).CRcells = CRcells_tosave;
    DataImagingE(fl).Tcells = Tcells_tosave;
    DataImagingE(fl).Fcells = Fcells_tosave;
    if ~isempty(ProbeConca)
        DataImagingE(fl).HcellsP = HcellsP_tosave;
        DataImagingE(fl).McellsP = McellsP_tosave;
        DataImagingE(fl).FAcellsP = FAcellsP_tosave;
        DataImagingE(fl).CRcellsP = CRcellsP_tosave;  
        DataImagingE(fl).TcellsP = TcellsP_tosave;
        DataImagingE(fl).FcellsP = FcellsP_tosave;  
    else
        DataImagingE(fl).HcellsP = [];
        DataImagingE(fl).McellsP = [];
        DataImagingE(fl).FAcellsP = [];
        DataImagingE(fl).CRcellsP = []; 
        DataImagingE(fl).TcellsP = []; 
        DataImagingE(fl).FcellsP = [];  
    end
    DataImagingE(fl).nPlanes =  nPlanes;
    DataImagingE(fl).nNeurosperplane =  nNeurosperplane;
    DataImagingE(fl).totalneurons =  totalneurons;
    DataImagingE(fl).NeuronIndless60 =  NeuronIndless60;
    DataImagingE(fl).NeuronIndmore60 =  NeuronIndmore60;
    
    save([ behaviour_folder filesep 'summarydata\DataImagingE'],'DataImagingE');
% 
 end

