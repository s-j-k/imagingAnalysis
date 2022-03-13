% % %run_pipeline_2_bbs.m
% % BBS 1/8/21
% % %This function performs the following operations
% % %1) motion corrects and registers the data collected at 820 to that collected at 920
% % %2) computes the average flourescence for each cell using the masks from
% % %   the full imaging session
% % before running this script make sure you have read the readme file and
% % have run the script run_pipeline BBS
% %% STEP 1: REGISTER THE GCAMP CHANNEL RECORDED AT 820nm
% % this uses the template created at 920 to register 820
% clear;
% 
% % STEP 1A: LOAD THE TEMPLATE MATRIX FROM THE .mat FILE YOU SAVED AFTER RUNNING run_pipeline_bbs, change the
% % path as needed:
path_for_920_data='\Users\sjkim1\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker\raw data\0104data\gcamp\data920';
% load(path_for_920_data,'template')
% % STEP 1B: CHANGE DIRECTORY TO THE FOLDER CONTAINING THE 820 TIF STACK
path_for_820_data='\Users\sjkim1\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker\raw data\0104data';
% cd(path_for_820_data)
% % STEP 1C: IDEFITY THE 820 TIF STACK
files='Full820.tif';
% % STEP 1D: RUN MOTION CORRECTION USING THE TEMPLATE FROM THE 920 DATA SET.

FOV = size(read_file(files));
motion_correct = true;                            % perform motion correction
non_rigid = true;                                 % flag for non-rigid motion correction
output_type = 'mat';                               % format to save registered files
append = '_rig';        % use this to save motion corrected files
%col_shift = [];
options_mc = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'grid_size',[128,128],'init_batch',200,...
    'overlap_pre',32,'mot_uf',4,'bin_width',200,'max_shift',24,'max_dev',8,'us_fac',50,...
    'output_type',output_type);
fullname = files;
[folder_name,file_name,~] = fileparts(fullname);
output_filename = fullfile(folder_name,[file_name,append,'.',output_type]);
options_mc = NoRMCorreSetParms(options_mc,'output_filename',output_filename,'h5_filename','','tiff_filename',''); % update output file name

template = []; 
[M_green_final,shifts,template,options_mc,col_shift] = normcorre_batch_even(fullname,options_mc,template);

% %% OPTIONAL, You can plot the registered FOVs.
figure;
subplot(1,2,1); imagesc(template); axis off;  title('ex920')
g820=squeeze(mean(M_green_final,3));
subplot(1,2,2); imagesc(g820); axis off;  title('ex820')
% 
% %% STEP 2: APPLY THE MASKS FROM CALMAN TO EXTRACT THE 820 FLUORESCENCE
load(path_for_920_data,'A_keep')
[~,N]=size(A_keep);
M_green_final=double(M_green_final);
load(path_for_920_data,'F_dff')
M_green_final=M_green_final(:,:,2001:5000);%reduce to end of imaging session
ISO_F=NaN(N,length(M_green_final));
%%
tic
figure;
time=(1:length(M_green_final))/30;
for n=1:N
    mask1=reshape(A_keep(:,n),256,256);
    mask2=full(mask1);
    mask3=repmat(mask2,1,1,length(M_green_final));
    trace=mean(mean(M_green_final.*mask3));
    ISO_F(n,:)=(trace-mean(trace))/mean(trace);
    toc
    
    plot(time(1:length(F_dff)),F_dff(n,:))
    hold on
    plot(time,ISO_F(n,:))    
    ylim([-0.5 3])
    xlabel('Sec')
    ylabel('DF/F')
    title(sprintf('Cell %d',n))
    legend('920','820')
    drawnow
    hold off
end

save(['DATA820.mat'])


%% OPTIONAL, You can plot the components to see what cells were selected.
%
% load(path_for_920_data)
sizY = data.sizY;

patch_size = [40,40];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 7;                                            % number of components to be found
tau = 8;                                          % std of gaussian kernel (half size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold


options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'deconv_method','constrained_foopsi',...    % neural activity deconvolution method
    'p',p,...                                   % order of calcium dynamics
    'ssub',2,...                                % spatial downsampling when processing
    'tsub',2,...                                % further temporal downsampling when processing
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'max_size_thr',300,'min_size_thr',10,...    % max/min acceptable size for each component
    'spatial_method','regularized',...          % method for updating spatial components
    'df_prctile',50,...                         % take the median of background fluorescence to compute baseline fluorescence 
    'fr',fr/tsub,...                            % downsamples
    'space_thresh',0.35,...                     % space correlation acceptance threshold
    'min_SNR',2.0,...                           % trace SNR acceptance threshold
    'cnn_thr',0.2,...                           % cnn classifier acceptance threshold
    'nb',1,...                                  % number of background components per patch
    'gnb',3,...                                 % number of global background components
    'decay_time',0.5...                         % length of typical transient for the indicator used
    );


figure;
subplot(1,2,1)
plot_contours(A_keep(:,58),g820,options,0,[],Coor_k,[],1); title('820','fontweight','bold','fontsize',14);
subplot(1,2,2)
plot_contours(A_keep(:,58),template,options,0,[],Coor_k,[],1); title('920','fontweight','bold','fontsize',14);
