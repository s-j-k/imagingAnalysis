%Guide to running Calman on the 920 and 820 data
%BBS 1/6/221

%Download matlab
%Parallel processing toolbox
%Signal processing toolbox

%Install CALMAN:
%Install NORMCORRE: https://github.com/flatironinstitute/NoRMCorre

%STEP1 Open 920 and 820 stacks sequence with FIJI, see note 1 below
%STEP2 Image-> adjust -> size ->256 by 256 (downsample to 256 to save space)
%STEP3 Save 820 and 920 tif in new folder, 920 should be in a new folder (see step 4).

%STEP4 Run run_pipeline on the new folder
%STEP5 Run run_pipeline2 on a folder containing data from 820, this will
%extract the fluorescence data from the mask generated by CALMAN see note 2 below




%Note 1, related to STEP1:
%you should be able to run this pipeline on a 256 by 256 movie up to 30K frames, but start with 10K
%for debugging.  If you go over 30K frames an error occurs running normcorre:
%"Error using normcorre_batch_even (line 111)
%Movie appears to have only one frame. Use the function normcorre instead"
%seems to be a problem with mfinfo when file gets too big. probably best to
%run on 30000K frames. this is 10000 seconds or ~16min 40 sec, which seems
%fine.

%Note 2, related to STEP5:
% to get the mask for the nth cell you need to do the following:
%mask=reshape(A_keep(:,n),256,256);

