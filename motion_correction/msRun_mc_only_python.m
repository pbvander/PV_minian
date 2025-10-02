%%Foreword by PV
%This code is meant to be run in Python as part of the MinAn pipeline
%(replacing the native MiniAn motion correction step)
%%
function msRun_mc_only_python(dpath, vid_path, isnonrigid, grid_size, spatial_downsampling, bin_width, mot_uf, correct_bidir, overlap_pre, overlap_post, max_shift);
%Set defaults for arguments
arguments
   dpath; %no default setting
   vid_path = [dpath '\\' 'msRun' '\\' 'msvideo.avi'] % default [dpath '\\' 'msRun' '\\' 'msvideo.avi']
   isnonrigid = false; %default false
   grid_size = []; %default []   %Changing gride size can cause errors. [128,128] works with overlap=32 but many other combinations causes errors in NormCorre
   spatial_downsampling = 1;
   bin_width = 50;
   mot_uf = 4;
   correct_bidir = false;
   overlap_pre = 32; %default 32
   overlap_post = 32; %default 32
   max_shift = 20;
end
%% msRun2018
% Version 1.0 GE
% Updated version of the msRun script originally proposed by Daniel B
% Aharoni to analyse miniscope 1p calcium imaging data.
% This version is build on top of the original package to maximize compatibility.
% It includes NormCorre for image registration, CNMF-E for source extraction,
% and CellReg for chronic registration across sessions. It also includes
% custom written scripts to explore the data (eg spatial firing, transients
% properties visualization)

% Copyright (C) 2017-2018 by Guillaume Etter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or any
% later version.
% Contact: etterguillaume@gmail.com

%Add path to files (If this doesn't add folders to path correctly just do it manually)
addpath(genpath('C:/Users/General Correa Lab/Documents/GitHub/PV_minian/motion_correction/')); 

%
%% Auto-detect operating system
%clear
if ispc
    separator = '\'; % For pc operating systems
else
    separator = '/'; % For unix (mac, linux) operating systems
end
parpool_size = 8;

%% Parameters
%spatial_downsampling = 1; % (Recommended range: 2 - 4. Downsampling significantly increases computational speed, but verify it does not
%isnonrigid = True; % If true, performs non-rigid registration (slower). If false, rigid alignment (faster).
analyse_behavior = false;
copy_to_googledrive = false;
if copy_to_googledrive
    copydirpath = uigetdir([],'Please select the root folder in which files will be copied');
end
filepath = [dpath '\\'];
% Generate timestamp to save analysis
script_start = tic;
analysis_time =strcat(date,'_', num2str(hour(now)),'-',num2str(minute(now)),'-',num2str(floor(second(now))));

%% 1 - Create video object and save into matfile
disp('Step 1: Create video object');
ms = msGenerateVideoObj(filepath);
ms.analysis_time = analysis_time;
ms.ds = spatial_downsampling;
mkdir(strcat(filepath,separator,'msRun'));
save([ms.dirName separator 'ms.mat'],'ms');

%% 2 - Perform motion correction using NormCorre
disp('Step 2: Motion correction');

%Parameters for NormCorre
normcorre_options=struct('bin_width',bin_width, 'grid_size',grid_size, ... %Changing gride size to empty (default) fixed an error here!
    'mot_uf',mot_uf,'correct_bidir',correct_bidir, ...
    'overlap_pre',overlap_pre,'overlap_post',overlap_post,'max_shift',max_shift, 'vid_path',vid_path);

ms = msNormCorre(ms,isnonrigid, normcorre_options); %PV altered msNormCorre.m script so that it worked with different grid sizes
save([ms.dirName separator 'ms.mat'],'ms');
% sys=system('python "C:/Users/General Correa Lab/Box/correalab/Member Folders/Paul Vander/Code/SendTextMessage.py" 7349685923 verizon "Motion correction is done!"'); %sends a text when code is done running

%Write shifts to a file (PV added this 3/7/25)
for vid_i = 1:size(ms.shifts,2) %iterate through all videos
    vid_struct = ms.shifts{1,vid_i}; %get shift struct for the video
    vid_shifts = {vid_struct.shifts}; %extract shifts column from vid_struct
    num_frames = size(vid_shifts, 2); %get size of dimension 2 (frames) of vid_shifts array
    if vid_i==1 
        motion = cell(0, 2); %initialize array, if not already done
        motion{1,1} = 'width'; %add header
        motion{1,2} = 'height'; %add header
    end
    for frame_i = 1:num_frames %iterate through each frame
        motion{frame_i + ((vid_i-1)*1000) + 1, 1} = vid_shifts{frame_i}(2); %horizontal/width shift
        motion{frame_i + ((vid_i-1)*1000) + 1, 2} = vid_shifts{frame_i}(1); %vertical/height shift
    end
end
vid_path_split = strsplit(vid_path, "\\");
writetable(cell2table(motion), strjoin([vid_path_split(1:length(vid_path_split)-1) 'motion.csv'],"\\"), "WriteVariableNames",false) %write to a csv file
end 