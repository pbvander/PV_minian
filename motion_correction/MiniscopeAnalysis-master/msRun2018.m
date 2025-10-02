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

addpath(genpath('C:/Users/Paul/Box/correalab/Member Folders/Paul Vander/Code/')); % add this one first to use the right version of msGenerateVideoObj
addpath(genpath('D:\Miniscope\DataAnalysis\Emily_Wu\Miniscope\CNMF_E-master'));
addpath(genpath('D:\Miniscope\DataAnalysis\Emily_Wu\Miniscope\ffmpeg-r8'));
addpath(genpath('D:\Miniscope\DataAnalysis\Emily_Wu\Miniscope\Fiji.app'));
addpath(genpath('D:\Miniscope\DataAnalysis\Emily_Wu\Miniscope\MiniscopeAnalysis-master'));
addpath(genpath('D:\Miniscope\DataAnalysis\Emily_Wu\Miniscope\NoRMCorre-master'));
addpath('D:\Miniscope\DataAnalysis\Emily_Wu\Miniscope');

%% Auto-detect operating system
clear
if ispc
    separator = '\'; % For pc operating systems
else
    separator = '/'; % For unix (mac, linux) operating systems
end
parpool_size = 8;

%% Parameters
spatial_downsampling = 1; % (Recommended range: 2 - 4. Downsampling significantly increases computational speed, but verify it does not
isnonrigid = false; % If true, performs non-rigid registration (slower). If false, rigid alignment (faster).
analyse_behavior = true;
copy_to_googledrive = false;
if copy_to_googledrive
    copydirpath = uigetdir([],'Please select the root folder in which files will be copied');
end
filepath = 'D:\Miniscope\Prosocial_ACC\MZ67\20210716\ConcatenatedCrop';
% Generate timestamp to save analysis
script_start = tic;
analysis_time = 'msCam';

%% 1 - Create video object and save into matfile
disp('Step 1: Create video object');
ms = msGenerateVideoObj(filepath,'msCam');
ms.analysis_time = analysis_time;
ms.ds = spatial_downsampling;
mkdir(strcat(filepath,separator,analysis_time));
save([ms.dirName separator 'ms.mat'],'ms');

%% 2 - Perform motion correction using NormCorre
disp('Step 2: Motion correction');
ms = msNormCorre(ms,isnonrigid);
save([ms.dirName separator 'ms.mat'],'ms');

%% 3 - For ACC/PFC data, run dFFNormalize.m separately (cannot be run within another function) before running CNMFE (not needed for MeA data)

% %% Alternative: temporal downsampling to reduce video size - take the mean of every two frames
% temporal_ds(ms);
% ms.tds = 2; %temporal downsampling factor
% save([ms.dirName separator 'ms.mat'],'ms');

%% Cut movies for testing purposes
% video_cut(1, 18000, 'D:\Miniscope\Prosocial_ACC\MZ67\20210716\ConcatenatedCrop\msCam\dFF.avi');

%% 4 - Perform CNMFE
disp('Step 4: CNMFE');
load([filepath separator 'ms.mat']);

parpool_size = 2;
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(parpool_size)
elseif poolobj.NumWorkers ~= parpool_size
    delete(poolobj)
    parpool(parpool_size)
end

% Parameters for CNMFE (most likely to affecy results: merge_thr_spatial, min_pnr, gSiz, tsub)  
pars_envs = struct('memory_size_to_use', 12, ...   % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 0.6, ...   % GB, space for loading data within one patch
    'patch_dims', [64, 64]);  %patch size

deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'constrained', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100);    % maximum decay time (unit: frame);

CNMFE_options = struct(...
    'movie', get_fullname([ms.dirName separator ms.analysis_time separator 'dFF.avi']), ...
    'pars_envs', pars_envs, ...
    'include_residual', false,...
    'gSig', 3,... % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
    'gSiz', 12,... % pixel, neuron diameter
    'ssub', 1, ...
    'with_dendrites', false, ...
    'updateA_bSiz', 3, ...
    'updateA_dist', 3, ...
    'Fs', 30,... % frame rate
    'tsub', 5,... % temporal downsampling factor
    'deconv_flag', true, ...
    'deconv_options', deconv_options, ...
    'nk', 3,...
    'detrend_method', 'spline', ...
    ...% background model
    'bg_model', 'ring',... % model of the background {'ring', 'svd'(default), 'nmf'}
    'nb', 1,...             % number of background sources for each patch (only be used in SVD and NMF model)
    'ring_radius', 16,...  % when the ring model used, it is the radius of the ring used in the background model.
    ...% merge
    'show_merge', false, ...
    'merge_thr', 0.65,...
    'merge_thr_spatial', [0.425,0.1,-Inf],...% thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
    'dmin', 3,... % minimum distances between two neurons. it is used together with merge_thr
    'dmin_only', 2, ...
    'method_dist', 'max', ...
    ...% initialize
    'min_corr', 0.75,... % minimum local correlation for a seeding pixel, default 0.8, cmk 0.75
    'min_pnr', 21,... % minimum peak-to-noise ratio for a seeding pixel, cmk 21, gaba 12
    ...% residual
    'min_corr_res', 0.7,... % cmk 0.7 gaba 0.7
    'min_pnr_res', 19); % cmk 19 gaba 10

ms = msRunCNMFE_large(ms, CNMFE_options);
analysis_duration = toc(script_start);
ms.analysis_duration = analysis_duration;

ms.CNMFE_options = CNMFE_options;
save([ms.dirName separator 'ms.mat'],'ms','-v7.3');
disp(['Data analyzed in ' num2str(analysis_duration) 's']);

% msExtractSFPs(ms); % Extract spatial footprints for subsequent re-alignement
msPlotSFPs(ms); % Plot SFP outlines filled with different colors, generate mask in ImageJ to remove bright spots on the periphery mis-identified as cells
% save image with exactly the same dimension in order to create mask in ImageJ
I = getframe;
imwrite(I.cdata, [ms.dirName '\msCam\dFF_source_extraction\ROI_thr0.6.jpg']);

% remove bright spots on the periphery mis-identified as cells
mask = imread([ms.dirName separator 'msCam\dFF_source_extraction_1\Mask_ROI_thr0.6.jpg']); mask = mask(:, :, 1);
% mask = imread('D:\Miniscope\Prosocial_ACC\MZ67\20210716\ConcatenatedCrop\msCam\dFF_f1-18000_source_extraction\Crop_ROI_f1-18000.jpg');
mask = mask == 255; % 0
ms = rm_periphery(ms, mask);
title([num2str(size(ms.SFPs, 3)), ' cells to ', num2str(sum(ms.SFPs_good ==1)), ' cells']);
fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPositionMode = 'auto'; fig.PaperPosition = [0 0 5 5]; 
print([ms.dirName '\msCam\dFF_source_extraction\ROI_thr0.6_rm-periphery.png'], '-dpng', '-r300');

% % Plot ROIs overlayed with video
% video = 'D:\Miniscope\Prosocial_ACC\MZ67\20210716\ConcatenatedCrop\msCam\dFF_f1-18000.avi';
% SFPs = ms.SFPs(:, :, ms.SFPs_good == 1);
% msPlotSFPs_video(SFPs, video);

ms.SFPs_0 = ms.SFPs; ms.SFPs = ms.SFPs(:, :, ms.SFPs_good == 1);
ms.RawTraces_0 = ms.RawTraces; ms.RawTraces = ms.RawTraces(:, ms.SFPs_good == 1);
ms.FiltTraces_0 = ms.FiltTraces; ms.FiltTraces = ms.FiltTraces(:, ms.SFPs_good == 1);

save([ms.dirName separator 'ms.mat'],'ms','-v7.3');

% Use CellScreener (by Xinjian Zhang) to remove bad cells
CellScreener;

% if copy_to_googledrive
%     destination_path = char(strcat(copydirpath, separator, ms.Experiment));
%     mkdir(destination_path);
%     copyfile('ms.mat', [destination_path separator 'ms.mat']);
%     copyfile('SFP.mat', [destination_path separator 'SFP.mat']);
%     disp('Successfully copied ms and SFP files to GoogleDrive');
%     try % This is to attempt to copy an existing behav file if you already analyzed it in the past
%             copyfile([ms.dirName separator 'behav.mat'], [destination_path separator 'behav.mat']);
%         catch
%             disp('Behavior not analyzed yet. No files will be copied.');
%     end
% end
% 
% %% 4 - Cleanup temporary files
% rmdir([ms.dirName separator ms.analysis_time], 's');
% 
% %% 5 - Analyse behavior (optional)
% if analyse_behavior
%     behav = msGenerateVideoObj(pwd,'behavCam');
%     behav = msSelectPropsForTracking(behav);
%     trackLength = 95; %cm
%     behav = msExtractBehavoir(behav, trackLength);
%     save([ms.dirName separator 'behav.mat'],'behav','-v7.3');
%     
%     if copy_to_googledrive;
%         destination_path = char(strcat(copydirpath, separator, ms.Experiment));
%         copyfile('behav.mat', [destination_path separator 'behav.mat']);
%         disp('Successfully copied behav file to GoogleDrive');
%     end
% end