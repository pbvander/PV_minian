%% msRunCNMFE_large
% CNMFE source extraction, based on original work by Pengcheng Zhou. This
% version allows you to perform unsupervised analysis. It is recommended to
% first run a few manual analyses to establish your parameters, then use
% the CNMFE_large (patches analysis) with established parameters.
% Based on original script by Pengcheng Zhou, edited by Guillaume Etter

function ms = msRunCNMFE_large(ms, options)

%% Auto-detect operating system
if ispc
    separator = '\'; % For pc operating systems
else
    separator = '/'; % For unix (mac, linux) operating systems
end

%% choose data
neuron = Sources2D();
% nam = get_fullname([ms.dirName separator ms.analysis_time separator 'msvideo.avi']);
nam = options.movie;
nam = neuron.select_data(nam);  %if nam is [], then select data interactively

%% parameters
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = options.pars_envs
% -------------------------      SPATIAL      -------------------------  %
include_residual = options.include_residual; % If true, look for neurons in the residuals
gSig = options.gSig;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
gSiz = options.gSiz;          % pixel, neuron diameter
ssub = options.ssub;           % spatial downsampling factor
with_dendrites = options.with_dendrites;   % with dendrites or not
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = options.updateA_bSiz;
    updateA_dist = neuron.options.dist;
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; %#ok<UNRCH>
    updateA_dist = options.updateA_dist;
    updateA_bSiz = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
spatial_algorithm = 'hals_thresh';

% -------------------------      TEMPORAL     -------------------------  %
Fs = options.Fs;             % frame rate
tsub = options.tsub;           % temporal downsampling factor
deconv_flag = options.deconv_flag; % Perform deconvolution
deconv_options = options.deconv_options;

nk = options.nk;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
% when changed, try some integers smaller than total_frame/(Fs*30)
detrend_method = options.detrend_method;  % compute the local minimum as an estimation of trend.

% -------------------------     BACKGROUND    -------------------------  %
bg_model = options.bg_model;  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = options.nb;             % number of background sources for each patch (only be used in SVD and NMF model)
ring_radius = options.ring_radius;  % when the ring model used, it is the radius of the ring used in the background model.
%otherwise, it's just the width of the overlapping area
num_neighbors = []; % number of neighbors for each neuron

% -------------------------      MERGING      -------------------------  %
show_merge = options.show_merge;  % if true, manually verify the merging step
merge_thr = options.merge_thr;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = options.method_dist;   % method for computing neuron distances {'mean', 'max'}
dmin = options.dmin;       % minimum distances between two neurons. it is used together with merge_thr
dmin_only = options.dmin_only;  % merge neurons if their distances are smaller than dmin_only.
merge_thr_spatial = options.merge_thr_spatial;  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)

% -------------------------  INITIALIZATION   -------------------------  %
K = [];             % maximum number of neurons per patch. when K=[], take as many as possible.
min_corr = options.min_corr;     % minimum local correlation for a seeding pixel, default 0.8
min_pnr = options.min_pnr;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames
save_initialization = true;    % save the initialization procedure as a video. (changed to true)
use_parallel = true;    % use parallel computation for parallel computing
show_init = options.show_init;   % show initialization results %PV changed this so it can be set in msRun (was false)
choose_params = false; % manually choose parameters
center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
% set the value as false when the background fluctuation is small (2p)

% -------------------------  Residual   -------------------------  %
min_corr_res = options.min_corr_res; % Default 0.7
min_pnr_res = options.min_pnr_res;
seed_method_res = 'auto';  % method for initializing neurons from the residual
update_sn = true;

% ----------------------  WITH MANUAL INTERVENTION  --------------------  %
with_manual_intervention = false;

% -------------------------    UPDATE ALL    -------------------------  %
neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
    'gSiz', gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'spatial_constraints', spatial_constraints, ...
    'spatial_algorithm', spatial_algorithm, ...
    'tsub', tsub, ...                       % -------- temporal --------
    'deconv_flag', deconv_flag, 'deconv_options', deconv_options, ...
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background --------
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'num_neighbors', num_neighbors, ...
    'merge_thr', merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', min_corr, ...               % ----- initialization -----
    'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, ...
    'bd', bd, ...
    'center_psf', center_psf);
neuron.Fs = Fs;

%% distribute data and be ready to run source extraction
neuron.getReady(pars_envs);

%% initialize neurons from the video data within a selected temporal range
if choose_params
    % change parameters for optimized initialization
    [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
end

[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel,options.use_prev); %pv added use_prev option
neuron.compactSpatial();
if show_init
    figure();
    ax_init= axes();
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
    drawnow
    
    figure;
    imagesc(PNR);
    drawnow
end

%% estimate the background components
neuron.update_background_parallel(use_parallel);
neuron_init = neuron.copy();

%%  merge neurons and update spatial/temporal components
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, merge_thr_spatial);

%% update spatial components

%% pick neurons from the residual
if include_residual;
    [center_res, Cn_res, PNR_res] = neuron.initComponents_residual_parallel([], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
    if show_init
        figure
        imagesc(Cn_res, [0, 1]); colormap gray; hold on;
        plot(center_res(:, 2), center_res(:, 1), '.g', 'markersize', 10);
        drawnow
        
        figure;
        imagesc(PNR_res);
        drawnow
    end
    neuron_init_res = neuron.copy();
end

%% udpate spatial&temporal components, delete false positives and merge neurons
% update spatial
if update_sn
    neuron.update_spatial_parallel(use_parallel, true);
    udpate_sn = false;
else
    neuron.update_spatial_parallel(use_parallel);
end
% merge neurons based on correlations
neuron.merge_high_corr(show_merge, merge_thr_spatial);

for m=1:2
    % update temporal
    neuron.update_temporal_parallel(use_parallel);
    
    % delete bad neurons
    neuron.remove_false_positives();
    
    % merge neurons based on temporal correlation + distances
    neuron.merge_neurons_dist_corr(show_merge);
end

%% add a manual intervention and run the whole procedure for a second time
neuron.options.spatial_algorithm = 'nnls';
if with_manual_intervention
    show_merge = true;
    neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
    neuron.viewNeurons([], neuron.C_raw);
    
    % merge closeby neurons
    neuron.merge_close_neighbors(true, dmin_only);
    
    % delete neurons
    tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
    ids = find(tags>0);
    if ~isempty(ids)
        neuron.viewNeurons(ids, neuron.C_raw);
    end
end
%% run more iterations
neuron.update_background_parallel(use_parallel);
neuron.update_spatial_parallel(use_parallel);
neuron.update_temporal_parallel(use_parallel);

K = size(neuron.A,2);
tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
neuron.remove_false_positives();
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, merge_thr_spatial);

if K~=size(neuron.A,2)
    neuron.update_spatial_parallel(use_parallel);
    neuron.update_temporal_parallel(use_parallel);
    neuron.remove_false_positives();
end

%% save the workspace for future analysis
neuron.orderROIs('snr');
%cnmfe_path = neuron.save_workspace();

%% show neuron contours
ms.Options = neuron.options;
ms.Centroids = center;
ms.CorrProj = Cn;
ms.PeakToNoiseProj = PNR;

if include_residual;
ms.CentroidsRes = center_res;
ms.CorrProjRes = Cn_res;
ms.PeakToNoiseProjRes = PNR_res;
end

ms.FiltTraces = neuron.C';
ms.RawTraces = neuron.C_raw';
ms.SFPs = neuron.reshape(neuron.A, 2);
ms.numNeurons = size(ms.SFPs,3);

end