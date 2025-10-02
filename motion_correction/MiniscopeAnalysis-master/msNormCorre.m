function ms = msNormCorre(ms,isnonrigid,options);
% Performs fast, rigid registration (option for non-rigid also available).
% Relies on NormCorre (Paninski lab). Rigid registration works fine for
% large lens (1-2mm) GRIN lenses, while non-rigid might work better for
% smaller lenses. Ideally you want to compare both on a small sample before
% choosing one method or the other.
% Original script by Eftychios Pnevmatikakis, edited by Guillaume Etter


warning off all

%% Auto-detect operating system
if ispc
    separator = '\'; % For pc operating systems
else
    separator = '/'; % For unix (mac, linux) operating systems
end

%% Filtering parameters
gSig = 7/ms.ds;
gSiz = 17/ms.ds;
psf = fspecial('gaussian', round(2*gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;
bound = round(ms.height/(2*ms.ds));

template = [];

%Set video output path
disp(options.vid_path) %PV added this to display output path
vid_path = options.vid_path; %PV added this to set custom video location in msRun
writerObj = VideoWriter(vid_path,'Grayscale AVI'); %PV changed this to use vid_path variable
open(writerObj);

%Do motion correction
ms.shifts = [];
ms.meanFrame = [];

for video_i = 1:ms.numFiles;
    name = [ms.vidObj{1, video_i}.Path separator ms.vidObj{1, video_i}.Name];
    disp(['Registration on: ' name]);
    
    % read data and convert to single
    Yf = read_file(name);
    Yf = single(Yf);
    Yf = downsample_data(Yf,'space',1,ms.ds,1);
    
    Y = imfilter(Yf,psf,'symmetric');
    [d1,d2,T] = size(Y);
      
    % Setting registration parameters (rigid vs non-rigid) %PV altered so that most parameters are set in msRun2018
    if isnonrigid
        disp('Non-rigid motion correction...');
    options = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',options.bin_width, 'grid_size',options.grid_size, ... 
    'mot_uf',options.mot_uf,'correct_bidir',options.correct_bidir, ...
    'overlap_pre',options.overlap_pre,'overlap_post',options.overlap_post,'max_shift',options.max_shift);
    else
        disp('Rigid motion correction...');
    options = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);
    end
    
    %% register using the high pass filtered data and apply shifts to original data
    if isempty(template);
        if isnonrigid; %PV added this
            [M1,shifts1,template] = normcorre(Y,options); % register filtered data %PV altered, Original code which "exclude boundaries due to high pass filtering effects": normcorre(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options);
        else
            [M1,shifts1,template] = normcorre(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options); %Original code that PV changed for nonrigid analysis (new code for nonrigid doesn't work for rigid)
        end
    else
        if isnonrigid; %PV added this
            [M1,shifts1,template] = normcorre(Y,options,template); % register filtered data %PV altered, Original code which "exclude boundaries due to high pass filtering effects": normcorre(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options,template);
        else
            [M1,shifts1,template] = normcorre(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options,template); %Original code that PV changed for nonrigid analysis (new code for nonrigid doesn't work for rigid)
        end
    end
    
    Mr = apply_shifts(Yf,shifts1,options,options.overlap_post(1)/2,options.overlap_post(1)/2); % apply shifts to full dataset %PV altered, original code: apply_shifts(Yf,shifts1,options,bound/2,bound/2);
    % apply shifts on the whole movie    
    
    writeVideo(writerObj,uint8(Mr));
    
    %% compute metrics
    [cYf,mYf,vYf] = motion_metrics(Yf,options.max_shift); 
    [cM1f,mM1f,vM1f] = motion_metrics(Mr,options.max_shift);
    
    if video_i == 1;
    ms.meanFrame = mM1f;
    else;
    ms.meanFrame = (ms.meanFrame + mM1f)./2;
    end
    corr_gain = cYf./cM1f*100;
    
    ms.shifts{video_i} = shifts1;
       
end

close(writerObj);

end