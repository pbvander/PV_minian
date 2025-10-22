%%Foreword by PV
%This code is meant to combine output videos from
%msRun_mc_only_python_chunk for easy visualization in CScreener

function combine_chunks(vid_dir, vid_name);
%% Auto-detect operating system
%clear
if ispc
    separator = '\'; % For pc operating systems
else
    separator = '/'; % For unix (mac, linux) operating systems
end

%% Get file list
pattern = 'msvideo_[a-z]{1,2}[0-9]{1,4}.avi';
files={dir(vid_dir).name};
matchingFilesIndices = ~cellfun('isempty', regexp(files, pattern, 'once'));
filteredFiles = sort_nat(files(matchingFilesIndices)); %download sort_nat code from https://www.mathworks.com/matlabcentral/fileexchange/10959-sort_nat-natural-order-sort and add the folder to MATLAB path

%% Set up video reading/writing
out_path = [vid_dir separator vid_name];
disp(['Output to ' out_path]);
outputVideo = VideoWriter(out_path,'Grayscale AVI');
initialReader = VideoReader(join([vid_dir, separator, filteredFiles{1}]));
outputVideo.FrameRate = initialReader.FrameRate;

%% Combine videos
open(outputVideo);

for i = 1:length(filteredFiles)
    in_path = [vid_dir separator filteredFiles{i}];
    disp(['Reading ' in_path])
    inputVideo = VideoReader(in_path);
        while hasFrame(inputVideo)
            frame = readFrame(inputVideo);
            writeVideo(outputVideo, frame);
        end
end

close(outputVideo);