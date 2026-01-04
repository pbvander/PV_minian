%%Foreword by PV
%This code is meant to filter bad frames from output of msNormCorre for use in CScreener (CScreener does not align Minian data based on frame index so
%you must remove bad frames from video or there will be mis-match).
%Specifically, it creates a copy of the input video, but with bad frames removed (skipped over when writing output video)

function filter_bad_frames(in_path, out_path, bad_frames);
%% Set up video reading/writing
outputVideo = VideoWriter(out_path,'Grayscale AVI');
initialReader = VideoReader(in_path);
outputVideo.FrameRate = initialReader.FrameRate;

%% Combine videos
open(outputVideo);

inputVideo = VideoReader(in_path);
frame_num=1;
disp(join(["First bad frame:", bad_frames(1,1)]))
disp(join(["Last bad frame:", bad_frames(1,end)]))
while hasFrame(inputVideo)
    frame = readFrame(inputVideo);
    if ~ismember(frame_num, bad_frames)
        writeVideo(outputVideo, frame);
    else
        if rem(frame_num,100) == 0 || frame_num == bad_frames(1,1)
            disp(join(["Skipping", frame_num]))
        end
    end
    if rem(frame_num,5000)==0
        disp(join([frame_num, "frames completed"]))
    end
    frame_num=frame_num+1;
end

close(outputVideo);