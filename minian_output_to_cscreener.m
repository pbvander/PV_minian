%% %%%%%%%%----------------README-----------------%%%%%%%%%%%%%%%%%%
%This function was written by Paul Vander to convert Miniscope data analyzed in MiniAn (https://github.com/denisecailab/minian) to a format compatible with CScreener (https://github.com/hongw-lab/CScreener)
%Specifically, this code takes the long, tabular csvs written in the MiniAn notebook, and converts them back to 2/3D arrays (pivoting data to be "wider") in the format the CScreener expects the data to be in
%The final output of this code is the "ms.mat" file that is written at the end - this should be able to be loaded directly into CScreener

%%
function minian_output_to_cscreener (path);
%% Import data and add to ms variable
cd (path)

%A (ms.SFPs)
a = readtable("A.csv"); %Read in data as a table
numNeurons = length(unique(a.unit_id)); %Get the number of unique values in "unit_id" and assign this to numNeurons
vid_h=length(unique(a.height));
vid_w=length(unique(a.width));
a = sortrows(a, ["unit_id" "width" "height"]); %get data sorted in the correct way so that it gets reshaped correctly
A = reshape(a.A, vid_h, vid_w, numNeurons); %reshape data to a 3D array (order of argument here should be opposite the sort step above)
ms = struct('SFPs', A);

%C (ms.FiltTraces)
c = readtable("C.csv");
if length(unique(c.unit_id)) ~= numNeurons %Display warning if number of unit_ids is different
    disp("Warning! Different number of unit_ids")
end
vid_length=length(unique(c.frame));
c = sortrows(c, ["unit_id" "frame"]);
C = reshape(c.C, vid_length, numNeurons);
ms.FiltTraces = C;

%YrA (ms.RawTraces)
yra = readtable("YrA.csv");
if length(unique(yra.unit_id)) ~= numNeurons %Display warning if number of unit_ids is different
    disp("Warning! Different number of unit_ids")
end
if length(unique(yra.frame)) ~= vid_length
    disp("Warning! Different number of frames")
end
yra = sortrows(yra, ["unit_id" "frame"]);
YrA = reshape(yra.YrA, vid_length, numNeurons);
ms.RawTraces = YrA;

%S (ms.S)
s = readtable("S.csv");
if length(unique(s.unit_id)) ~= numNeurons %Display warning if number of unit_ids is different
    disp("Warning! Different number of unit_ids")
end
if length(unique(s.frame)) ~= vid_length
    disp("Warning! Different number of frames")
end
s = sortrows(s, ["unit_id" "frame"]);
S = reshape(s.S, vid_length, numNeurons);
ms.S = S;

%numNeurons
ms.numNeurons = numNeurons;

%cell_label (recommended when using .mat v7.3 or later)
ms.cell_label = ones(numNeurons,1);

%% Save data
save("ms.mat",'ms','-v7.3');