%% %%%%%%%%----------------README-----------------%%%%%%%%%%%%%%%%%%
%This function was written by Paul Vander to convert Miniscope data
%previous converted by minian_output_to_cscreener.m to a format compatible
%with CellReg (https://github.com/zivlab/CellReg)

%%
function cscreener (path);
ms=load(path);
A = permute(ms.ms.SFPs, [3,1,2]);
[dir_path, ~,~] = fileparts(path);
save(fullfile(dir_path,"cellreg.mat"),"A");
