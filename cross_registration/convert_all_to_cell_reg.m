%% %%%%%%%%----------------README-----------------%%%%%%%%%%%%%%%%%%
%This function was written by Paul Vander to convert all ms.mat files for
%each recording session into a CellReg-compatible format using
%cscreener_to_cell_reg.m (see that file for details on the function). Parts
%of this code were generated using Claude

%% Things to set manually
exp_dir = 'C:\Users\paulv\Box\correalab\Member Folders\Paul Vander\Experiments';   % Root directory to search
dirs = ["250417_circulating_E2_torpor_miniscope/pre-OVX_torpor",
        "250417_circulating_E2_torpor_miniscope/post-OVX_torpor",
        "251013_circulating_E2_torpor_miniscope/pre-ovx_torpor",
        "251013_circulating_E2_torpor_miniscope/post-ovx_torpor",
        "260108_circulating_E2_torpor_miniscope/pre-ovx_torpor",
        "260108_circulating_E2_torpor_miniscope/post-ovx_torpor"];

%% Find ms.mat files
allFiles = [];   % will accumulate dir structs across all dirs

for s = 1 : numel(dirs)
    searchRoot = fullfile(exp_dir, dirs{s});

    if ~exist(searchRoot, 'dir')
        fprintf('  [WARN] Subdirectory not found, skipping:\n    %s\n', searchRoot);
        continue
    end
    
    found = dir(fullfile(searchRoot, '*', '*', '*', 'concatenated','minian', 'ms.mat')); % One '*' per known level: mouseID / start_date / sessionID
    fprintf('  %s  ->  %d file(s) found\n', dirs{s}, numel(found));
    allFiles = [allFiles; found];
end

%% Convert files
for i= 1:numel(allFiles)
    name = [allFiles(i).folder '\ms.mat'];
    disp(name)
    cscreener_to_cell_reg(name)
end