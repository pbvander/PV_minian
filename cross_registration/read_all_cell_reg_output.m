%% %%%%%%%%----------------README-----------------%%%%%%%%%%%%%%%%%%
%This function was written by Paul Vander. It is used to quickly inspect
%and sanity check the alignment parameters across cross-regresistrations,
%then output the mappings as a .csv file for use in R.
% Parts of this code were generated using Claude

%% Things to set manually
exp_dir = 'C:\Users\paulv\Box\correalab\Member Folders\Paul Vander\Experiments';   % Root directory to search
output_dir = 'C:\Users\paulv\Box\correalab\Member Folders\Paul Vander\Data\Torpor project cross-experiment analyses\Miniscope'
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
    
    found = dir(fullfile(searchRoot, '*', 'cellRegistered*.mat')); % One '*' per known level: mouseID / start_date / sessionID
    fprintf('  %s  ->  %d file(s) found\n', dirs{s}, numel(found));
    allFiles = [allFiles; found];
end

%% Read files
rotations = [];
x = [];
y = [];

for i= 1:numel(allFiles)
    name = [allFiles(i).folder '\' allFiles(i).name];
    disp(name)
    S = load(name);
    disp([S.cell_registered_struct.alignment_rotations(2) S.cell_registered_struct.alignment_x_translations(2) S.cell_registered_struct.alignment_y_translations(2)])
    rotations = [rotations; S.cell_registered_struct.alignment_rotations(2)];
    x = [x; S.cell_registered_struct.alignment_x_translations(2)];
    y = [x; S.cell_registered_struct.alignment_y_translations(2)];
    writematrix(S.cell_registered_struct.cell_to_index_map, fullfile(allFiles(i).folder,'cellreg_cell_to_index_map.csv'))
end