% Plot Compiled Correlations
% generate and save summary plots of correlations using Compiled_Table

% 1. load NetworkOps and Compiled_Table
    % example template of NetworkOps is here: 
    % ...GitHub\Calcium-Imaging-Analysis\Network Analysis\Network Analysis Ops Files\Network_Ops_Anne.m
% 2. specify folder name to save
folder_name = 'Summary_AveragedSoundStim'; 

%% save summary figures for pairwise cell correlations

    cd(NetworkOps.save_path);
    [~, ~, ~] = mkdir(NetworkOps.save_path,folder_name);

    for cell = 1:length(NetworkOps.Summary.all_celltypes) % loop through cell types
        NetworkOps.Summary.celltype = NetworkOps.Summary.all_celltypes{cell};
        cd(NetworkOps.save_path);
        [~, ~, ~] = mkdir(folder_name, NetworkOps.Summary.celltype);
        cd(NetworkOps.save_path)
        summary_path = fullfile(NetworkOps.save_path,folder_name,NetworkOps.Summary.celltype);
        PlotAllCorrelations(Compiled_Table,summary_path,NetworkOps.Summary);
    end

    cd(NetworkOps.save_path);
    [~, ~, ~] = mkdir(folder_name, 'AllCells');
    summary_path = fullfile(NetworkOps.save_path, folder_name,'AllCells');
    NetworkOps.Summary.celltype = [];
    PlotAllCorrelations(Compiled_Table,summary_path,NetworkOps.Summary);
