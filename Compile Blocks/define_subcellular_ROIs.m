function [block] = define_subcellular_ROIs(block)
% If subcellular ROIs are delineated, pull out info about each ROI
% 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes:
% Created on 1/2/2024 to incorporate specific information about subcellular
% ROIs. 


%% Skip if no Subcellular Data

if block.setup.subcellular_ROIs == 0
    return
end

%% Setup

disp('Adding Subcellular data to block...');

%Accomodate multiplane data
if isfield(block, 'MultiplaneData')
    multiplaneData = true;
    nPlanes = block.setup.XML.nPlanes;
else
    multiplaneData = false;
    nPlanes = 1;
end

%Accommodate channel 2
if isfield(block, 'F_chan2')
    chan2_exists = true;
else
    chan2_exists = false;
end

%Single or multi-BOT file?
if ~isequal(block.setup.BOT_filename, 'NaN')
    singleBOT = true;
else
    singleBOT = false;
end

%% Add Subcellular Info to Block

if contains(block.setup.analysis_name, 'plane') % if analyzing each plane separately
    foldername = extractBefore(block.setup.analysis_name,'\plane');
    cd (foldername);
    load('snt_paths_data.mat'); % load path data
    block.subcellular.sntPaths = snt_paths_data; % save snt path data in block
    load('snt_suite2p_conversion_global_info.mat'); % snt load global info
    block.subcellular.globalInfo = snt_suite2p_conversion_global_info; % save snt global path in block
    cd(block.setup.analysis_name);
    plane_snt_data = load('snt_suite2p_conversion_data.mat');
    block.subcellular.planeData = plane_snt_data;
    cd(foldername);

else % if analyzing two planes together
    cd(block.setup.analysis_name);
    load('snt_paths_data.mat'); % load path data
    block.subcellular.sntPaths = snt_paths_data; % save snt path data in block
    load('snt_suite2p_conversion_global_info.mat'); % snt load global info
    block.subcellular.globalInfo = snt_suite2p_conversion_global_info; % save snt global path in block

    for n = 1:nPlanes
        planeName = strcat('plane', num2str(n - 1));
        cd(planeName);
        plane_snt_data = load('snt_suite2p_conversion_data.mat');
        block.subcellular.(planeName) = plane_snt_data;
        cd(block.setup.analysis_name);
    end
end
end