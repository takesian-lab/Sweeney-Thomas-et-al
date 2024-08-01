%% compile_blocks_from_info
%
%  This script saves a compiled 'block.mat' file for each row in Info.
%  A compiled block contains Suite2p, Bruker, and Tosca data.
%
%  The Info spreadsheet contains the information required to access the
%  data by specifying the filepaths to each.
%   
%  If data is missing, the script will still run in order to allow for the
%  user to look at components of the data individually.
%
%  Use visualize_block and visualize_cell to preview the block contents once they've been compiled.
%
% VERSION HISTORY
%  - V1 March 2020
%  - V2 December 2020 saves blocks with new naming system
%  - V3 October 2021 adds Fiber Photometry
%  - V4 February 2023, Added user ops files
%% [OpsFile, OpsPath] = uigetfile;

%Show user dialog box to prompt them to select Ops.m file
% see default template for ops files in: ...GitHub\Calcium-Imaging-Analysis\Compile Blocks\CompileBlocks Ops Files    

message = sprintf('Select a CompBlks Ops file, check settings, and hit space bar');
    uiwait(msgbox(message));
    
[CompBlk_OpsFile, CompBlk_OpsPath] = uigetfile;

cd(CompBlk_OpsPath)     %Go to file location
edit(CompBlk_OpsFile)   %Open the file for the user to edit and pause while they do son
pause;                  %Wait for user to hit enter in command line
run(CompBlk_OpsFile)    %Run the Ops.m file once they do so
[MissingOps] = TAK_Lab_Ops_check(CompBlk_OpsFile,CompBlk_OpsPath);

if ~isempty(MissingOps)
    display(MissingOps)
    error('User Ops file does not match template. Missing variables listed in Missing Ops')
end

if ~exist('subcellular_ROIs', 'var')
    subcellular_ROIs = 0;
end

%% Test info and save paths
if ~isfolder(info_path)
    error('Your info path is incorrect.')
end

if ~isfolder(save_path)
    error('Your save path is incorrect.')
end

%% Read Info table
cd(info_path)
Info = tak_read_info_table(info_filename, 'StimProtocol', stim_protocol, 'StimName', stim_name, 'StripFields', false);

%% Load user_ops variable to check Suite2p Fall.ops against user-specified ops

if checkOps
    load(ops_path); %Loads user_ops variable
    user_ops.checkOps = 1;
else
    user_ops.checkOps = 0;
end

%% Prepare to compile blocks for each row of Info

for i = 1:size(Info,1)
    
    %Save all necessary info about the block in setup
    setup = table2struct(Info(i,:));
    setup.constant = constant;
    setup.UpsampleTo30 = upsampleTo30;
    setup.subcellular_ROIs = subcellular_ROIs;
    setup.skip_avi = skip_avi; % new parameter as of 12/1/23. If you are receiving an error, update your ops sheet from template!
    setup.substitute_ToscaTimes_for_VR = substitute_ToscaTimes_for_VR; %new parameter as of march 2024, see template ops
    setup.redo_loco_check = redo_loco_check; %new parameter as of 4/3/24, see template ops

    
    %Update after_stim value if user overwrote with Info sheet
    if isfield(setup, 'after_stim')
        if ~ismissing(setup.after_stim)
            setup.constant.after_stim = setup.after_stim;
        end
    end
    
    %Update baseline_length value if user included and overwrote with Info sheet
    if isfield(setup, 'baseline_length')
        if ~ismissing(setup.baseline_length)
            setup.constant.baseline_length = setup.baseline_length;
        end
    end

    %Generate block filename
    [setup.block_filename, setup.block_supname] = generate_block_filename_from_Info(Info(i,:));
    
    %Filename created
    disp('Processing...');
    disp(setup.block_filename);

    %Skip files that have previously been compiled
    if ~recompile
        cd(save_path)
        if isfile(strcat(setup.block_filename, '.mat'))
            disp('Skipping (already compiled)');
            continue
        end
    end
    
    %% Establish and test paths, allowing for paths to be missing
    
    %TOSCA PATH
    if ismissing(setup.Tosca_session)
        setup.Tosca_path = nan;
        setup.Tosca_session = nan;
        setup.Tosca_run = nan;
    else
        setup.Tosca_path = strcat(setup.pathname, '\', setup.mousename, '\Tosca_', setup.mousename, {'\Session '}, num2str(setup.Tosca_session));
        if ~isfolder(setup.Tosca_path)
            %Accommodate datasets where Tosca path mousename was saved with/without dashes 
            new_mousename = add_or_remove_dashes_from_mousename(setup.mousename);
            setup.Tosca_path = strcat(setup.pathname, '\', setup.mousename, '\Tosca_', new_mousename, {'\Session '}, num2str(setup.Tosca_session));
            
            %If it is still wrong, the user made an error
            if ~isfolder(setup.Tosca_path)
                error('Your Tosca path is incorrect.')
            end
        end
    end
    
    %BLOCK PATH
    if ismissing(setup.block_name)
        setup.block_path = nan;
    else
        if strcmp(setup.analysis_name, 'FibPhot')
            setup.block_path   = strcat(setup.pathname, '\', setup.mousename, '\', setup.expt_date);
        else
            setup.block_path   = strcat(setup.pathname, '\', setup.mousename, '\', setup.expt_date, '\', setup.block_name);
        end
        
        if ~isfolder(setup.block_path)
            error('Your block path is incorrect.')
        end
    end
    
    %VOLTAGE RECORDING (VR) PATH
    if ismissing(setup.VR_name)
        setup.VR_path = nan;
    else
        setup.VR_path  = strcat(setup.pathname, '\', setup.mousename, '\', setup.expt_date, '\', setup.VR_name);
        if ~isfolder(setup.VR_path)
            error('Your voltage recording path is incorrect.')
        end
    end
    
    %SUITE2P PATH
    if ismissing(setup.analysis_name) || strcmp(setup.analysis_name, 'FibPhot')
        setup.suite2p_path = nan;
    elseif isequal(setup.analysis_name, setup.block_name)
        setup.suite2p_path = setup.block_path;
        disp('Using BOTs instead of Suite2p data')
    else
        if isfolder(setup.analysis_name)% added by Maryse and Zahra for when analysis folder is seperate than data
            setup.suite2p_path = setup.analysis_name;
        else
            setup.suite2p_path = strcat(setup.pathname, '\', setup.analysis_name);
            if ~isfolder(setup.suite2p_path)
                error('Your Suite2p analysis path is incorrect.')
            end
        end
    end
    
    %% COMPILE BLOCK
    
    block = struct;
    block.setup = setup;
    
    %pull out the Tosca-derived, behaviorally relevant data (behavior, locomotor, pupillometry)
    [block] = define_behavior_singleblock(block); 
    [block] = define_pupillometry(block);
    if compute_whiskerpad
        [block] = define_whiskerpad(block);
    end
    
    % If running Fiber Photometry, pull out blue/green/red traces and behavioral timestamps
    if strcmp(setup.analysis_name, 'FibPhot')
        [block] = FiberPhotometryAnalysis_blocks(block,Ops.FibPhot);
        if ~any(contains(block.FibPhot.Ops.params.notes,'No fiber photometry data to analyze!!'))
            [block] = define_sound_fibPhot(block);
        end

        [block] = align_to_stim(block, 'EstimatedFrameRate', block.setup.framerate);
        
    elseif ismissing(block.setup.block_path) && ismissing(block.setup.VR_path) % no functional data
        [block] = align_to_stim(block, 'EstimatedFrameRate', EFR);
    else
        %For 2p data, pull out the Bruker-derived timestamps from BOTs and Voltage Recordings
        [block] = define_sound_singleblock(block, WFtrigChange);
        
        %pull out block-specific data from Fall.mat
        [block] = define_suite2p_singleblock(block, 'UserOps', user_ops);
        
        %if running online analysis, this will run instead of define_suite2p_singleblock
        [block] = define_BOTs_for_online_analysis(block);  
      
        %correct markpoints timestamp (for markpoints experiments only)
        [block] = define_markpoints_timestamp(block);
        
        %find the stim-aligned traces for S2P data
        [block] = align_to_stim(block);

        %add metadata for subcellular ROIs (for subcellular imaging experiments only)
        [block] = define_subcellular_ROIs(block);

        %define markpoints experiment details & confirm stimulated cells
        [block] = define_markpoints_experiment(block, 'LoadPreviousMarkpointsExperiment', look_for_previous_MarkpointsExperiments);

    end
    
    %Add custom parameters (v1, v2, v3, stimLength) from info sheet
    [block] = add_custom_parameters_to_block(block);
    
    %% Save block_
    disp('Saving...');
    save(strcat(save_path, '\', setup.block_filename), 'block', '-v7.3');
end

cd(save_path)

disp('Finished compiling all blocks.');
 