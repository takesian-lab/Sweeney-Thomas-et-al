%% recompile_with_markpoints
% Loop through already compiled blocks and run define_markpoints_experiment

%% Paths

recompile = 0; %0 to skip blocks found in the save_blocks_path, 1 to write over
align_type = 'RedoMarkpoints'; %RedoSuite2p

%define_markpoints_experiment options
look_for_previous_MarkpointsExperiments = 0; %1 = look for prev saved data, 0 = redo, 2 = load ensembles from prev saved data but redo everything else
save_fig_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\Activation\5. NDNF to PYR\Markpoints figures'; %'default' or add save path if you want to save all figures to a specific folder
UserOffset = []; %[0 0] to skip offset correction, [] to perform, otherwise input desired x/y offset

%load/save block options
info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\Activation\5. NDNF to PYR';
load_blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\Activation\5. NDNF to PYR\Compiled Blocks';
save_blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\Activation\5. NDNF to PYR\Compiled Blocks XY correct';
info_filename = 'Info_NDNFtoPYR_temp';
stim_protocol = [];
stim_name = '';
        
%% Load Info Sheet

cd(info_path)
Info = tak_read_info_table(info_filename, 'StimProtocol', stim_protocol, 'StimName', stim_name);

%Figure out which blocks have already been realigned
if ~recompile
    cd(save_blocks_path)
    [filenames, ~] = generate_block_filename_from_Info(Info);
    skip = zeros(size(filenames));
    for f = 1:length(filenames)
        if isfile(strcat(filenames(f), '.mat'))
            skip(f) = 1;
            continue
        end
    end
    Info(skip == 1,:) = [];
end

%Return if Info is empty
if isempty(Info)
    disp('No blocks found to realign')
    return;
end

%Make list of remaining blocks
[BlockInfo, ~, ~] = fillSimpleDataTableFromInfo(Info, load_blocks_path, 'StimProtocol', stim_protocol, 'PrintBlockName', 0);

%% Loop through all blocks to recompile

cd(save_blocks_path)
blocknames = BlockInfo.Block;

for b = 1:length(blocknames)
    block_filename = blocknames{b};
    load([load_blocks_path '\' block_filename '.mat']);
    fs = block.setup.framerate;
        
    %Filename created
    disp('Processing...');
    disp(block_filename);
    
    % Redefine parameters based on align_type
    switch align_type

        case 'RedoSuite2p'
            %Redo define_suite2p
            block = define_suite2p_singleblock(block);
            cd(save_blocks_path) %Necessary because define_suite2p moves you to suite2p folder
            
        case 'RedoMarkpoints'
            block = define_markpoints_experiment(block,'LoadPreviousMarkpointsExperiment', look_for_previous_MarkpointsExperiments, 'SaveFigPath', save_fig_path, 'UserOffset', UserOffset);
            save(block_filename, 'block', '-v7.3');
                
        otherwise
            disp('align_type not defined yet')
            return
    end
    
    %Save recompiled block to folder
    save(block_filename, 'block', '-v7.3');
end

disp('All done!')