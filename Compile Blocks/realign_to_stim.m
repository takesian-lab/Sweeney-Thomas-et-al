%% realign_to_stim
% Loop through already compiled blocks and realign to stim according to align_type

%% Paths

recompile = 0; %Redo previous realignment
align_type = 'Start'; %Start, Target, ReactionTime, RandomLicks, Loco, UserShift
info_path = 'Z:\Carolyn\Fiber Photometry\Behavior\5HT sensor\Info sheets';
load_blocks_path = 'C:\Users\sweenec\Desktop\Local Male 5HT';
save_blocks_path = 'C:\Users\sweenec\Desktop\Local Male 5HT\Realigned no shift use z';
info_filename = 'Info_MalesJune2023';
stim_protocol = [7]; %Leave empty to include all
stim_name = ''; %Leave empty to include all
        
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

%Make list of remaining blocks
[BlockInfo, ~, ~] = fillSimpleDataTableFromInfo(Info, load_blocks_path, 'StimProtocol', stim_protocol, 'PrintBlockName', 0);

%% Loop through all blocks to realign to stim

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

        case 'Start'
            %Trials are aligned to start time
            block.setup.constant.baseline_length = 1; %s
            block.setup.constant.after_stim = 5; %s
            block.setup.constant.locowindow = 5; %s
            sound_time_delay = zeros(size(block.start_time)); %Use start_time since even blocks without functional data will have this

        case 'Target'
            %Trials are aligned to response window 
            block.setup.constant.baseline_length = 1; %s
            block.setup.constant.after_stim = 5; %s
            block.setup.constant.locowindow = 1; %s
            sound_time_delay = block.holdingPeriod + block.waitPeriod;

        case 'ReactionTime'
            %Trials are aligned to response window 
            block.setup.constant.baseline_length = 1; %s
            block.setup.constant.after_stim = 5; %s
            block.setup.constant.locowindow = 1; %s
            reaction_time = block.rxn_time/1000; %convert to s
            reaction_time(reaction_time < 0) = 0; %set trials with no lick to 0s
            sound_time_delay = block.holdingPeriod + block.waitPeriod + reaction_time;
            
        case 'RandomLicks'
            %Trials are aligned to start of a random lick bout from between sound presentations
            block.setup.constant.baseline_length = 1; %s
            block.setup.constant.after_stim = 5; %s
            block.setup.constant.locowindow = 1.5; %s
            [RandomLickTime, LickBoutFound] = align_to_RandomLicks_MaryseBehavior(block);
            sound_time_delay = block.parameters.stimLength + RandomLickTime;
            block.aligned_stim.LickBoutFound = LickBoutFound'; %Store LickBoutFound in aligned_stim
            %figure; imagesc(block.aligned_stim.licks(LickBoutFound,:)); vline(30); %TO VISUALIZE

        case 'Loco'
            %Trials are aligned to start of locomotor bouts
            block.setup.constant.baseline_length = 0.5; %s
            block.setup.constant.after_stim = 5; %s
            block.setup.constant.locowindow = 1.5; %s
            block = align_to_loco(block);
            save(block_filename, 'block', '-v7.3');
            continue; %Skip rest of loop for this type sign align_to_stim happens in aligh_to_loco
            
        case 'Redo'
            %Do nothing, just rerun align_to_stim below
            sound_time_delay = zeros(size(block.start_time));

        case 'UserShift'
            %Trials are aligned to start time
            block.setup.constant.UserShiftTime = -1; %s
            block.setup.constant.baseline_length = 0.5; %s
            block.setup.constant.after_stim = 4; %s
            block.setup.constant.locowindow = 1.5; %s
            sound_time_delay = zeros(size(block.start_time)); %Use start_time since even blocks without functional data will have this
            sound_time_delay = sound_time_delay +  block.setup.constant.UserShiftTime;

        otherwise
            disp('align_type not defined yet')
            return
    end
    
    %Realign to stim
    block = align_to_stim(block, 'SoundTimeDelay', sound_time_delay, 'EstimatedFrameRate', fs);

    %Save realigned block to folder
    save(block_filename, 'block', '-v7.3');
end

disp('All done!')
