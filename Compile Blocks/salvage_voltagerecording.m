%% salvage_voltagerecording
% If 2P voltage recording file has been lost or corrupted, this script can attempt to salvage the data
% First, compile all of your blocks with the option substitute_ToscaTimes_for_VR = 1
% That will set the SoundTimes == Tosca sound times, relative to the start of Tosca
% Next, run this script to either manually enter or automatically detect
% the best offset time to try to align SoundTimes to the start of the Bruker recording
% The manual process could be done by trial and error or through a particularly sound-responsive cell
% The automatic process is done by correlating the Loco trace with the Suite2p xoff/yoff, after which
% you may still want to do some additional manual adjustment
% *Needs Tosca data and Suite2p data for this process to work

% WARNING: even once you have done this, your sound times will still not be
% perfectly aligned to Bruker due to drifts in the Tosca PXI clock related to
% the Bruker computer clock. This will be especially bad for longer recordings

% This is meant to be used in times of dire need only, and is not a perfect fix!
% The only way to know if the realignment worked is if you have very convincing sound-responsive cells 

%% Paths

plot_figure = 0;
recompile = 1; %0 to skip blocks found in the save_blocks_path, 1 to write over
align_type = 'Automatic'; %Manual or Automatic
offset = 0; %For Manual only: time in s to offset Tosca data by
info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\NxDL041223M1';
load_blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\NxDL041223M1\Compiled Blocks';
save_blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\NxDL041223M1\Compiled Blocks\Corrected';
info_filename = 'Info_NxDL041223M1_Maryse';
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

        case 'Manual'
            %Realign Tosca to SoundTime + offset
            block.Sound_Time = block.Sound_Time(:,1) + offset;
            [block] = align_Tosca_to_Bruker(block, block.Sound_Time);

        case 'Automatic'
            %Figure out offset
            loco_trace = block.loco_activity;
            loco_time = block.loco_times;
            xoff = double(block.ops.xoff);
            yoff = double(block.ops.yoff);
            XY = sqrt(xoff.^2 + yoff.^2); %proxy for motion in the FOV: compute euclidean distance from x/y shift
            timestamp = block.timestamp;
            
            %Pad the end of loco with zeros so that match_fluor_loco does not trim
            Z = timestamp(end);
            if Z < loco_time(end)
                error('Code assumes loco is always shorter than XY')
            end
            loco_time = loco_time - loco_time(1); %Start at 0
            loco_fr = loco_time(end)/length(loco_time);
            padding = (loco_time(end)+loco_fr):loco_fr:Z;
            if padding(end) < Z
                padding = [padding Z];
            elseif padding(end) > Z
                padding(end) = Z;
            end
            loco_time = [loco_time, padding];
            loco_trace = [loco_trace, zeros(size(padding))];
            
            %Resample loco to match Bruker framerate
            [x, ~, new_loco] = match_fluor_loco_timestamps(XY, loco_trace, timestamp, loco_time, 1, plot_figure);

            %Compute xcorr between traces
            maxlagtime = 60;
            [cc, lag_in_s, MotorCorrData, ~] = tak_compute_xcorr(XY, new_loco, maxlagtime, block.setup.framerate, 'nShuffles', 500, 'pvalue', 0.05, 'PlotFigures', plot_figure);
            offset = MotorCorrData.lag_max;
            
            if (offset+block.loco_times(end)) > timestamp(end)
                
                %Recompute with shorter lag time
                maxlagtime = 20;
                [cc, lag_in_s, MotorCorrData, ~] = tak_compute_xcorr(XY, new_loco, maxlagtime, block.setup.framerate, 'nShuffles', 500, 'pvalue', 0.05, 'PlotFigures', plot_figure);
                offset = MotorCorrData.lag_max;
                
                %If that still didn't work, do manually
                if (offset+block.loco_times(end)) > timestamp(end)
                    error('This offset is impossible. Try adjusting lag time above.')
                end
            end
            
            %Realign Tosca to Sound_Time + offset
            block.Sound_Time = block.Sound_Time(:,1) + offset;
            [block] = align_Tosca_to_Bruker(block, block.Sound_Time);
               
            %Figure
            if plot_figure
                figure; hold on
                
                xmax = max([timestamp(end), loco_time(end)]);
                
                subplot(4,1,1); hold on
                plot(timestamp,XY, 'k')
                xlim([0 xmax])
                title('Suite2p X/Y shift')
                
                subplot(4,1,2);
                plot(loco_time, loco_trace, 'r');
                xlim([0 xmax])
                title('Loco')
                
                subplot(4,1,3);
                plot(block.locomotion_trace, block.loco_activity, 'b')
                xlim([0 xmax])
                title(['Loco corrected +' num2str(offset) ' seconds'])
               
                subplot(4,1,4); hold on
                plot(timestamp,XY, 'k')
                plot(block.locomotion_trace, block.loco_activity, 'b')
                xlim([0 xmax])
                title(['Loco corrected +' num2str(offset) ' seconds'])
                
                sgtitle(block.setup.block_supname)
            end

        otherwise
            disp('align_type not defined yet')
            return
    end
    
    %Redo align to stim
    [block] = align_to_stim(block);

    %Save recompiled block to folder
    save(block_filename, 'block', '-v7.3');
end

disp('All done!')