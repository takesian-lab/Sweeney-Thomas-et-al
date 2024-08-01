function [TrialData, StimInfo] = makeTrialList(BlockInfo, blocks, save_full_trace)
%Compile all single trial data for list of blocks

%% Optional

if nargin < 3
    save_full_trace = 0; %Option to save full F7, loco, and timestamp
    %This will be done automatically for spontaneous data (stim protocol 12)
end

%% Find Fiber Photometry data

hasFibPhot = any(strcmp('FibPhot', BlockInfo.AnalysisPath));
channels = {'blue', 'green', 'red'};

%% Check frame rates and trial durations before continuing
%If blocks have different frame rates or trial durations (e.g. baseline)
%then we will have a problem because the traces won't match up.

unique_baseline = unique(BlockInfo.Baseline);
unique_after_stim = unique(BlockInfo.After_Stim);
unique_fs = unique(BlockInfo.Framerate);

if length(unique_baseline) > 1 || length(unique_after_stim) > 1
    error('Blocks have different trial durations. Check and recompile blocks to match.')
end

if length(unique_fs) == 1
    %This will be stored in StimInfo
    fs = unique_fs;
    baseline = unique_baseline;
    after = unique_after_stim;
    nBaselineFrames = round(baseline*fs);
    after_inFrames = round(after*fs);
    nFrames = nBaselineFrames + after_inFrames;
else
    error('Blocks have different framerates.')
end

%% Compute variables that may not have been computed during compile blocks

%Loop through each block and add missing variables
fields = fieldnames(blocks);

for b = 1:length(fields)
    block = blocks.(fields{b});
     if strcmp(block.setup.block_path,'cat_block')
         cat_block =1;
     else cat_block = 0;
     end
    
    %There are some old 2P blocks with no loco data
    if ~hasFibPhot && ~isfield(block, 'loco_activity') && ~isfield(block, 'locomotion_trace')
        block.locomotion_trace = block.timestamp; %Use 2P timestamp for loco timestamp
        block.loco_activity = nan(size(block.timestamp)); %Make loco trace with all NaNs
        block.aligned_stim.velocity = nan(size(block.aligned_stim.zscore,2), size(block.aligned_stim.zscore,3));
    end
            
    %Old blocks without functional data won't have active_trials variable
    if ~isfield(block, 'active_trials')
        stopwin = block.setup.constant.locowindow*fs;
        if isfield(block, 'aligned_stim')
            velocity = block.aligned_stim.velocity(:,nBaselineFrames+1:nBaselineFrames+stopwin);
            block.active_trials = any(velocity');
        else
            block.active_trials = nan(size(block.New_sound_times))';
        end
    end
    
    %Old blocks without functional data won't have the Sound_Time variable
    if ~isfield(block, 'Sound_Time')
        block.Sound_Time = get_Sound_Time_behavior(block);
    end
    
    %Z-scored FibPhot trials -> put into same format as 2P data
    if hasFibPhot && isfield(block.aligned_stim, 'F_stim')
        block.aligned_stim.zscore = block.aligned_stim.F_stim_zscore;
    end
    
    %Full trimmed fluoresence and resampled loco, whisker, and pupil trace
    %**Loco will be upsampled, whisker and pupil will be downsampled to match F**
    if save_full_trace || block.setup.stim_protocol == 12 %Always save for spontaneous data
        if isfield(block, 'FibPhot')
            %Stack channels in matrix; if channel doesn't exist placeholder will be NaNs
            F = nan(length(channels),length(block.FibPhot.timestamp_ds));
            for c = 1:length(channels)
                if isfield(block.FibPhot, [channels{c} '_z']) %TODO: decide what to do with dff 
                    F(c,:) = block.FibPhot.([channels{c} '_z']);

                elseif isfield(block.FibPhot, [channels{c} '_dff']) 
                    F(c,:) = block.FibPhot.([channels{c} '_dff']);
                end
            end
            [block.full_timestamp, block.full_F, block.full_loco] = match_fluor_loco_timestamps(F, block.loco_activity, block.FibPhot.timestamp_ds, block.locomotion_trace);
        
        elseif isfield(block, 'F7') %2P
            [block.full_timestamp, block.full_F, block.full_loco] = match_fluor_loco_timestamps(block.F7, block.loco_activity, block.timestamp, block.locomotion_trace,0,0,cat_block);

            %Save suite2p x/y correction and zcorr if they exist, otherwise leave empty
            [block.full_xoff, block.full_yoff, block.full_zcorr] = deal([]);
            if length(block.ops.xoff) == length(block.timestamp)
                [~, block.full_xoff(1,:), ~] = match_fluor_loco_timestamps(double(block.ops.xoff), block.loco_activity, block.timestamp, block.locomotion_trace,0,0,cat_block);
                [~, block.full_yoff(1,:), ~] = match_fluor_loco_timestamps(double(block.ops.yoff), block.loco_activity, block.timestamp, block.locomotion_trace,0,0,cat_block);
            else
                warning('Could not store xoff and yoff. Recompile block to get the latest define_suite2p_singleblock corrections')
            end
            
            %Save zcorr (only in block if user computed zcorr with suite2p)
            if isfield(block,'zcorr')
                [~, block.full_zcorr, ~] = match_fluor_loco_timestamps(block.zcorr, block.loco_activity, block.timestamp, block.locomotion_trace);
            end
        end
        
   
        %Whisker and pupil trace
        if isfield(block, 'full_F')
            if isfield(block, 'concat_whisker'); [~, ~, block.full_whisker] = match_fluor_loco_timestamps(block.full_F, block.concat_whisker, block.full_timestamp, block.pupil_timestamp, 1,0,cat_block); end
            if isfield(block, 'concat_pupil');   [~, ~, block.full_pupil] = match_fluor_loco_timestamps(block.full_F, block.concat_pupil, block.full_timestamp, block.pupil_timestamp, 1,0,cat_block); end
        else
            %If there is no fluoresence data, save behavior traces anyway
            %Upsample loco to match whisker/pupil timestamps (alternatively, we could resample both to match 2P or FibPhot framerate)
            if isfield(block, 'concat_whisker')
                 if isfield(block,'loco_times')
                [block.full_timestamp, block.full_whisker, block.full_loco] = match_fluor_loco_timestamps(block.concat_whisker, block.loco_activity, block.pupil_timestamp, block.loco_times,0,0,cat_block);
                 else
                     [block.full_timestamp, block.full_whisker, block.full_loco] = match_fluor_loco_timestamps(block.concat_whisker, block.loco_activity, block.pupil_timestamp, block.locomotion_trace,0,0,cat_block);
                 end 
                if isfield(block, 'concat_pupil'); [~, ~, block.full_pupil] = match_fluor_loco_timestamps(block.full_whisker, block.concat_pupil, block.full_timestamp, block.pupil_timestamp, 1,0,cat_block); end
            elseif isfield(block,'loco_times')
                block.full_timestamp = block.loco_times;
                block.full_loco = block.loco_activity;
            else
                if isfield(block, ' pupil_timestamp')
                block.full_timestamp = block.pupil_timestamp;
                end
            end
        end
    end
    
    blocks.(fields{b}) = block;
end

%% Make TrialData

TrialData = struct;

for b = 1:length(fields)
    block = blocks.(fields{b});
    
    TrialData(b).Block = block.setup.block_filename;
    
    %STIM DATA
    [v1, v2, v3, stim_length, blank_ind, trialsToIgnore, StimInfo] = simple_prepare_variables(block);
    
    T = ~trialsToIgnore; %An index of trials to keep in analysis (for most protocols this will be all trials)
    
    if any(StimInfo.StimProtocol == [21 211])
        % cant be single
         TrialData(b).Stim_V1 = v1(T);
        TrialData(b).Stim_V2 = v2(T);
        if ~isempty(v3)
            TrialData(b).Stim_V3 = v3(T);
        end
    else
        TrialData(b).Stim_V1 = single(v1(T));
        TrialData(b).Stim_V2 = single(v2(T));
        if ~isempty(v3)
            TrialData(b).Stim_V3 = single(v3(T));
        end
    end
    
    TrialData(b).Stim_Length = single(stim_length(T));
    TrialData(b).IsBlank = logical(blank_ind(T));
    TrialData(b).Sound_Time = single(block.Sound_Time(T)');
    
    %BEHAVIORAL DATA
    TrialData(b).Loco = single(block.aligned_stim.velocity(T,:));
    TrialData(b).IsRunning = single(block.active_trials(T)'); %Convert to integer
    TrialData(b).IsRunning(all(isnan(block.aligned_stim.velocity(T,:))')) = 9; %Use 9 to represent NaN trials where we didn't collect loco data
    if isfield(block.aligned_stim, 'licks'); TrialData(b).Licks = sparse(block.aligned_stim.licks(T,:)); end
    if isfield(block.aligned_stim, 'pupil'); TrialData(b).Pupil = single(block.aligned_stim.pupil(T,:)); end
    if isfield(block.aligned_stim, 'whisker'); TrialData(b).Whisker = single(block.aligned_stim.whisker(T,:)); end
        
    %FUNCTIONAL DATA
    if isfield(block.aligned_stim, 'zscore'); TrialData(b).Trials = single(block.aligned_stim.zscore(:,T,:)); end
    if isfield(block.aligned_stim, 'spks_stim'); TrialData(b).Spikes = single(block.aligned_stim.spks_stim(:,T,:)); end

    %FULL TRACES
    if isfield(block, 'full_timestamp')
        block.full_loco(block.full_loco < block.setup.constant.locoThresh) = 0; %correct for noise floor of wheel (just in case this was not done before now)
        TrialData(b).Timestamp = single(block.full_timestamp);
        TrialData(b).Full_Loco = single(block.full_loco);
        if isfield(block, 'full_F');         TrialData(b).F7 = single(block.full_F); end
        if isfield(block, 'concat_whisker'); TrialData(b).Full_Whisker = single(block.full_whisker); end
        if isfield(block, 'concat_pupil');   TrialData(b).Full_Pupil = single(block.full_pupil); end
        if isfield(block, 'full_xoff');      TrialData(b).Full_xoff = single(block.full_xoff); end
        if isfield(block, 'full_yoff');      TrialData(b).Full_yoff = single(block.full_yoff); end
        if isfield(block, 'full_zcorr');     TrialData(b).Full_zcorr = single(block.full_zcorr); end
        if isfield(block, 'concat_licks');    TrialData(b).Full_licks = single(block.concat_licks);end
    end
    
    %MARKPOINTS
    if isfield(block, 'MarkpointsExperiment')
       TrialData(b).MarkpointsExperiment = block.MarkpointsExperiment.StimTable(T,:);
    end
    
    %PROTOCOL SPECIFIC VARIABLES
    switch block.setup.stim_protocol
        
        case 13 %Maryse behavior stim
            TrialData(b).StimLevel = single(block.stim_level);
            TrialData(b).Outcomes = single(block.Outcome)';
            TrialData(b).ReactionTime = single(block.rxn_time)';
            TrialData(b).HoldingPeriod = single(block.holdingPeriod)';
            TrialData(b).WaitPeriod = single(block.waitPeriod)';
            TrialData(b).WaitCondition = single(block.waitCondition)';
            
            %Remove error trials from results
            block.errorData.result(block.errorData.error_trials) = [];
            block.errorData.result_orig(block.errorData.error_trials) = [];
            TrialData(b).Result = block.errorData.result';
            TrialData(b).Result_orig = block.errorData.result_orig';
        
            %For trials that have been aligned to lick bouts
            if isfield(block.aligned_stim, 'LickBoutFound'); TrialData(b).LickBoutFound = block.aligned_stim.LickBoutFound; end
            
            %Water catch trials
            if isfield(block, 'water_delivery'); TrialData(b).WaterDelivery = single(block.water_delivery)'; end
            
            %Try making V2 Go/No-Go
            StimInfo.Parameters{2} = 'Outcome';
            StimInfo.Units{2} = 'Go/No-Go';
            outcomes = TrialData(b).Outcomes;
            gonogo = nan(size(outcomes));
            gonogo(outcomes == 1) = 1; %Hit
            gonogo(outcomes == 0) = 0; %Miss
            gonogo(outcomes == 3) = 0; %Withhold
            gonogo(outcomes == 4) = 1; %FP
            gonogo(isnan(outcomes)) = 9; %Timeouts
            gonogo(block.water_delivery == 0) = 5; %Catch
            TrialData(b).Stim_V2 = single(gonogo);
            
        case 7 % Frequency discrimination task
            TrialData(b).Outcomes = single(block.Outcome)';
          
            %Water catch trials
            if isfield(block, 'water_delivery'); TrialData(b).WaterDelivery = single(block.water_delivery)'; end
            
            %Try making V2 Go/No-Go
            StimInfo.Parameters{2} = 'Outcome';
            StimInfo.Units{2} = 'Outcome Go/No Go';
            outcomes = TrialData(b).Outcomes;
            gonogo = nan(size(outcomes));
            gonogo(outcomes == 1) = 1; %Hit
            gonogo(outcomes == 0) = 0; %Miss
            gonogo(outcomes == 3) = 3; %Withhold
            gonogo(outcomes == 4) = 4; %FP
            gonogo(block.water_delivery == 0) = 5; %Catch hit
            TrialData(b).Stim_V2 = single(gonogo);
%             TrialData(b).PrepBlock = block.prepTrial;
%             TrialData(b).HitRate = block.HitRate;
    end
end

%% Make list of unique V1, V2 for all blocks together and store in StimInfo

fields = {'Stim_V1', 'Stim_V2'}; %Does not do mapping for V3 yet

concat_Combined = [];
for b = 1:size(TrialData,2)
    blank_ind = TrialData(b).IsBlank;
    temp_Combined = [];
    for f = 1:length(fields)
        if isfield(TrialData, fields{f})
            temp_Combined = [temp_Combined, TrialData(b).(fields{f})];
        end
    end
    temp_Combined(blank_ind,:) = []; %Don't include blank trials as part of variable list
    concat_Combined = [concat_Combined; temp_Combined];
end

%Final StimInfo will reflect data for all blocks combined
StimInfo = add_RF_mapping_to_StimInfo(StimInfo, concat_Combined);
TrialData = add_RF_mapping_to_TrialData(StimInfo, TrialData); %Store StimID for each block

%% Store baseline and frame info

StimInfo.fs = fs;
StimInfo.baseline = baseline;
StimInfo.after = after;
StimInfo.nBaselineFrames = nBaselineFrames;
StimInfo.after_inFrames = after_inFrames;
StimInfo.nFrames = nFrames;
