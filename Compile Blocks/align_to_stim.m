function [block] = align_to_stim(block, varargin)
% This function pulls out the trial windows for each sound presentation
% and stores stim-aligned neural and behavioral data in block
% 
% Argument(s): 
%   block (struct)
%   varargin: 'SoundTimeDelay'
%   varargin: 'EstimatedFrameRate'
% 
% Returns: 
%   block (struct)
% 
% Notes:
%
% TODO:
% Search 'TODO'
%
% VERSIONS:
% v1 - archived May 13, 2023
% v2 - current version: combines align_to_stim, align_to_stim_FibPhot, and align_to_stim_behavior

%% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% ADDITIONAL PARAMETERS
% List of times (s) to delay each trial start time by (should be same size as Sound_Time)
addOptional(p, 'SoundTimeDelay', []);

% If no functional data, provide desired framerate (WILL NOT AFFECT BLOCKS WITH FUNCTIONAL DATA)
addOptional(p, 'EstimatedFrameRate', 30);

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% Figure out if FibPhot or Suite2p data and get timestamp

multiplaneData = false;
nPlanes = 1;

if ~ismissing(block.setup.VR_name) & isempty(block.setup.analysis_name) % WF but NOT WF online analysis.
    disp('Aligning widefield data to stim...');
    FunctionalDataType = "Widefield";
    
    %For widefield, behavioral data has already been aligned to Bruker timestamps but there is no functional data
    timestamp{1} = block.timestamp;
    
elseif ismissing(block.setup.analysis_name)
    disp('No Fiber Photometry or 2P data found. Aligning behavioral variables to estimated sound times');
    FunctionalDataType = "None";
    
    %Generate variables not contained in block already:
    block.setup.framerate = ops.EstimatedFrameRate;
    block.Sound_Time = get_Sound_Time_behavior(block);
    if isfield(block, 'loco_times')
         block.locomotion_trace = block.loco_times;
    end
   
    if isfield(block, 'concat_pupil') || isfield(block, 'concat_whisker')
        [block] = generate_pupil_timestamp(block); 
    end
    
    %Generate a fake timestamp with desired framerate
    MinDuration = block.Sound_Time(end) + block.setup.constant.after_stim; %Min duration of timestamp needed to include all trials
    if ~isempty(ops.SoundTimeDelay); MinDuration = MinDuration + ops.SoundTimeDelay(end); end %Don't forget to accommodate for sound time delay
    timestamp{1} = 0:(1/ops.EstimatedFrameRate):ceil(MinDuration); %generate timestamp
    
elseif strcmp(block.setup.analysis_name, 'FibPhot')
    disp('Aligning FibPhot data to stim...');
    FunctionalDataType = "FibPhot";
    timestamp{1} = block.FibPhot.timestamp_ds;
else
    disp('Aligning 2P data to stim...');
    FunctionalDataType = "2P";
    
    %Determine if we will need to remove NaNs from MarkPoints data to have pseudo baseline (as close as possible)
    RemoveNansFromMarkPointsData = isfield(block.setup, 'hasMarkPoints') && isfield(block, 'timestamp_idx');
    
    %Accommodate multiplane data
    if isfield(block, 'MultiplaneData')
        multiplaneData = true;
        nPlanes = block.setup.XML.nPlanes;
        timestamp = cell(1,nPlanes);
        for n = 1:nPlanes
            planeName = strcat('plane', num2str(n - 1));
            timestamp{n} = block.timestamp.(planeName);
            if RemoveNansFromMarkPointsData
                idx = block.timestamp_idx.(planeName);
                timestamp{n} = timestamp{n}(idx);
            end
        end
    else
        timestamp{1} = block.timestamp;
        if RemoveNansFromMarkPointsData
            idx = block.timestamp_idx;
            timestamp{1} = timestamp{1}(idx);
        end
    end
end

%% Establish Sound_Time to align to

if block.setup.stim_protocol == 9
    [Sound_Time,block] = align_to_water(block,1); % no sound - this variable is the time that the mouse actually licks the uncued water
    block.Sound_Time = Sound_Time;
elseif block.setup.stim_protocol == 7
    [lick_times,block] = align_to_water(block,2); %
    block.lick_onset_time = lick_times;
end

Sound_Time = block.Sound_Time;

%Add sound_time_delay if needed (used for Maryse behavior stim)
if ~isempty(ops.SoundTimeDelay)
    Sound_Time = Sound_Time + ops.SoundTimeDelay;
end

%Catch errors where Sound_Time is empty
if isempty(Sound_Time)
    warning('No sound times found to align to')
    block.aligned_stim = struct;
    block.active_trials = [];
    return
end

%% FIND EACH TRIAL START AND END

fs = round(block.setup.framerate);
baseline_inFrames = round(block.setup.constant.baseline_length*fs);
after_inFrames = round(block.setup.constant.after_stim*fs);
duration_inFrames = baseline_inFrames + after_inFrames;
locowindow_in_frames = round(block.setup.constant.locowindow*fs);

% loop through each stim-presenation
FrameInd = cell(1,nPlanes); %Store timestamp frame index closest to the start and end of each trial
TrialInd = cell(1,nPlanes); %Store local trial index for the start and end of each trial

for n = 1:nPlanes
    FrameInd{n} = nan(length(Sound_Time),2);
    TrialInd{n} = nan(length(Sound_Time),2);

    for time=1:length(Sound_Time)

        sound = Sound_Time(time);
        if isnan(sound)
            error('No sound time found')
        end
        
        %Align everything to closest frame to sound start
        [~, closest_frame_sound] = min(abs(timestamp{n}(:)-sound));
        A = closest_frame_sound - baseline_inFrames;
        B = closest_frame_sound + after_inFrames - 1;
        a = 1;
        b = duration_inFrames;

        %If user-defined baseline is before the beginning of the block
        %recording, set A = 1 and the beginning of the block will be nan
        if A < 1
            a = abs(A) + 2;
            A = 1;
        end

        %If user-defined trial is longer than block recording, take portion
        %of trial up to the end of recording, the rest of the frames will be nan
        if B > length(timestamp{n})
            B = length(timestamp{n});
            b = length(A:B);
        end

        %Catch problems
        if A > B
            error('A should not be greater than B. Check timestamps')
        end

        FrameInd{n}(time,:) = [A,B];
        TrialInd{n}(time,:) = [a,b];
    end
end

%% ALIGN FIBPHOT DATA

if strcmp(FunctionalDataType, "FibPhot")

    % pre allocate space for 3 channels of FibPhot data
    [F, F_dff] = deal(nan(3,length(timestamp{1})));

    % fill with blue, green, and red data if we have it
    if isfield(block.FibPhot,'blue_F')
        F(1,:) = block.FibPhot.blue_F;
        F_dff(1,:) = block.FibPhot.blue_dff;
    end

    if isfield(block.FibPhot,'green_F')
        F(2,:) = block.FibPhot.green_F;
        F_dff(2,:) = block.FibPhot.green_dff;
    end

    if isfield(block.FibPhot,'red_F')
        F(3,:) = block.FibPhot.red_F;
        F_dff(3,:) = block.FibPhot.red_dff;
    end

    %Preallocate and fill trials
    [F_stim, F_stim_dff] = deal(nan(3,length(Sound_Time),duration_inFrames));
    for t = 1:length(Sound_Time)
        
        a = TrialInd{1}(t,1);
        b = TrialInd{1}(t,2);
        A = FrameInd{1}(t,1);
        B = FrameInd{1}(t,2);
        
        F_stim(:,t,a:b) =  F(:,A:B);
        F_stim_dff(:,t,a:b) =  F_dff(:,A:B);
    end

    Fstim_baseline = F_stim_dff(:,:,1:baseline_inFrames);
    zscore = (F_stim_dff - mean(Fstim_baseline,3,'omitnan'))./std(Fstim_baseline,[],3,'omitnan');

    % SAVE ALIGNED DATA TO BLOCK
    block.aligned_stim.F_stim = F_stim;
    block.aligned_stim.F_stim_dff =  F_stim_dff;
    block.aligned_stim.F_stim_zscore =  zscore;
end

%% ALIGN 2P DATA

if strcmp(FunctionalDataType, "2P")
    
    for n = 1:nPlanes
        if multiplaneData
            planeName = strcat('plane', num2str(n - 1));
            F = block.F.(planeName);
            Fneu = block.Fneu.(planeName);
            F7 = block.F7.(planeName);
            spks = block.spks.(planeName);
            if isfield(block, 'F_chan2')
                F_chan2 = block.F_chan2.(planeName);
                Fneu_chan2 = block.Fneu_chan2.(planeName);
                F7_chan2 = block.F7_chan2.(planeName);
            end
            if isfield(block, 'zcorr'); zcorr = block.zcorr.(planeName); end
            if RemoveNansFromMarkPointsData; idx = block.timestamp_idx.(planeName); end
        else
            F = block.F;
            Fneu = block.Fneu;
            F7 = block.F7;
            spks = block.spks;
            if isfield(block, 'F_chan2')
                F_chan2 = block.F_chan2;
                Fneu_chan2 = block.Fneu_chan2;
                F7_chan2 = block.F7_chan2;
            end
            if isfield(block, 'zcorr'); zcorr = block.zcorr; end
            if RemoveNansFromMarkPointsData; idx = block.timestamp_idx; end
        end

        if RemoveNansFromMarkPointsData
            F = F(:,idx);
            Fneu = Fneu(:,idx);
            F7 = F7(:,idx);
            spks = spks(:,idx);
            if isfield(block, 'F_chan2')
                F_chan2 = F_chan2(:,idx);
                Fneu_chan2 = Fneu_chan2(:,idx);
                F7_chan2 = F7_chan2(:,idx);
            end
        end

        %Skip if no cells exist in this plane
        if isempty(F)
            continue
        end

        %Preallocate trial variables
        nanMat = nan(size(F,1),length(Sound_Time),duration_inFrames);
        [F_stim, Fneu_stim, F7_stim, spks_stim] = deal(nanMat);
        if isfield(block, 'F_chan2'); [F_stim_chan2, Fneu_stim_chan2, F7_stim_chan2] = deal(nanMat); end
        if isfield(block, 'zcorr'); zcorr_stim = nan(size(zcorr,1),length(Sound_Time),duration_inFrames); end

        %Fill trials
        for t = 1:length(Sound_Time)
            a = TrialInd{n}(t,1);
            b = TrialInd{n}(t,2);
            A = FrameInd{n}(t,1);
            B = FrameInd{n}(t,2);

            F_stim(:,t,a:b) =  F(:,A:B);
            Fneu_stim(:,t,a:b) = Fneu(:,A:B);
            F7_stim(:,t,a:b) = F7(:,A:B);
            spks_stim(:,t,a:b) = spks(:,A:B);

            if isfield(block, 'F_chan2')
                F_stim_chan2(:,t,a:b) =  F_chan2(:,A:B);
                Fneu_stim_chan2(:,t,a:b) = Fneu_chan2(:,A:B);
                F7_stim_chan2(:,t,a:b) = F7_chan2(:,A:B);
            end

            if isfield(block, 'zcorr')
                zcorr_stim(:,t,a:b) = zcorr(:,A:B);
            end
        end

        %Compute df_f and zscore based on LOCAL baseline
        F7_baseline = F7_stim(:,:,1:baseline_inFrames);
        df_f = (F7_stim - mean(F7_baseline,3,'omitnan'))./mean(F7_baseline,3,'omitnan');
        zscore = (F7_stim - mean(F7_baseline,3,'omitnan'))./std(F7_baseline,[],3,'omitnan');

        if isfield(block, 'F_chan2')
            F7_baseline_chan2 = F7_stim_chan2(:,:,1:baseline_inFrames);
            df_f_chan2 = (F7_stim_chan2 - mean(F7_baseline_chan2,3,'omitnan'))./mean(F7_baseline_chan2,3,'omitnan');
            zscore_chan2 = (F7_stim_chan2 - mean(F7_baseline_chan2,3,'omitnan'))./std(F7_baseline_chan2,[],3,'omitnan');        
        end

        % SAVE ALIGNED DATA TO BLOCK
        if multiplaneData
             block.aligned_stim.F_stim.(planeName) = F_stim;
             block.aligned_stim.Fneu_stim.(planeName) = Fneu_stim;
             block.aligned_stim.F7_stim.(planeName) = F7_stim;
             block.aligned_stim.spks_stim.(planeName) = spks_stim;
             block.aligned_stim.df_f.(planeName) = df_f;
             block.aligned_stim.zscore.(planeName) = zscore;
             if isfield(block, 'F_chan2')
                 block.aligned_stim.F_stim_chan2.(planeName) = F_stim_chan2;
                 block.aligned_stim.Fneu_stim_chan2.(planeName) = Fneu_stim_chan2;
                 block.aligned_stim.F7_stim_chan2.(planeName) = F7_stim_chan2;
                 block.aligned_stim.df_f_chan2.(planeName) = df_f_chan2;
                 block.aligned_stim.zscore_chan2.(planeName) = zscore_chan2;
             end
             if isfield(block, 'zcorr')
                block.aligned_stim.zcorr_stim.(planeName) = zcorr_stim;
            end
        else
             block.aligned_stim.F_stim = F_stim;
             block.aligned_stim.Fneu_stim = Fneu_stim;
             block.aligned_stim.F7_stim = F7_stim;
             block.aligned_stim.spks_stim = spks_stim;
             block.aligned_stim.df_f = df_f;
             block.aligned_stim.zscore = zscore;
             if isfield(block, 'F_chan2')
                 block.aligned_stim.F_stim_chan2 = F_stim_chan2;
                 block.aligned_stim.Fneu_stim_chan2 = Fneu_stim_chan2;
                 block.aligned_stim.F7_stim_chan2 = F7_stim_chan2;
                 block.aligned_stim.df_f_chan2 = df_f_chan2;
                 block.aligned_stim.zscore_chan2 = zscore_chan2;
             end
             if isfield(block, 'zcorr')
                block.aligned_stim.zcorr_stim = zcorr_stim;
             end
        end
    end
end

%% ALIGN BEHAVIORAL DATA
% Behavioral data are on different timescales/framerates
% To align them, we are finding the closet behavior datapoint to each frame in our FibPhot or Suite2p data

for n = 1:nPlanes
    
    %Preallocate for behavior traces
    nanMat = nan(length(Sound_Time),duration_inFrames);

    behaviorField = {'loco_activity', 'concat_licks', 'concat_pupil', 'concat_whisker'};
        
    for B = 1:length(behaviorField)
        
        %Skip behavior data that does not exist in the block
        if ~isfield(block, behaviorField{B})
            continue;
        end
        
        if B == 1 %LOCO
            behavior_timestamp = block.locomotion_trace;
            behavior_trace = abs(block.loco_activity);
            behavior_trace(behavior_trace < block.setup.constant.locoThresh) = 0; %Special for loco, remove noise floor
            
            %Special for loco, upsample to be on the same timescale as fluoresence data
            if strcmp(FunctionalDataType, "2P")
                [behavior_timestamp, ~, behavior_trace] = match_fluor_loco_timestamps(F7, behavior_trace, timestamp{1}, behavior_timestamp); %F7 and timestamp are just placeholders to get the code to run
            elseif strcmp(FunctionalDataType, "FibPhot")
                [behavior_timestamp, ~, behavior_trace] = match_fluor_loco_timestamps(F, behavior_trace, timestamp{1}, behavior_timestamp); %F and timestamp are just placeholders to get the code to run
            else
                %do nothing
            end

        elseif B == 2 %LICKS
            behavior_timestamp = block.concat_times;
            behavior_trace = block.concat_licks;
            
        elseif B == 3 %PUPIL
            behavior_timestamp = block.pupil_timestamp;
            behavior_trace = block.concat_pupil;
            
        elseif B == 4 %WHISKER
            behavior_timestamp = block.pupil_timestamp;
            behavior_trace =  block.concat_whisker;
        else
            error('Add info for new behavior trace')
        end

        %Catch potential problems
        %You can temporarily comment this section out to compile your block, but it means there's a bigger problem elsewhere we should address
        
        if length(behavior_timestamp) ~= length(behavior_trace)
            %This means your loco traces are not correctly concatenated
            error('Traces do not match. Ask 2P slack for help');
        end 
        
        if any(isnan(behavior_timestamp))
            %Most likely to happen with video data
            error('There are NaNs in behavior timestamp. Ask 2P slack for help');
        end 
        
        %Go through each frame of functional data and find closest behavioral datapoint
        behaviorMat = nanMat;
        for t = 1:length(Sound_Time)
            frames = FrameInd{n}(t,1):FrameInd{n}(t,2);
            positions = TrialInd{n}(t,1):TrialInd{n}(t,2);

            for f = 1:length(frames)
                [~, closest_frame] = min(abs(behavior_timestamp(:) - timestamp{n}(frames(f))));
                behaviorMat(t,positions(f)) = behavior_trace(closest_frame);
            end
        end
        
        if B == 1 %LOCO
            alignedVelocity = behaviorMat;
            
            % Determine if mouse is running or not during each trial
            % Set the sound "window" in which to look for locomotor activity. We will use the sound
            % to end of constant.locowindow. These two numbers are defined at the top of compile_blocks_from_info
            behaviorMat_locowindow = behaviorMat(:,baseline_inFrames+1:baseline_inFrames+locowindow_in_frames);
            active_trials = any(behaviorMat_locowindow,2)'; %Using 'any' works as long as noise floor is removed above
            
        elseif B == 2 %LICKS
            alignedLicks = behaviorMat;
        elseif B == 3 %PUPIL
            alignedPupil = behaviorMat;
        elseif B == 4 %WHISKER
            alignedWhisker = behaviorMat;
        end
    end
    
    % SAVE ALIGNED DATA TO BLOCK
    if multiplaneData
        %Store data for each plane separately because closest behavior frames will be different for each plane
        planeName = strcat('plane', num2str(n - 1));
        if isfield(block, 'loco_activity'); block.active_trials.(planeName) = active_trials; end
        if isfield(block, 'loco_activity'); block.aligned_stim.velocity.(planeName) = alignedVelocity; end
        if isfield(block, 'concat_licks'); block.aligned_stim.licks.(planeName) = alignedLicks; end
        if isfield(block, 'concat_pupil'); block.aligned_stim.pupil.(planeName) = alignedPupil; end
        if isfield(block, 'concat_whisker'); block.aligned_stim.whisker.(planeName) = alignedWhisker; end
    else
        if isfield(block, 'loco_activity'); block.active_trials = active_trials; end
        if isfield(block, 'loco_activity'); block.aligned_stim.velocity = alignedVelocity; end
        if isfield(block, 'concat_licks'); block.aligned_stim.licks = alignedLicks; end
        if isfield(block, 'concat_pupil'); block.aligned_stim.pupil = alignedPupil; end
        if isfield(block, 'concat_whisker'); block.aligned_stim.whisker = alignedWhisker; end
    end
end

end %end function
