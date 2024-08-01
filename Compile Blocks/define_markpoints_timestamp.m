function [block] = define_markpoints_timestamp(block)
% Estimate the timing of markpoints stimulation with respect to Bruker timestamp
% If we have voltage recording and Tosca trials as the markpoints trigger, this will be accurately recorded
% If not, we will make the best estimate that we can using other information
% During markpoints stimulation, imaging stops to protect the PMTs:
% Where data does not exist because the PMTs were closed, adjust the
% timestamp to add NaNs and correct suite2p data if present
% 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes:
%
% VERSION HISTORY:
% - v1: archived May 16, 2023
% - v2: current version, updates corresponding to XML improvement of define_sound_singleblock

%% Skip if no MarkPoints data

if ~isfield(block.setup, 'hasMarkPoints')
    return
end

%% Setup

disp('Correcting MarkPoints Timestamp...');

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

%% Create NaN-adjusted timestamp

for n = 1:nPlanes
    
    if multiplaneData
        planeName = strcat('plane', num2str(n - 1));
        timestamp = block.timestamp.(planeName); %BOT timestamp
    else
        timestamp = block.timestamp; %BOT timestamp
    end
    
    %Sometimes, no TIFFs are recorded with markpoints, so they are not represented in the timestamp
    %If this is the case, find number of frames that would have been recorded if imaging continued during markpoints
    frameDiff = diff(timestamp); %time (s) between each frame
    frameMode = mode(frameDiff); %most common frame period
    nPeriodsPerFrame = [1; round(frameDiff/frameMode)]; %add 1 at the beginning because you lose a frame with diff
    timestamp_idx = cumsum(nPeriodsPerFrame);

    %create new timestamp where missing periods will be filled by NaNs
    new_timestamp = nan(sum(nPeriodsPerFrame),1);
    new_timestamp(timestamp_idx) = timestamp;

    %If imaging did take place continuosly, this fix will not affect timestamp
    %new_timestamp == timestamp; %FYI
    
    %Save new timestamp
    %block.timestamp(block.timestamp_idx) will get you back the original timestamp
    if multiplaneData
        block.timestamp.(planeName) = new_timestamp;
        block.timestamp_idx.(planeName) = timestamp_idx;
    else
        block.timestamp = new_timestamp;
        block.timestamp_idx = timestamp_idx;
    end
end

%% If Suite2p data is available, NaN-adjust that as well

if ~ismissing(block.setup.suite2p_path)
    
    for n = 1:nPlanes
        if multiplaneData
            planeName = strcat('plane', num2str(n - 1));
            idx = block.timestamp_idx.(planeName);
            nCells = length(block.cell_number.(planeName));
            nanmat = nan(nCells, length(block.timestamp.(planeName)));
            
            block.F.(planeName) = nanAdjust(block.F.(planeName), idx, nanmat);
            block.Fneu.(planeName) = nanAdjust(block.Fneu.(planeName), idx, nanmat);
            block.F7.(planeName) = nanAdjust(block.F7.(planeName), idx, nanmat);
            block.df_f.(planeName) = nanAdjust(block.df_f.(planeName), idx, nanmat);
            block.spks.(planeName) = nanAdjust(block.spks.(planeName), idx, nanmat);
            block.ops.xoff.(planeName) = nanAdjust(block.ops.xoff.(planeName), idx, nanmat(1,:));
            block.ops.yoff.(planeName) = nanAdjust(block.ops.yoff.(planeName), idx, nanmat(1,:));
            
            if isfield(block, 'zcorr')
                zcorr_nanmat = nan(size(block.zcorr.(planeName), 1), length(block.timestamp.(planeName)));
                block.zcorr.(planeName) = nanAdjust(block.zcorr.(planeName), idx, zcorr_nanmat);
            end
            
            if chan2_exists
                block.F_chan2.(planeName) = nanAdjust(block.F_chan2.(planeName), idx, nanmat);
                block.Fneu_chan2.(planeName) = nanAdjust(block.Fneu_chan2.(planeName), idx, nanmat);
                block.F7_chan2.(planeName) = nanAdjust(block.F7_chan2.(planeName), idx, nanmat);
                block.df_f_chan2.(planeName) = nanAdjust(block.df_f_chan2.(planeName), idx, nanmat);
            end
        else
            idx = block.timestamp_idx;
            nCells = length(block.cell_number);
            nanmat = nan(nCells, length(block.timestamp));
           
            block.F = nanAdjust(block.F, idx, nanmat);
            block.Fneu = nanAdjust(block.Fneu, idx, nanmat);
            block.F7 = nanAdjust(block.F7, idx, nanmat);
            block.df_f = nanAdjust(block.df_f, idx, nanmat);
            block.spks = nanAdjust(block.spks, idx, nanmat);
            block.ops.xoff = nanAdjust(block.ops.xoff, idx, nanmat(1,:));
            block.ops.yoff = nanAdjust(block.ops.yoff, idx, nanmat(1,:));
            
            if isfield(block, 'zcorr')
                zcorr_nanmat = nan(size(block.zcorr, 1), length(block.timestamp));
                block.zcorr = nanAdjust(block.zcorr, idx, zcorr_nanmat);
            end
            
            if chan2_exists
                block.F_chan2 = nanAdjust(block.F_chan2, idx, nanmat);
                block.Fneu_chan2 = nanAdjust(block.Fneu_chan2, idx, nanmat);
                block.F7_chan2 = nanAdjust(block.F7_chan2, idx, nanmat);
                block.df_f_chan2 = nanAdjust(block.df_f_chan2, idx, nanmat);
            end
        end   
    end
end
    
%% Estimate MarkPoints times

if length(block.setup.XML.markPoints.filenames) > 1
    %T-SERIES MARKPOINTS
    %T-Series markpoints have a separate markpoints file for each activation

    %Markpoint sequences do not come with an absolute or relative time in XML (unlike BOTs)
    %They have a save time, which represents when the sequence was reached in T-series and saved as MarkPoints.xml
    %This is most often NOT the time of activation (especially if it was waiting for a sound trigger)

    XML = block.setup.XML;
    
    if multiplaneData
        timestamp = block.timestamp.combined;
    elseif ~multiplaneData && XML.nPlanes > 1
        %Situation where we have multiplane data that was run through Suite2p as single plane
        %Figure out which index the current plane corresponds to
        timestamp_difference = nan(1,XML.nPlanes);
        for n = 1:XML.nPlanes
            currentTime = XML.relativeTime(XML.index == n);
            minLength = min([length(currentTime), length(timestamp)]); %If the planes have a different number of frames this is probably already enough info for us to know which is which, but in some cases they might not
            timestamp_difference(n) = sum(timestamp(1:minLength) - currentTime(1:minLength));
        end
        [~, index] = min(abs(timestamp_difference));

        %Adjust sequenceType, parameterSet, and savetime to have the correct number of frames to match this plane
        startSequenceInd = length(XML.sequenceType) - sum(XML.index == index) + 1;
        XML.sequenceType = XML.sequenceType(startSequenceInd:end);
        XML.parameterSet = XML.parameterSet(startSequenceInd:end);
        XML.savetime = XML.savetime(startSequenceInd:end);
        
        %Now update sequence, absoluteTime, relativeTime, and index knowing the correct plane
        XML.sequence = XML.sequence(XML.index == index);
        XML.absoluteTime = XML.absoluteTime(XML.index == index);
        XML.relativeTime = XML.relativeTime(XML.index == index);
        XML.index = XML.index(XML.index == index);
    end

    sequence = XML.sequence; %tells us which sequence each BOT frame corresponds to
    sequenceType = XML.sequenceType; %FYI: A sequence = a BOT, markpoints, or other segment of a T-series script

    %convert save times to seconds
    savetime = XML.savetime; %save times (one per sequence)
    savetimevec = datevec(savetime);
    savetimeInSeconds = savetimevec(:,4)*3600 + savetimevec(:,5)*60 + savetimevec(:,6);
    savetimeInSeconds = savetimeInSeconds - savetimeInSeconds(1);

    %adjust save times with respect to first sequence
    %SOMETIMES the first sequence is not a BOT (it could be a markpoints)
    %If it is not a BOT, adjust all times to be relevant to the first BOT
    firstBOTSequence = sequence(1); %Sequence corresponding to first BOT frame
    adjustedSaveTime = savetimeInSeconds - savetimeInSeconds(firstBOTSequence); %Set time relative to first BOT frame
    adjustedSaveTime = adjustedSaveTime + timestamp(1); %Add any timestamp offset

    %Use first BOT timestamps as actual times
    BOT_indices = [0; find(diff(sequence))] + 1; %add 1 at the beginning because you lose a frame with diff
    BOT_time = timestamp(BOT_indices);

    %Estimate markpoints times based on BOT times and markpoints duration
    time = nan(size(adjustedSaveTime));
    BOT_ind = find(strcmp(sequenceType, 'BOT') + strcmp(sequenceType, 'Combined'));
    time(BOT_ind) = BOT_time;
    block.setup.XML.time = time; %Save in XML.markPoints
    
    markpoints_seq = find(strcmp(sequenceType, 'MarkPoints'));
    if ~isempty(markpoints_seq)
        for m = 1:length(markpoints_seq)
            ind = markpoints_seq(m);
            markpointsDuration = (XML.markPoints.Duration(m)/1000); %convert to ms
            if ind == length(time)
                time(ind) = time(ind-1) + frameMode; %If there is no BOT after an MP sequence, we cannot know the actual time, so just add one frame
            else
                time(ind) = time(ind+1) - markpointsDuration;
            end
        end
        block.setup.XML.markPoints.time = time(strcmp(sequenceType, 'MarkPoints'));
    end

    % Sanity check figure
    plot_figure = 0;

    if plot_figure
        figure; hold on
        subplot(3,1,1); hold on
        title('Timestamp (Green = BOT, Red = MarkPoints)')
        plot(timestamp, rand(size(timestamp)), 'Linewidth', 0.5)
        ylim([-0.5 1.5])
        vline(time(strcmp(sequenceType, 'BOT'))', 'g')
        vline(adjustedSaveTime(strcmp(sequenceType, 'MarkPoints'))', 'r')
        legend({'Fake trace'})

        subplot(3,1,2); hold on
        title('NaN-Adjusted Timestamp (Green = BOT, Red = MarkPoints)')
        plot(new_timestamp, rand(size(new_timestamp)), 'Linewidth', 0.5)
        ylim([-0.5 1.5])
        vline(time(strcmp(sequenceType, 'BOT'))', 'g')
        vline(adjustedSaveTime(strcmp(sequenceType, 'MarkPoints'))', 'r')
        xlabel('Time (s)')
        legend({'Fake trace'})
        
        subplot(3,1,3); hold on
        title('NaN-Adjusted Timestamp with Estimated MarkPoints (Green = BOT, Red = MarkPoints)')
        plot(new_timestamp, rand(size(new_timestamp)), 'Linewidth', 0.5)
        ylim([-0.5 1.5])
        vline(BOT_time', 'c')
        vline(time(strcmp(sequenceType, 'BOT'))', 'g')
        vline(time(strcmp(sequenceType, 'MarkPoints'))', 'r')
        xlabel('Time (s)')
        legend({'Fake trace'})
    end
    
    %Find whether uncaging shutter was open or closed using parameter set
    parameterSet = block.setup.XML.parameterSet;
    if ~isempty(markpoints_seq) %Trial based markpoints
        uncagingShutter = nan(length(markpoints_seq),1);
        for i = 1:length(markpoints_seq)
            %Find parameter set associated with BOT right after each markpoints
            currentSequenceIdx = markpoints_seq(i)+1;
            if isempty(currentSequenceIdx)
                %If markpoints is the last sequence of the experiment, use previous value instead
                currentSequenceIdx = markpoints_seq(i)-1;
            end

            currentParameterSet = parameterSet(currentSequenceIdx);
            if isempty(currentParameterSet)
                error('Could not find ParameterSet. Ask 2P slack for help.') 
            elseif strcmp(currentParameterSet, 'CurrentSettings') %CurrentSettings is the default
                uncagingShutter(i) = 0; %Closed
            else
                %The parameter set for the uncaging shutter will have
                %different names based on what the user called it. Ex. Parameter Set 1
                uncagingShutter(i) = 1; %Open
            end
        end
        uncagingShutter = uncagingShutter(~isnan(uncagingShutter))';
        trialsToIgnore = zeros(size(uncagingShutter));
    else %Continuous BOT markpoints
        nMarkPoints = length(block.setup.XML.markPoints.Iteration);
        uncagingShutter = nan(nMarkPoints,1)'; 
        uncagingShutter(:) = block.setup.XML.uncagingShutter;
        trialsToIgnore = zeros(size(uncagingShutter));
    end
else
    %BOT MARKPOINTS
    %BOT markpoints have one markpoints file for all activations
    
    time = block.setup.XML.markPoints.time;
    nMarkPoints = length(time);
    uncagingShutter = nan(nMarkPoints,1)'; 
    uncagingShutter(:) = block.setup.XML.uncagingShutter;
    trialsToIgnore = zeros(size(uncagingShutter));
   
    %Figure to check timestamps before Tosca correction
    plot_figure = 0;

    if plot_figure            
        figure; hold on
        subplot(2,1,1); hold on
        title('Timestamp')
        plot(timestamp, rand(size(timestamp)), 'Linewidth', 0.5)
        vline(time, 'g')

        subplot(2,1,2); hold on
        title('NaN-Adjusted Timestamp')
        plot(new_timestamp, rand(size(new_timestamp)), 'Linewidth', 0.5)
        vline(time, 'g')
        xlabel('Time (s)')
    end
end

%% Estimate Sound Times
% This portion is meant to mimic the part of define_sound_singleblock where
% we put Tosca data and Bruker data on the same timescale
% In this case we may or may not have a voltage recording, but we do know the order
% of sound presentation and we know that MarkPoints was triggered by
% sounds, so we can use the Markpoints times as an estimate for sound time

nanadjustBOT = 0;

if ismissing(block.setup.Tosca_path) || any(ismissing(block.setup.VR_filename)) %No sound or loco information, so replace sound info with markpoints info
    block.Sound_Time = block.setup.XML.markPoints.time;
    block.parameters.variable1 = 0; %Like spontaneous data
    block.parameters.variable2 = 0; %Like spontaneous data
    block.parameters.stimLength = block.setup.XML.markPoints.Duration;
    block.parameters.uncagingShutter = uncagingShutter;
    block.parameters.trialsToIgnore = trialsToIgnore;
else

    %If there are VR files, use them to get the actual sound time, otherwise bruker_trial_time estimated above will be used
    VR_files = block.setup.VR_filename;
    disp(['Loading ' num2str(length(VR_files)) ' VoltageRecording files'])

    if length(VR_files) == 1 %One VR file with all triggers

        %Prepare to nan-adjust Suite2p data later if needed
        if ~ismissing(block.setup.suite2p_path) && length(unique(diff(timestamp_idx))) == 1
            nanadjustBOT = 1;
        end

        %Load VR file
        M = csvread(strcat(block.setup.block_path, '\', VR_files),1,0);

        %Triggers are 5V signals that are active for 19-20ms
        Bruker_trial_triggers = M(M(:,4) > 3)./1000; %Detect all timepoints in the VR that are > 3V and convert to seconds
        diffTrigger = diff([0; Bruker_trial_triggers]); %append 0 to the start so that we won't be off by 1 when using diff function
        Bruker_trial_time = Bruker_trial_triggers(diffTrigger > 1);

        %FYI: Error trials will be removed from Bruker_trial_time in align_Tosca_to_Bruker

    else %One VR file per trigger
        tosca_trigger_time = nan(length(VR_files),1);
        for v = 1:length(VR_files)
            M = csvread(strcat(block.setup.block_path, '\', VR_files{v}),1,0);

            % find the start of each trial, and align times to it
            first_trigger = M(find(M(:,2) > 3, 1, 'first') )./1000;
            if ~isempty(first_trigger)
                tosca_trigger_time(v) = first_trigger;
            end
        end

        %Remove error trials
        if ~isempty(block.errorData.error_trials)
            error('First time having an error trial with this type of data. Add code')
        end

        %Because we cannot record any triggers DURING activation in the trial based T-Series,
        %Bruker_trial_triggers is the time when the VR file recorded the
        %first trigger AFTER the start of the trial. This means that if we
        %subtract the duration between trial start to trigger 1 we will
        %know more accurately when the sound played
        sound_state_duration = block.errorData.target_state_duration';
        actual_sound_time = tosca_trigger_time - sound_state_duration;
        if ~all(block.New_sound_times == 0); error('If these do not equal 0 this may be a new paradigm'); end

        %Correct markpoints times using actual_sound_time + MarkPoints duration
        markpoints_seq = find(strcmp(sequenceType, 'MarkPoints'));
        for m = 1:length(markpoints_seq)
            ind = markpoints_seq(m);
            markpointsDuration = (XML.markPoints.Duration(m)/1000); %convert to s
            block.setup.XML.time(ind) = block.setup.XML.time(ind+1) + markpointsDuration + actual_sound_time(m);
        end
        Bruker_trial_time = block.setup.XML.time(markpoints_seq); %TODO: Where to remove error trials?
        block.setup.XML.markPoints.time = Bruker_trial_time; %TODO: Where to remove error trials?

        %FOR SOUND COMBINED WITH MARKPOINTS T-SERIES ONLY:
        %When there was sound stimulation, we marked the shutter opening with a Tosca trigger
        %The only stimulus happening during this trial should be the sound of the shutter itself
        %This should be the first trial every time the shutter status changes
        if any(isnan(uncagingShutter)) %Uncaging shutter wasn't recorded
            trialsToIgnore = zeros(size(uncagingShutter));
        else %Uncaging shutter opens and closes
            trialsToIgnore = abs([1, diff(uncagingShutter)]);
        end
    end

    block.Sound_Time = Bruker_trial_time;
end

% Align Tosca variables (Loco, Licks) to Bruker timestamp as in define_sound_singleblock
if ~ismissing(block.setup.Tosca_path)
    [block] = align_Tosca_to_Bruker(block, block.Sound_Time); %block.Sound_Time will be updated in this function
    block.setup.XML.markPoints.time = block.Sound_Time; 
    block.parameters.uncagingShutter = uncagingShutter;
    block.parameters.trialsToIgnore = trialsToIgnore;
end


%% NAN-adjust in the opposite direction for data recorded using BOT
% For BOTs, data is continuously recorded but PMT closes during stim
% Replace this data with NANs to treat data the same as T-Series

if nanadjustBOT
    trialDuration = block.setup.XML.markPoints.Duration/1000;
    for n = 1:nPlanes
        
        %Find closest frames to each Sound_Time and compile list of frames where we expect PMTs to be zeroed
        frameStart = nan(size(block.Sound_Time));
        frameEnd = nan(size(block.Sound_Time));
        zeroInd = [];
            
        if multiplaneData
            planeName = strcat('plane', num2str(n - 1));

            for t = 1:length(block.Sound_Time)
                [~, frameStart(t)] = min(abs(block.timestamp.(planeName) - block.Sound_Time(t)));
                [~, frameEnd(t)] = min(abs(block.timestamp.(planeName) - (block.Sound_Time(t) + trialDuration(t))));
                zeroInd = [zeroInd, frameStart(t):frameEnd(t)];
            end

            %Replace zeroInd values with NaN
            block.timestamp.(planeName)(zeroInd) = NaN;
            block.timestamp_idx.(planeName)(zeroInd) = [];
            block.F.(planeName)(:,zeroInd) = NaN;
            block.Fneu.(planeName)(:,zeroInd) = NaN;
            block.F7.(planeName)(:,zeroInd) = NaN;
            block.df_f.(planeName)(:,zeroInd) = NaN;
            block.spks.(planeName)(:,zeroInd) = NaN;
            block.ops.xoff.(planeName)(1,zeroInd) = NaN;
            block.ops.yoff.(planeName)(1,zeroInd) = NaN;

            if isfield(block, 'zcorr')
                block.zcorr.(planeName)(:,zeroInd) = NaN;
            end

            if chan2_exists
                block.F_chan2.(planeName)(:,zeroInd) = NaN;
                block.Fneu_chan2.(planeName)(:,zeroInd) = NaN;
                block.F7_chan2.(planeName)(:,zeroInd) = NaN;
                block.df_f_chan2.(planeName)(:,zeroInd) = NaN;
            end
        else

            for t = 1:length(block.Sound_Time)
                nBlankFrames = ceil(trialDuration(t)/(1/block.setup.framerate))+1;
                [~, frameStart(t)] = min(abs(block.timestamp - block.Sound_Time(t)));
                [~, frameEnd(t)] = min(abs(block.timestamp - (block.Sound_Time(t) + trialDuration(t))));
                %zeroInd = [zeroInd, frameStart(t):frameEnd(t)];
                zeroInd = [zeroInd, (frameStart(t)):(frameStart(t)+nBlankFrames)]; %This is a more conservative blank frame estimate than what we are doing for multiplane data, because we have higher framerate here
            end

            %Replace zeroInd values with NaN
            block.timestamp(zeroInd) = NaN;
            block.timestamp_idx(zeroInd) = [];
            block.F(:,zeroInd) = NaN;
            block.Fneu(:,zeroInd) = NaN;
            block.F7(:,zeroInd) = NaN;
            block.df_f(:,zeroInd) = NaN;
            block.spks(:,zeroInd) = NaN;
            block.ops.xoff(1,zeroInd) = NaN;
            block.ops.yoff(1,zeroInd) = NaN;

            if isfield(block, 'zcorr')
                block.zcorr(:,zeroInd) = NaN;
            end

            if chan2_exists
                block.F_chan2(:,zeroInd) = NaN;
                block.Fneu_chan2(:,zeroInd) = NaN;
                block.F7_chan2(:,zeroInd) = NaN;
                block.df_f_chan2(:,zeroInd) = NaN;
            end
        end   
    end
    
    % Sanity check figure
    plot_figure = 0;

    if plot_figure
        figure; hold on
        subplot(2,1,1); hold on
        title('Timestamp')
        plot(timestamp, rand(size(timestamp)), 'Linewidth', 0.5)
        vline(block.Sound_Time, 'g')

        subplot(2,1,2); hold on
        title('NaN-Adjusted Timestamp')
        plot(block.timestamp, rand(size(block.timestamp)), 'Linewidth', 0.5)
        vline(block.Sound_Time, 'g')
        xlabel('Time (s)')
    end
end

%% Upsample 2p data to 30Hz if needed
%This has to be done after NaN-correcting for BOT files otherwise 0s from PMT
%shutter closing get upsampled with adjacent data

if block.setup.UpsampleTo30 && block.setup.framerate < 30
    %Remove timestamp_idx from block since upsampling will eliminate nans
    block = rmfield(block, 'timestamp_idx');
    
    disp('Upsampling data to 30Hz...')
    block = upsample_suite2p_data(block,30);
end

end %function end

function B = nanAdjust(A, idx, nanmat)
    B = nanmat;
    B(:,idx) = A;
end