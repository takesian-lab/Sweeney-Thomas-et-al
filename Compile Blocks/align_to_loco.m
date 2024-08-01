function [block] = align_to_loco(block, varargin)
% Find the times associated with the start of loco bouts

% Version history:
% v1 - current version

%% Setup

% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% ADDITIONAL PARAMETERS
% Discard trials if a sound was played anytime during loco trial
addOptional(p, 'RemoveSounds', 0);

% Sounds during baseline period up to RemoveSoundWindow will be discarded
addOptional(p, 'RemoveSoundWindow', block.setup.constant.after_stim);

% Add a positive or negative time delay for each trial start with respect to the loco bout
addOptional(p, 'OnsetShift', 0);

% Maximum running speed allowed to be included as a trial
addOptional(p, 'MaxLoco', inf);

% Minimum speed (cm/s) threshold to be considered running
addOptional(p, 'BoutThreshold', 2);

% Minimum duration (s) of a loco bout to be included as a trial
addOptional(p, 'MinBoutDuration', 1);

% Max duration (s) of a loco bout to be included as a trial
addOptional(p, 'MaxBoutDuration', inf);

% Plot figures
addOptional(p, 'PlotFigures', 1);

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% Prepare to remove sound_times

if block.setup.stim_protocol == 12
    ops.RemoveSounds = 0; %Dont remove sounds for spontaneous data
else
    %Get sound times
    Sound_Time = block.Sound_Time;
    
    %Remove blank trials
    [~, ~, ~, ~, blank_ind, ~, ~] = simple_prepare_variables(block);
    Sound_Time(blank_ind == 1) = [];
end

%% find the loco bouts [METHOD #1]

% TEST: Use raw loco trace (10fps) to detect bouts
% loco = block.loco_activity;
% timetrace = block.locomotion_trace;
% fs = round(1/((timetrace(end)-timetrace(1))/length(timetrace)));
% nBaselineFramesLoco = round(block.setup.constant.baseline_lengt*fs); %CAREFUL, LOCO FS IS DIFFERENT THAN 2P FS

% INSTEAD: Upsample raw loco trace to match Bruker timestamp. This will make align_to_stim more accurate
loco_activity = block.loco_activity;
loco_activity(loco_activity < block.setup.constant.locoThresh) = 0; %Remove the noise floor BEFORE smoothing, otherwise traces will go back to 0 if we correct after
[timetrace, ~, loco] = match_fluor_loco_timestamps(block.F7, loco_activity, block.timestamp, block.locomotion_trace);
timetrace = timetrace';
nBaselineFramesLoco = round(block.setup.framerate*block.setup.constant.baseline_length);
nAfterStimFramesLoco = round(block.setup.framerate*block.setup.constant.after_stim);
nOnsetShiftFramesLoco = round(block.setup.framerate*ops.OnsetShift);
nRemoveSoundWindowFramesLoco = round(ops.RemoveSoundWindow*block.setup.constant.after_stim);

%If user specified an ops.OnsetShift less than 0, extend the baseline period by this many seconds
if ops.OnsetShift < 0
    nBaselineFramesLoco = round(block.setup.framerate*(block.setup.constant.baseline_length + abs(ops.OnsetShift)));
end

raw_loco = loco; %Store for plotting
loco = smooth(loco); %Smooth a bit to make loco bouts "stick together" more *WARNING* more smoothing makes align_to_stim less accurate (loco onset time gets earlier)
isRunning = loco > ops.BoutThreshold;
bout_timestamp = timetrace(isRunning);

% loop through possible bouts and make sure there is enough baseline and after_stim
BoutTime = [];
BoutDuration = [];
IsBout = zeros(size(loco)); %Use to shade bout area in plot later

if any(bout_timestamp)
    bout_diff = [bout_timestamp(1); diff(bout_timestamp)]; %seconds between each frame where mouse is running > ops.BoutThreshold
    possible_bouts = bout_timestamp(bout_diff > ops.MinBoutDuration);
    
    for b = 1:length(possible_bouts)

        % time where bout reached minimum speed threshold
        bout_time = possible_bouts(b);
        [~, bout_ind] = min(abs(timetrace - bout_time));

        %find time when mouse started running to use as acutal bout start
        if bout_ind > 1 %If bout was not already ongoing at frame 1..
            bout_ind = find(loco(1:bout_ind) == 0, 1, 'last') + 1; %Find the last index where mouse was not running and set the next frame as the start of the bout
        end

        %find time when mouse stopped running to use as bout duration
        bout_end_ind = find(loco(bout_ind:end) == 0, 1, 'first') - 1 + bout_ind; %Find the first index after the mouse starts running where they stop again, and set the previous frame as the end of the bout
        if isempty(bout_end_ind)
            bout_end_ind = length(loco); %If no 0 found, it means mouse was running until loco recording ended
        end
        bout_duration = timetrace(bout_end_ind) - timetrace(bout_ind);

        %If the bout is not long enough, continue to next possible bout
        if bout_duration < ops.MinBoutDuration
            continue;
        end

        %If the bout is too long, continue to next possible bout
        if bout_duration > ops.MaxBoutDuration
            continue;
        end
        
        %If it is just right, determine loco-aligned trial start and end times
        baseline_ind = bout_ind - nBaselineFramesLoco + 1 + nOnsetShiftFramesLoco;
        after_stim_ind = bout_ind + nAfterStimFramesLoco + nOnsetShiftFramesLoco;
        if after_stim_ind > length(timetrace)
            after_stim_ind = length(timetrace);
        end
        sound_window_ind = bout_ind + nRemoveSoundWindowFramesLoco + nOnsetShiftFramesLoco;
        if sound_window_ind > length(timetrace)
            sound_window_ind = length(timetrace);
        end
        
        %Unfortunately, if we don't have enough baseline, continue to next possible bout
        if baseline_ind <= 0
            continue;
        end

        %If ops.RemoveSounds, check that there is no sound in the running trial
        if ops.RemoveSounds
            
            T1 = timetrace(baseline_ind); %baseline time
            T2 = timetrace(sound_window_ind); %loco time + after_stim

            %Don't include as a bout if a sound was found
            if any(((Sound_Time >= T1) + (Sound_Time <= T2)) == 2)
                continue;
            end
        end
        
        %Check that MaxLoco is not exceeded
        if any(loco(baseline_ind:after_stim_ind) > ops.MaxLoco)
            continue;
        end
        
        %Finally, check that there is no running in the baseline period.
        if any(loco(baseline_ind:bout_ind) > ops.BoutThreshold)     
            continue;
        end
        
        %If you've made it this far, you get to count as a loco bout!
        BoutTime = [BoutTime, timetrace(bout_ind)];
        BoutDuration = [BoutDuration, bout_duration];
        IsBout(bout_ind:bout_end_ind) = 1;
        
    end
end

%% figure

if ops.PlotFigures
    figure;
    
    maxy = max(loco)+ 1;
    
    subplot(2,1,1); hold on
    plot(block.locomotion_trace, block.loco_activity(1:length(block.locomotion_trace)), 'k', 'Linewidth', 2)
    %plot(timetrace, raw_loco, 'k', 'Linewidth', 1)
    ylim([0 maxy]) %set before vline
    if ~isempty(BoutTime)
        vline(possible_bouts, 'r')
    end
    %plot(timetrace, loco, 'k', 'Linewidth', 2) %plot again 
    %legend({'Raw loco trace', 'Smoothed loco trace'})
    legend('Raw loco trace')
    xlim([timetrace(1) timetrace(end)])
    hline(ops.BoutThreshold)
    ylabel('speed (cm/s)')
    xlabel('time (s)')
    title('Method 1: Possible bouts')
    
    subplot(2,1,2); hold on
    ylim([0 maxy]) %set before vline
    if ~isempty(BoutTime)
        area(timetrace, IsBout*maxy, 'EdgeColor', 'none', 'Facecolor', [210/255, 248/255, 210/255])
        vline(BoutTime,'r')
    end
    plot(timetrace, loco, 'k', 'Linewidth', 2)
    legend({'Detected bout', 'Upsampled loco trace'})
    xlim([timetrace(1) timetrace(end)])
    hline(ops.BoutThreshold)
    ylabel('speed (cm/s)')
    xlabel('time (s)')
    ylim([0 maxy])
    title('Method 1: Final bouts')
    
    sgtitle(block.setup.block_supname)
end

%% Make BoutTime new SoundTime

if isempty(BoutTime)
    %If mouse is not running during block, this is what will be returned
    block.Sound_Time = [];
    block.parameters.stimLength = [];
    warning('No running bouts found in block')
else
    block.Sound_Time = BoutTime;
    block.parameters.stimLength = BoutDuration;
end

%Make this a "Loco block"
%block.locomotion_trace = loco; %do we maybe want to replace loco trace with smoothed one???
block.setup.stim_protocol = 50;
block.parameters.variable1 = ones(size(BoutTime));
block.parameters.variable2 = ones(size(BoutTime));
block.parameters.variable3 = [];
block.AlignedToLocoOps = ops;
sound_time_delay = zeros(size(block.Sound_Time)) + ops.OnsetShift;
block = align_to_stim(block, sound_time_delay); %<---- active_trials variable gets adjusted here. If there is no running, this will replace aligned_stim with an empty struct

%% Neural responses (check that align_to_stim worked well)
%After align_to_stim, loco has been put on same timescale as 2P data

nBaselineFrames = round(block.setup.framerate*block.setup.constant.baseline_length);

if ~isempty(BoutTime) && ops.PlotFigures
    figure;
    
    subplot(4,1,1); hold on
    for v = 1:size(block.aligned_stim.velocity,1)
        plot(block.aligned_stim.velocity(v,:), 'k');
    end
    plot(mean(block.aligned_stim.velocity,1),'k','Linewidth',2)
    vline(nBaselineFrames, 'r')
    xlim([1 size(block.aligned_stim.velocity,2)])
    title('Loco traces')
    
    subplot(4,1,2)
    imagesc(block.aligned_stim.velocity)
    set(gca, 'YTick', 1:size(block.aligned_stim.velocity,2))
    xlim([1 size(block.aligned_stim.velocity,2)])
%     try
%         ylim([0.5 size(block.aligned_stim.velocity,1)-0.5])
%     catch
%         %do nothing
%     end
    vline(nBaselineFrames,'r')
    %colorbar
    ylabel('Loco trials')
    
    subplot(4,1,3); hold on
    mean_zscore = squeeze(mean(block.aligned_stim.zscore,2));
    for z = 1:size(mean_zscore,1)
        plot(mean_zscore(z,:), 'g', 'Linewidth', 0.25);
    end
    plot(mean(mean_zscore,1),'k','Linewidth',2)
    xlim([1 size(mean_zscore,2)])
    vline(nBaselineFrames, 'r')

    title('Averaged cell traces')
    
    subplot(4,1,4); hold on
    imagesc(mean_zscore)
    xlim([1 size(mean_zscore,2)])
    ylim([0.5 size(mean_zscore,1)-0.5])
    vline(nBaselineFrames, 'r')
    %colorbar
    caxis([0 2])
    ylabel('Cells')
    xlabel('Frames')
    
    sgtitle(block.setup.block_supname)
end