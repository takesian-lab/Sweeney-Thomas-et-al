function [block] = align_to_silence(block, varargin)
% Find the times associated with any silent period in the block
% Use block.setup.constant to determine trial baseline and duration

% Version history:
% v1 - current version

%% Setup

% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% ADDITIONAL PARAMETERS
% Minimum amount of time to wait after sound presentation
addOptional(p, 'MinSeparationTime', 2.5);

% Plot figures
addOptional(p, 'PlotFigures', 1);

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% Prepare to remove sound_times
%For Markpoints experiment, even spontaneous trials will be associated with PMT shutter clicks, so continue to exclude them

if block.setup.stim_protocol == 12 && ~isfield(block,'MarkpointsExperiment')
    Sound_Time = [];
else
    %Get sound times
    Sound_Time = block.Sound_Time;
    
    %Remove blank trials
    if ~isfield(block,'MarkpointsExperiment')
        [~, ~, ~, ~, blank_ind, ~, ~] = simple_prepare_variables(block);
        Sound_Time(blank_ind == 1) = [];
    end
end
timestamp = block.timestamp;

%Make the beginning of the block a sound (imaging shutter opening)
Sound_Time = [timestamp(1), Sound_Time];

%Make binary timestamp for entire block representing sillence OK and not OK times
isSilent = ones(size(timestamp));
for t = 1:length(Sound_Time)
    [~,a] = min(abs(timestamp - Sound_Time(t))); %sound start
    [~,b] = min(abs(timestamp - (Sound_Time(t)+ops.MinSeparationTime))); %sound end
    isSilent(a:b) = 0;
end

%If PMTs are closed (for markpoints experiments) do not include as silence
isSilent(isnan(timestamp)) = 0;

%Get silent trial duration
baseline = block.setup.constant.baseline_length;
after_stim = block.setup.constant.after_stim;
trial_duration_in_s = baseline + after_stim;
trial_duration_in_fr = round(trial_duration_in_s*round(block.setup.framerate));
MinSeparationTime_in_fr = round(ops.MinSeparationTime*round(block.setup.framerate));

%% Make new silent trial times

%Find each sound stop (silence start) time to maximize silent intervals
isSilent_diff = [0; diff(isSilent)];
SilenceBouts = [find(isSilent_diff == 1); length(timestamp)];

Silence_Time = [];
for s = 1:length(SilenceBouts)-1
    A = SilenceBouts(s);
    B = SilenceBouts(s+1)-2-MinSeparationTime_in_fr;
    
    %Not enough frames for full trial
    if (B-A) < trial_duration_in_fr
        continue;
    end
    
    ticker = A;
    
    %Find as many silent intervals as possible within ISI
    while ticker < B
    
        a = ticker;
        b = ticker + trial_duration_in_fr;

        %if trial exceeds duration of ISI, don't include
        if b < B
            if all(isSilent(a:b))
                Silence_Time = [Silence_Time, timestamp(a)];
            end
        end
    
        ticker = ticker+trial_duration_in_fr+1;
    end
end

Silence_Time = Silence_Time + baseline; %Offset each trial time by baseline length

%% figure

if ops.PlotFigures
    figure; hold on
    area(timestamp, isSilent, 'EdgeColor', 'none', 'Facecolor', [210/255, 248/255, 210/255])
    vline(Silence_Time - baseline,'r')
    legend('Silent periods')
    xlim([timestamp(1) timestamp(end)])
    xlabel('time (s)')    
    sgtitle(block.setup.block_supname)
end

%% Make Silence_Time new SoundTime

if isempty(Silence_Time)
    %If there is not enough silence, this is what will be returned
    block.Sound_Time = [];
    block.parameters.stimLength = [];
    warning('No silent bouts found in block')
else
    block.Sound_Time = Silence_Time;
    block.parameters.stimLength = zeros(size(Silence_Time));
end

%Make this a "Silence block"
block.setup.stim_protocol = 50; %Use same # as for realigned Loco right now
block.AlignedToSilenceOps = ops;
block.parameters.variable1 = ones(size(Silence_Time));
block.parameters.variable2 = ones(size(Silence_Time));
block.parameters.variable3 = [];

%Accommodate markpoints blocks
if isfield(block, 'MarkpointsExperiment')
    block.parameters = rmfield(block.parameters, 'trialsToIgnore');
    block.parameters = rmfield(block.parameters, 'uncagingShutter');
    block = rmfield(block, 'MarkpointsExperiment');
end

%Realign stim
block = align_to_stim(block); %<---- active_trials variable gets adjusted here. If there is no silence_time, this will replace aligned_stim with an empty struct

%% Neural responses (check that align_to_stim worked well)
%After align_to_stim, loco has been put on same timescale as 2P data

nBaselineFrames = round(round(block.setup.framerate)*block.setup.constant.baseline_length);

if ~isempty(Silence_Time) && ops.PlotFigures
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