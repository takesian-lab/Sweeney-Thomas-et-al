function [PeakData, fig1] = simple_check_if_responsive(trials, nBaselineFrames, fs, Z_level, AUC_level, smTime, plot_figure, varargin)
% This will check if a cell is significantly active/responsive in response to sound stimuli
% by detecting significant peaks and troughs in the response relative to Z_level
%
% Argument(s): 
%   trials (T x nFrames) matrix of stim trials
%   nBaselineFrames (int) number of frames to consider as baseline
%   fs (int) framerate
%   Z_level (double) Z-score threshold to count as a significant response
%   AUC_level (double) area under the curve threshold to count as a significant response
%   smTime - moving window (s) for smoothing raw traces. Set to [] for no smoothing
%   plot_figure (0 or 1)
% 
% Returns:
%   Peak_Data (table)
%       - IsResponsive  - 0 or 1 considering both peak and trough
%       - Response_Type - activated, prolonged, suppressed, or none
%
%   The following values exist for both Peak and Trough
%       - HasPeak       - 0 or 1 based on Z_level, AUC_level, and minimum_peak_duration
%       - Peak          - peak or trough value
%       - Peak_3fr      - peak or trough value averaged with the frames around it
%       - Peak_win      - average of defined response_window (same for every trace)
%       - Peak_avg      - average of peak response from onset to offset
%       - Onset         - onset latency of response in s
%       - Latency       - peak latency in s
%       - Width         - width from peak onset to offset
%       - AUC           - area under the curve for peak or area above trough

% Authors: Takesian Lab, 2022, anne_takesian@meei.harvard.edu 

%% Settings

% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% ADDITIONAL PARAMETERS
% Set minimum peak/trough duration (s) to be considered significant
addOptional(p, 'MinimumPeakDuration', 0.2);

% Time (s) after which to consider response (peak) as prolonged
addOptional(p, 'LatePeak', 1); %This is time NOT including baseline

% Time (s) after which to consider response (onset) as prolonged
addOptional(p, 'LateOnset', 1); %Set to equal LatePeak to have no effect

% Set standard window of time(s) to average for Peak_win and Trough_win
addOptional(p, 'ResponseWindow', 1);

% PLOTTING
% Set [m,n] figure dimensions for tiledlayout
addOptional(p, 'FigureDimensions', []); 

% Apply subtitles to figure axes (number of titles = number of plots)
addOptional(p, 'Subtitles', []); 

% Smooth traces for plotting only
addOptional(p, 'SmoothFigure', 0); 

% Handle to figure & axes for plotting graph
addOptional(p, 'FigureHandle', []); 
addOptional(p, 'AxesHandle', []); 

% Plot previously extracted PeakData instead of computing new values
addOptional(p, 'Data', []); 

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

threshold = Z_level; %value to be considered a significant peak
threshold2 = 0; %secondary value to find start and end of peaks

%% Preallocate variables

%Check that trials is in the correct orientation
if size(trials,2) == 1
    trials = trials';
end

if isempty(ops.Data) %Skip all of this if PeakData was previously supplied
    fig1 = [];
    PeakData = table;

    %Preallocate
    PeakData.IsResponsive = false(size(trials,1),1); %initialize as logical using false fxn
    PeakData.ResponseType = strings(size(trials,1),1);
    PeakData.ResponseType(:) = 'undetermined';
    PeakData.HasPeak = false(size(trials,1),1);
    PeakData.Peak = nan(size(trials,1),1);
    PeakData.Peak_3fr = nan(size(trials,1),1);
    PeakData.Peak_win = nan(size(trials,1),1);
    PeakData.Peak_avg = nan(size(trials,1),1);
    PeakData.Peak_Onset = nan(size(trials,1),1);
    PeakData.Peak_Latency = nan(size(trials,1),1);
    PeakData.Peak_Width = nan(size(trials,1),1);
    PeakData.Peak_AUC = nan(size(trials,1),1);
    PeakData.HasTrough = false(size(trials,1),1);
    PeakData.Trough = nan(size(trials,1),1);
    PeakData.Trough_3fr = nan(size(trials,1),1);
    PeakData.Trough_win = nan(size(trials,1),1);
    PeakData.Trough_avg = nan(size(trials,1),1);
    PeakData.Trough_Onset = nan(size(trials,1),1);
    PeakData.Trough_Latency = nan(size(trials,1),1);
    PeakData.Trough_Width = nan(size(trials,1),1);
    PeakData.Trough_AUC = nan(size(trials,1),1);

    %If trials are all NaNs, all zeroes, or empty -- this will skip function
    if ~any(trials)
        return
    end

    %Record trials where all values == NaN. These will be classified as undetermined later
    all_nan_ind = all(isnan(trials),2);
end

%Convert to frames
response_window = ceil(ops.ResponseWindow/(1/fs));
minimum_peak_duration = ceil(ops.MinimumPeakDuration/(1/fs));

%Optional smooth traces with a moving window
if smTime > 0
    smoothval = ceil(smTime/(1/fs)); %window size
    trials = single(smoothdata(trials, 2,'movmean',smoothval)); %TODO: implement zero-phase filter to prevent time shifts (filtfilt)
end

%% PEAK RESPONSES

if isempty(ops.Data)
    
types = {'Peak', 'Trough'};
response = trials(:,nBaselineFrames+1:end);
win = mean(response(:,1:response_window),2,'omitnan'); %Same for Peak and Trough

for h = 1:2 %Peak and Trough

    %Compute peaks or troughs
    [peak, latency, peak_3fr] = getPeakResponse(trials, nBaselineFrames, 'Type', types{h});
    
    if h == 1
        %Make sure peaks are positive
        pos_peak = peak > 0;
        
    elseif h == 2
        %Make sure troughs are negative
        pos_peak = peak < 0;
        
        %Invert sign to measure trough responses
        response = -response;
    end

    %Preallocate and reset values to nan
    [onset, offset, avg, auc, width] = deal(nan(size(peak)));
    traces = cell(size(peak));
    
    for i = 1:size(peak,1)
        if ~pos_peak(i)
            continue
        end

        %Onset is the last point before the peak to cross the threshold
        temp_onset = find(response(i,1:latency(i)) < threshold2, 1, 'last') + 1;

        if isempty(temp_onset)
            onset(i) = 1;
        else
            onset(i) = temp_onset;
        end

        %Offset is the first point after the peak to drop below threshold
        temp_offset = find(response(i, latency(i):end) < threshold2, 1, 'first') - 2;

        if isempty(temp_offset)
            offset(i) = size(response,2);
        else
            offset(i) = temp_offset + latency(i);
        end

        %Width
        width(i) = offset(i) - onset(i);
        if width(i) == 0
            %Offset and onset are the same, set width to 1 frame
            width(i) = 1;
        elseif width(i) < 0
            error('Negative width')
        end
        
        %AUC
        trace = response(i,onset(i):offset(i));
        %If onset and offset are the same, trace will be a single point
        %Measure avg and AUC by bounding point by threshold2
        if length(trace) == 1
            trace = [threshold2 trace threshold];
        end
        avg(i) = mean(trace,'omitnan'); %mean of peak from onset to offset
        trace(isnan(trace)) = threshold2; %trapz function does not work on nans
        auc(i) = trapz(trace - threshold2); %Area under peak above threshold

        %Store peak trace for figure
        if h == 1
            traces{i} = trace;
        elseif h == 2
            traces{i} = -trace;
        end
    end

    %APPLY PEAK CONDITIONS
    
    %Check if peak exceeds Z-score threshold
    thresh_peak = abs(peak) >= threshold;

    %Combine this with pos_peak
    sig_peak = (thresh_peak + pos_peak) == 2;
    
    %Update sig_peak if width < minimum_peak_duration
    sig_peak(width < minimum_peak_duration) = 0;
    
    %Update sig_peak if AUC < AUC_level
    sig_peak(auc < AUC_level) = 0;
   
    %Store peak data
    if h == 1
        PeakData.HasPeak = logical(sig_peak);
        PeakData.Peak = single(peak);
        PeakData.Peak_3fr = single(peak_3fr);
        PeakData.Peak_win = single(win);
        PeakData.Peak_avg = single(avg);
        PeakData.Peak_Onset = single(onset*(1/fs)); %Seconds
        PeakData.Peak_Latency = single(latency*(1/fs)); %Seconds
        PeakData.Peak_Width = single(width*(1/fs)); %Seconds
        PeakData.Peak_AUC = single(auc);
        peak_traces = traces;
    elseif h == 2
        PeakData.HasTrough = logical(sig_peak);
        PeakData.Trough = single(peak);
        PeakData.Trough_3fr = single(peak_3fr);
        PeakData.Trough_win = single(win);
        PeakData.Trough_avg = single(-avg); %Invert trough average
        PeakData.Trough_Onset = single(onset*(1/fs)); %Seconds
        PeakData.Trough_Latency = single(latency*(1/fs)); %Seconds
        PeakData.Trough_Width = single(width*(1/fs)); %Seconds
        PeakData.Trough_AUC = single(auc); %Don't invert AUC
        trough_traces = traces;
    end
end

end

%% Auto-determined response type (activated/prolonged/suppressed)

if isempty(ops.Data)
    
for i = 1:size(trials,1)

    %If no peak or trough detected, activity = none
    if ~PeakData.HasPeak(i) && ~PeakData.HasTrough(i)
        PeakData.ResponseType{i} = 'none';
        continue;
    end
        
    %If there is only a trough and no peak, activity = suppressed
    if ~PeakData.HasPeak(i) && PeakData.HasTrough(i)
        PeakData.ResponseType{i} = 'suppressed';
        continue;
    end
    
    %If there is only a peak, activity = activated or prolonged
    if PeakData.HasPeak(i) && ~PeakData.HasTrough(i)
        PeakData.ResponseType{i} = 'excited';
        continue;
    end
    
    %If there is a peak and a trough, classify based on which has the greater AUC
    if PeakData.HasPeak(i) && PeakData.HasTrough(i)
        
        %If AUC is equal, classify based on whichever peak came first
        if PeakData.Peak_AUC(i) == PeakData.Trough_AUC(i)
            if PeakData.Trough_Latency(i) < PeakData.Peak_Latency(i)
                PeakData.ResponseType{i} = 'suppressed';
            else
                PeakData.ResponseType{i} = 'excited';
            end
        elseif PeakData.Trough_AUC(i) > PeakData.Peak_AUC(i)
            PeakData.ResponseType{i} = 'suppressed';
        else
            PeakData.ResponseType{i} = 'excited';  
        end
    end
end

%Further divide excitatory as activated or prolonged
for i = 1:size(trials,1)
    if strcmp(PeakData.ResponseType{i}, 'excited')
        if PeakData.Peak_Onset(i) < ops.LateOnset
            if PeakData.Peak_Latency(i) > ops.LatePeak
                PeakData.ResponseType{i} = 'prolonged';
            else
                PeakData.ResponseType{i} = 'activated';
            end
        else
            PeakData.ResponseType{i} = 'prolonged';
        end
    end
end

%Catch any cases that don't get classified
if any(strcmp(PeakData.ResponseType{i}, 'undetermined'))
    error('Trial did not get classified -> Figure out why')
end
    
PeakData.IsResponsive(~strcmp(PeakData.ResponseType, 'none')) = 1; %Record responsive traces

%Go back and call trials that were all NaNs undetermined
PeakData.ResponseType(all_nan_ind) = 'undetermined';
PeakData.Peak_Latency(all_nan_ind) = NaN; %fix latency
PeakData.Trough_Latency(all_nan_ind) = NaN; %fix latency

end

%% Plot

if plot_figure
    
    %If data previously supplied, use for plotting
    if ~isempty(ops.Data)
        PeakData = ops.Data;
    end
    
    if isempty(ops.FigureDimensions)
        %Set tiled layout dimensions
        if size(trials,1) > 15
            m = ceil(size(trials,1)/8);
            n = ceil(size(trials,1)/m);
        else
            m = 1;
            n = size(trials,1);
        end
    else
        m = ops.FigureDimensions(1);
        n = ops.FigureDimensions(2);
    end
    
    nPlots = m*n;
    
    %Get XTick labels in seconds
    nFrames = size(trials,2);
    nFramesInS = nFrames*(1/fs);
    XTickLabel = 0:0.5:nFramesInS;
    XTick = round(XTickLabel/(1/fs));
    
    %Find max and min to set ylim
    maxY = max(max(trials));
    minY = min(min(trials));
    
    if maxY < threshold
        maxY = threshold + 0.5;
    end
    
    if minY > -threshold
        minY = -threshold - 0.5;
    end
    
    %Plot figure
    if isempty(ops.FigureHandle)
        %Default plotting
        fig1 = figure('units','normalized','outerposition',[0 0 1 1]); hold on
        tiledlayout(m,n)
    else
        %If user supplied Figure and Axes handles with varargin
        set(ops.FigureHandle, 'CurrentAxes', ops.AxesHandle)
    end
    
    for i = 1:size(trials,1)
        
        y = trials(i,:);
        
        if ops.SmoothFigure
            y = smooth(y,3);
        end
        
        p = nan(1,3);
        t = nan(1,3);
        
        %Adjust for baseline
        p(1) = PeakData.Peak_Onset(i);
        p(2) = PeakData.Peak_Latency(i);
        p(3) = PeakData.Peak_Onset(i) + PeakData.Peak_Width(i);
        t(1) = PeakData.Trough_Onset(i);
        t(2) = PeakData.Trough_Latency(i);
        t(3) = PeakData.Trough_Onset(i) + PeakData.Trough_Width(i);        
    
        %Convert back into frames and add baseline
        p = floor(nBaselineFrames + p*fs);
        t = floor(nBaselineFrames + t*fs);
        
        p(p > length(y)) = length(y);
        t(t > length(y)) = length(y);
        
        if isempty(ops.FigureHandle)
            nexttile; hold on
        end
        plot(y, 'LineWidth', 1.5)
        
        for j = 1:length(p)
            if ~isnan(p(j))
                scatter(p(j), y(p(j)), 'o', 'r', 'LineWidth', 1.5)
            end
        end
        
        for j = 1:length(t)
            if ~isnan(t(j))
                scatter(t(j), y(t(j)), 'o', 'c', 'LineWidth', 1.5)
            end
        end
        
        %Axes
        ylim([minY maxY])
        xlim([1 nFrames])
        set(gca, 'XTick', XTick)
        set(gca, 'XTickLabel', num2str(XTickLabel'))
        
        %Add lines (do this after axes to fill the whole graph)
        hline(mean(y(1:nBaselineFrames),'omitnan'), 'k')
        hline(threshold, 'r')
        hline(-threshold, 'c')
        vline(nBaselineFrames, 'k')
        
        %Bottom row xlabel
        if i > nPlots - n
            xlabel('Seconds');
        end
    
        %First column ylabel
        if ismember(i,1:n:nPlots)
            ylabel('Z-score')
        end
        
        %Subtitles
        if ~isempty(ops.Subtitles)
            subtitle(ops.Subtitles(i), 'FontSize', 8);
        end
        
        %Titles
        if strcmp(PeakData.ResponseType{i}, 'activated') || strcmp(PeakData.ResponseType{i}, 'prolonged')
            title([PeakData.ResponseType{i} ' -  AUC: ' num2str(PeakData.Peak_AUC(i))])
            plot(p(1):p(3), y(p(1):p(3)), 'g', 'LineWidth', 1.5)
        elseif strcmp(PeakData.ResponseType{i}, 'suppressed')
            title([PeakData.ResponseType{i} ' - AUC: '  num2str(PeakData.Trough_AUC(i))])
            plot(t(1):t(3), y(t(1):t(3)), 'g', 'LineWidth', 1.5)
        else
            %Don't add title for None
            %title(PeakData.ResponseType{i})
        end
    end
end
