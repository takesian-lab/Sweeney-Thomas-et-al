function [CalciumData, summary_data, fig1] = simple_Ca_detection_spontaneous(timestamp, F7, loco, fs, cellIndex, Loco, plot_graphs)
% Detect and measure calcium transients in spontaneous data

% Version log:
%   - V1 = current version based off of Ca_detection_spontaneous

% Argument(s): 
%   timestamp, F7, loco - already trimmed and downsampled using match_fluor_loco_timestamps
%   fs = framerate
%   cellIndex (int) = cell order in F7 matrix
%   Loco = 'All', 'Running', 'NotRunning'
%   plot_graphs = 0 or 1 to plot figure
%
% Returns:
%   CalciumData (struct)
%   summary_data (table)
%     block_duration - duration (s) of block that was included
%     loco_percent - percentage of block when mouse is running (in case this affects FR)
%     N - number of transients
%     tr_rate - Rate of transients (Hz)
%     tr_med - Median transient amplitude
%     tr_mean - Mean transient amplitude
%     tr_STD - STD of transient amplitudes
%     tr_duration_med - Median transient duration (s)
%     tr_duration_mean - Mean transient duration (s)
%     tr_duration_STD - STD of transient duration (s)
%   fig1 - figure handle for trace + fits
%   fig2 - figure handle for summary statistics
% 
% Notes:
%
% TODO: 
% Search 'TODO'

%% Ca Detect Parameters (change to adjust Ca detection probability)
% Maryse is using mwin 5, bottom 30, thres 60, SD peak 3, SD base 1, duration 1

if strcmp(Loco, 'All')
    loco_filter = 0;
elseif strcmp(Loco, 'Running')
    loco_filter = 1;
elseif strcmp(Loco, 'NotRunning')
    loco_filter = 2;
end
                
mwin = 5; % length of moving window in seconds
bottom_percentile = 30; %bottom percentile of fluoresence used for baseline
thresh_percentile = 60;
SD_peak = 3; %threshold for peaks
SD_base = 1; % threshold to return to baseline
maxRiseDuration = 1; %seconds

thA = 1; %1.5; %zero crossings for threshold SD, i.e. when does the trace cross the threshold SD
thB = .5; %1; %zero crossings for baseline SD, i.e. when does the trace cross the threshold SD

params = []; 
params.threshA = thA; 
params.threshB = thB; 

plot_test_graphs = 0; %For troubleshooting

%% Setup

%Figures will return as [] if plot_graphs is 0
[fig1] = deal([]);

%Variables will return as NaN if no transients can be detected
[tr_med, tr_mean, tr_STD, tr_duration_med, tr_duration_mean, tr_duration_STD] = deal(single(nan));

%Number of transients and transient frequency will return as 0 if no transients are detected
[N, tr_rate] = deal(single(0));

CalciumData = struct;
CalciumData.cellindex = cellIndex;
CalciumData.Trace = [];
CalciumData.EventTrace = [];
CalciumData.RiseTrace = [];
CalciumData.FitTrace = [];
    
summary_data = table; 
summary_data.block_duration = single(nan); 
summary_data.loco_percent = single(nan); 
summary_data.N = N;
summary_data.tr_rate = tr_rate; 
summary_data.tr_med = tr_med; 
summary_data.tr_mean = tr_mean; 
summary_data.tr_STD = tr_STD; 
summary_data.tr_duration_med = tr_duration_med; 
summary_data.tr_duration_mean = tr_duration_mean; 
summary_data.tr_duration_STD = tr_duration_STD; 

fr = round(fs,0); %framerate
maxRiseDuration = maxRiseDuration*fr; %Convert max rise duration to frames

%% Focus on just this cell

F7 = F7(cellIndex,:);

%% Adjust F7 so that there are no negative values

if any(F7 < 0)
    min_F7 = min(F7); %Find minimum value
    F7 = F7 + abs(min_F7); %Add abs(min value) to full trace such that min is now 0
end

%% OPTIONAL: Remove locomotor activity from traces (or look only at loco activity)

locoThresh = 0.7; %MAGIC NUMBER

%Determine when mouse is running
isLoco = loco > locoThresh;

if plot_test_graphs
    figure; hold on

    subplot(2,1,1); hold on
    title('Locomotor activity')
    ylabel('Activity (cm/s)') 
    plotmax = max(loco)+1; %Used in area
    area(timestamp,(isLoco > 0)*plotmax, 'EdgeColor', 'none', 'Facecolor', [210/255, 248/255, 210/255])
    plot(timestamp, loco, 'LineWidth', 0.25)
    xlim([timestamp(1) timestamp(end)])
    try ylim([0 plotmax]); catch;  end
    hline(locoThresh)

    subplot(2,1,2); hold on
    title('Normalized cell activity')
    ylabel('F7')
    xlabel('Seconds')
    plotmax = max(F7./max(F7))+.25; %Used in area
    area(timestamp,(isLoco > 0)*plotmax, 'EdgeColor', 'none', 'Facecolor', [210/255, 248/255, 210/255])
    plot(timestamp, F7./max(F7), 'LineWidth', 0.25)
    ylim([0 plotmax])
    xlim([timestamp(1) timestamp(end)])
end

%Calculate percentage of block spent running
loco_percent = (sum(isLoco)/length(isLoco))*100;

%Change F7 and timestamp only if loco_filter > 0
if loco_filter == 1 %Remove loco bouts
    F7 = F7(~isLoco);
    %timestamp = timestamp(~isLoco); (has gaps)
    timestamp = (1:length(F7))*(1/fr); %Create new artificial timestamp without gaps
elseif loco_filter == 2 %Keep only loco bouts
    F7 = F7(isLoco);
    %timestamp = timestamp(isLoco); (has gaps)
    timestamp = (1:length(F7))*(1/fr); %Create new artificial timestamp without gaps
end

%Return if there is less than 1 second of data to analyze
if isempty(timestamp) || (length(timestamp) <= fs )
    warning('not enough data to analyze for this loco type')
    return
end

%Compute block duration
block_duration = timestamp(end) - timestamp(1);

%% Start from the beginning of the Harnett lab code/Wisam's scripts...
%To compute DF/F, the baseline F was estimated as the 10th percentile of the fluorescence using a 160 s rolling-window.

% determine size of moving window
win = fr*mwin;
% windows work better on odd sized - windows, correct for that here
s = (-1)^win;
if s == 1
    win = win+1;
end

%Window size must be <= size of X and >=1
%Data may be shorter than we expect due to removing loco segments
%Allow to readjust the window size here (to represent full trace)
if win >= length(F7)
    win = length(F7) - 1;
end

%To compute DF/F, the baseline F was estimated as the 10th percentile of the fluorescence using a predetermined rolling-window.
[Fo] = running_percentile(F7, win, bottom_percentile); %(data, window size in frames, percentile)
dffo = double((F7-Fo')./Fo');

%% Remove transients from the trace by thresholding all values above certain percent dffo
d60 = running_percentile(dffo, win, thresh_percentile); %%(data, window size in frames, percentile)

th = nan(size(dffo));
for i = 1:length(dffo)
    if dffo(i) < d60(i)
        th(i) = dffo(i);
    end
end

if plot_test_graphs
    figure
    subplot(2,1,1); hold on
    plot(F7(1,:))
    plot(Fo')
    legend({'F7', 'Fo'})
    
    subplot(2,1,2); hold on
    plot(dffo)
    plot(d60)
    plot(th)
    legend({'dffo', 'd60', 'th'})
end

%% calculate the SD for thresholding

M = movstd(th,win,'omitnan');
range = M*SD_peak; 
ret = M*SD_base;
limPos = M+range;
% limNeg = M-range;
retBase = M+ret;

if plot_test_graphs
    figure; hold on
    plot(M)
    plot(retBase)
    plot(limPos)
    legend({'STD of threshold percentile', 'Threshold to return to baseline', 'Threshold for peaks'})
end

%% low pass filter the raw signal

% Transient filter cutoff frequency
CutoffFrequency = 2; % Hz

% Lowpass filter order (chosen arbitrarily)
LPFilterOrder = 2; %12

[b, a] = butter(LPFilterOrder,CutoffFrequency/(fr/2));

%The length of the input X must be more than three times the filter order, defined as max(length(B)-1,length(A)-1).
if length(dffo) <= 3*max(length(b)-1,length(a)-1)
    warning('trace not long enough to compute transients')
    return    
end
    
transientFiltered = filtfilt(b,a,dffo);
% differentiate
dDif = diff(transientFiltered);
transDiffFilt = filtfilt(b,a,dDif); %filtered derivative

%% segment the data using Wisam's code, updated by CGS March 2021
% run this for both threshold A (3sd) and threshold b (1 sd) - SD values are set at the top of this code
% detectThreshold_CGS used to be in a separate script. It is now embeddded in this funcion

transDiffFilt = diff(transientFiltered); 

%zero crossings for baseline SD, i.e. when does the trace cross the baseline SD
[thB_idx] = detectThreshold_spontaneous(transientFiltered, 'fixed', thB);

%zero crossings for threshold SD, i.e. when does the trace cross the threshold SD
[thA_idx] = detectThreshold_spontaneous(transientFiltered, 'fixed', thA);

% Calculate the zero-crossings on the differentiated transient signal

% These are zero crossing indices
[ZeroCrossingIndices] = detectThreshold_spontaneous(transDiffFilt, 'fixed', 0); % signal, threshold type, threshold (zero in this case)

%what does this look like?
if plot_test_graphs
    figure; hold on
    plot(transientFiltered)
    plot(M)
    plot(limPos,'-k')
    plot(retBase, '-g')
    scatter(thA_idx, transientFiltered(thA_idx), 'r', 'filled');
    scatter(thB_idx, transientFiltered(thB_idx), 'r');
    scatter(ZeroCrossingIndices, transientFiltered(ZeroCrossingIndices), 'MarkerEdgeColor', 'k')
end

%% Temporal Segmentation of the transients
% Here is where we define our "Ca events" and the "Ca onset transients" 
% We will assign the time stamps associated with each, and segment the
% original trace using them.

% This is the number of events that might qualify as Ca-event
numberOfCandidateBouts = length(thB_idx)-1;

%If no candidate events, end code here
if numberOfCandidateBouts <= 0
    CalciumData.Trace = transientFiltered;
    warning('no calcium events found in cell')
    return
else
    % SegmentTransients used to be in a separate script. It is now embeddded in this funcion
    CalciumData = SegmentTransients_spontaneous(thB_idx, thA_idx, ZeroCrossingIndices, transientFiltered, fr, maxRiseDuration);
    
     if plot_graphs
        fig1 = figure; hold on
        plot(timestamp,CalciumData.Trace, '-k');
        plot(timestamp,CalciumData.EventTrace,'-b','LineWidth',2);
        plot(timestamp,CalciumData.RiseTrace,'-r','LineWidth',2);
        plot(timestamp,CalciumData.FitTrace,'-g','LineWidth',2);
        xlabel('time (s)')
        ylabel('df/f')
     end
end       

%% Calculate summary statistics

%Amplitude and frequency of transients
amplitudes = [CalciumData.Amplitude{:}];

if isempty(amplitudes)
    warning('no calcium events found in cell')
    return
end

N = length(amplitudes); %Number of transients
tr_med = median(amplitudes); %Median amplitude
tr_mean = mean(amplitudes); %Mean amplitude
tr_STD = std(amplitudes); %STD of amplitudes
tr_rate = N/timestamp(end); %Frequency (Hz)

%Duration of transients in seconds
durations = nan(size(amplitudes));
for i = 1:N
    bout = CalciumData.motorBoutTimeStamps{i};
    bout_in_s = timestamp(bout);
    durations(i) = bout_in_s(end) - bout_in_s(1);
end
tr_duration_med = median(durations);
tr_duration_mean = mean(durations);
tr_duration_STD = std(durations);

summary_data.block_duration =  single(block_duration); 
summary_data.loco_percent = single(loco_percent); 
summary_data.N = single(N);
summary_data.tr_rate = single(tr_rate); 
summary_data.tr_med = single(tr_med); 
summary_data.tr_mean = single(tr_mean); 
summary_data.tr_STD = single(tr_STD); 
summary_data.tr_duration_med = single(tr_duration_med); 
summary_data.tr_duration_mean = single(tr_duration_mean); 
summary_data.tr_duration_STD = single(tr_duration_STD); 


if plot_test_graphs
    figure;
    subplot(1,2,1); hold on
    bar(tr_mean)
    errorbar(tr_mean, tr_STD)
    scatter(ones(size(amplitudes)), amplitudes)
    ylabel('Amplitude (df/f)')
    title(['N = ' num2str(N) ' transients'])
    
    subplot(1,2,2); hold on
    bar(tr_duration_mean)
    errorbar(tr_duration_mean, tr_duration_STD)
    scatter(ones(size(durations)), durations)
    ylabel('Duration (s)')
    title(['Frequency = ' num2str(tr_rate) ' Hz'])
end
     
end

function crossings = detectThreshold_spontaneous(signal, thresholdType, threshold)
% detectThreshold finds the domain index where the input signal crosses a user defined threshold.
% 
% Written by Wisam Reid - June 2020 - wisam@g.harvard.edu
%  
% This function locates both upward and downward crossings by looking for a 
% sign change in the threshold shifted signal (thresholdShiftedSignal). 
% Where, thresholdShiftedSignal = signal - threshold 
% 
% This function is vectorized to optimize speed.
%
% Arguments:       
%                  signal: (double) vector
%           thresholdType: (string) options: {'STD','Fixed'} 
%               threshold: (double) multiplier of the std of signal or a fixed threshold
%
% Returns:      
%               crossings: (double) 2 column array where the first column contains 
%                          upward crossings and the second contains downward crossings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Example Use:
%              
% dt = 0.001; % (1 ms in seconds)
% time = 0:dt:20; % (seconds)
% triangle = [0:dt:10 flip(0:dt:10-dt)];
% threshold = 3.01;
% crossings = detectThreshold(triangle, 'Fixed', threshold);
% figure; % plot demo
% plot(time,triangle,'-b','DisplayName','Signal','LineWidth',2)
% hold on
% yline(threshold,'-g','DisplayName','Threshold','LineWidth',2);
% xline(time(crossings(1)),'-r','DisplayName','Positive Crossing','LineWidth',2);
% xline(time(crossings(2)),'-k','DisplayName','Negative Crossing','LineWidth',2);
% p = plot(time(crossings(1)),threshold,'ro','MarkerFaceColor','Red');
% p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p = plot(time(crossings(2)),threshold,'ko','MarkerFaceColor','Black');
% p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% xlim([time(1) time(end)]);
% title('detectThreshold Demo')
% xlabel('Time (seconds)') 
% ylabel('Amplitude (generic units)')
% legend('Location','south')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Are we using a standard deviation based or fixed threshold?
switch thresholdType
    case 'STD'
        % The threshold as a factor of the std of the input signal
        threshold = threshold*std(signal);
    case 'Fixed'
        threshold = threshold;
end

% Compute the thresholded shifted signal. 
% We will use this to look for sign changes.
thresholdShiftedSignal = signal - threshold;

% Now we will locate the crossings by looking for where consecutive elements of
% signal change sign 
% i.e. their product is negative
crossings = find(thresholdShiftedSignal.*circshift(thresholdShiftedSignal,1) < 0);

% Now we will pair the upward and downward crossing into a cell array:

% 1) Reshape the array into a 2 column array where the first column
% contains upward crossings and the second contains downward crossings
% crossings = reshape(crossings,2,[])';

% 2) Combine the upward and downward crossings into each element of the cell array 
% crossings = num2cell(crossings,2);

end

function TempCalciumData = SegmentTransients_spontaneous(thB_idx, thA_idx, ZeroCrossingIndices, transientFiltered, fr, maxRiseDuration)
%Based off SegmentTransients by CGS

numberOfCandidateBouts = length(thB_idx)-1;

frIdx = 1:length(transientFiltered);

%preallocate
TempCalciumData = struct;
TempCalciumData.motorBoutTimeStamps = {};
TempCalciumData.motorBouts = {};
TempCalciumData.riseEventTimeStamps = {};
TempCalciumData.riseEvents = {};
TempCalciumData.Amplitude = {};
            
% Loop over candidate motor bouts
count = 1; %count will be for the successful bouts
for i = 1:numberOfCandidateBouts

    % Grab the frames for all threshold A crossings within the current
    % candidate motor bout time interval. Find the frames that cross
    % threshold A and are between thresholdB and thresholdB+1
    thA_fr = thA_idx(thA_idx > thB_idx(i) & thA_idx < thB_idx(i+1));
    
    % Is this a Calcium event?
    if ~isempty(thA_fr)
        % This is now considered to be a Ca event
        
        % This is the starting frame for the current event
        Ca_start = thB_idx(i);
        
        % This is the ending frame for the Ca event
        Ca_end = thB_idx(i+1);
        
        % Grabbing all frames where the increase in Ca activity is zero and
        % the amplitude is above threshold A (or threshold B?)
        riseZero_fr = sort(ZeroCrossingIndices(ZeroCrossingIndices >= Ca_start(1) & ZeroCrossingIndices <= thA_fr(end)),'ascend');
       
        % set rise event end to the max value in bout:
        [maxPoint,maxIDX] = max(transientFiltered(Ca_start:Ca_end));
        bout = Ca_start:Ca_end;
        riseEventEndTS = bout(maxIDX);
        
        %Check rise and decay durations to eliminate bouts that are too broad
        riseDuration = length(Ca_start:riseEventEndTS);
        fallDuration = length(riseEventEndTS:Ca_end);
        
        minEventDuration = 2;  %This is a very lenient condition, only set to eliminate bouts that might be an artifact
        
        if riseDuration < minEventDuration || riseDuration > maxRiseDuration
            continue
        end
        
        if fallDuration < minEventDuration
            continue
        end
        
        %If we made it to this point, this will be considered a successful bout -> Record data
        TempCalciumData.motorBoutTimeStamps{count} = (Ca_start:Ca_end);
        TempCalciumData.motorBouts{count} = transientFiltered(1,Ca_start:Ca_end);
        TempCalciumData.riseEventTimeStamps{count} = frIdx(frIdx >= Ca_start & frIdx <= riseEventEndTS);
        TempCalciumData.riseEvents{count} = transientFiltered(riseZero_fr:riseEventEndTS);
        TempCalciumData.Amplitude{count} = maxPoint;
        
        count = count + 1;
    end
end
%% does a segmented curve fit the kinetics of gcamp?

% can we calculate this parameter solely on the decay of the candidate calcium event?
% first, let's pull out the decay of each Ca event:

CaDecay = cell(size(TempCalciumData.motorBoutTimeStamps));

for i = 1:length(TempCalciumData.motorBoutTimeStamps)
    CaDecay{1,i} = [];
    decayFit.fitTrace{1,i} = [];
    decayFit.gof{1,i} = [];
    decayFit.tau{1,i} = [];
    decayFit.t_half{1,i} = [];
    if ~isempty(TempCalciumData.motorBoutTimeStamps{1,i})
        RiseEnd = length(TempCalciumData.riseEventTimeStamps{1,i});
        DecStart= RiseEnd +1;
        Amplitude = TempCalciumData.motorBouts{1,i}(RiseEnd);
        CaDecay{1,i} = TempCalciumData.motorBoutTimeStamps{1,i}(DecStart:end); % store falling phase
        x = 1:length(TempCalciumData.motorBoutTimeStamps{1,i}(DecStart:end));
        fit_vector = TempCalciumData.motorBouts{1,i}(DecStart:end); %the falling phase
        
        %In case the data can't be fit for some reason, don't worry for now
        try
            [fitfallingPhase,gof] = fit(x',fit_vector','exp1'); %fit data
            realfitfallingPhase = squeeze(fitfallingPhase(x));
            datatauFall1 = (1/(fitfallingPhase.b*-1))./fr;
            tHalf = datatauFall1*log(2);
            TempCalciumData.decayFit.fitTrace{i} = realfitfallingPhase;
            TempCalciumData.decayFit.gof{i} = gof;
            TempCalciumData.decayFit.tau{i} = datatauFall1;
            TempCalciumData.decayFit.t_half{i} = tHalf;
            TempCalciumData.EventAmplitude{i} =  Amplitude;
        catch
            TempCalciumData.decayFit.fitTrace{i} = nan;
            TempCalciumData.decayFit.gof{i} = nan;
            TempCalciumData.decayFit.tau{i} = nan;
            TempCalciumData.decayFit.t_half{i} = nan;
            TempCalciumData.EventAmplitude{i} =  Amplitude;
        end
    end
end

TempCalciumData.CaDecayTimeStamps = CaDecay;

%% restructure for single trace
% 'TempCalciumData' is set up for a single event/bout, but I have multiple events
% per trace. Here, I restructure the outputs and concatenate the segmented traces

TempCalciumData.Trace = transientFiltered;
TempCalciumData.EventTrace = nan(1,length(transientFiltered));
TempCalciumData.RiseTrace = nan(1,length(transientFiltered));
TempCalciumData.FitTrace = nan(1,length(transientFiltered));

for i = 1:length(TempCalciumData.motorBoutTimeStamps)
    Aa = TempCalciumData.motorBoutTimeStamps{:,i}(:);
    for b = 1:length(Aa)
        TempCalciumData.EventTrace(Aa(b)) = TempCalciumData.motorBouts{:,i}(b);
    end
    
    Cc = TempCalciumData.riseEventTimeStamps{:,i}(:);
    for d2 = 1:length(Cc)
        TempCalciumData.RiseTrace(Cc(d2)) = TempCalciumData.motorBouts{:,i}(d2);
    end
    
    Ee = TempCalciumData.CaDecayTimeStamps{1,i}(:);
    for f = 1: length(Ee)
        TempCalciumData.FitTrace(Ee(f)) = TempCalciumData.decayFit.fitTrace{:,i}(f);
    end
end

end