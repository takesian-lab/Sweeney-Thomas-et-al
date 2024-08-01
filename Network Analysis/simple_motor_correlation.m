function    [fig1, MotorCorrData, cc, lag_in_s] = ...
    simple_motor_correlation(traces, ExtractedData, MotorCorrOps, tableRow, varargin)
% Correlation with running activity [AND OTHER BEHAVIOR TOO!]

% This function outputs cross correlations of full traces of neurons in
% a network with locomotor activity. Written by A. Takesian based on
% trace correlation and M. Thomas x_corr code.

% Argument(s): 
% traces - full raw traces from cell(s) in a block (n cells x frames)
% ExtractedData - output from simple_extract_data
%  Traces are already resampled and cropped to be the same # of frames in ExtractedData
% MotorCorrOps - parameters for computing correlation:
%  .smTime
%  .p_value
%  .maxlag_in_s
%  .shuffles
%  .zscore
%  .fs
% tableRow - index of TrialData row to analyze
% varargin - TraceToAnalyze: 'Loco' (default), 'XY', 'Residuals', 'zcorr', 'Whisker',  or 'Pupil'
% varargin - PlotFigure: 0 (default) or 1 to plot figure of motor traces

% Returns MotorCorrTraces with:
%   cc_max = max of the xcorr functions
%   cc_min = min of the xcorr functions
%   cc_sign = whether max/min is greater than 0 (largest abs. value is used to decide)
%   lag_max = lag time of max (s)
%   lag_min = lag time of min (s)
%   cc_zero = value at zero time lag of the xcorr funtions
%   cc_zero_sign = whether cc_zero is greater than 0
%   " " _z = Z-test result (0 or 1) for comparing each of the above values against control distribution
%   " " _p = p-value of the test
%   " " _r = z-score from the test

% VERSION HISTORY
% v1 archived 6/26/23 by MET
% v2 = current version using tak_compute_xcorr, improves how max and min significance is determined
% maryse added whisker and pupil 8/25/23

%% ------- parse varargin

p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% Choose which trace to correlate with neural activity: 'Loco', 'XY', 'Residuals', 'zcorr,' 'Whisker', or 'Pupil
addOptional(p, 'TraceToAnalyze', 'Loco');

% Plot figure of all movement traces?
addOptional(p, 'PlotFigure', 0);

% Plot figure of all xcorr traces?
addOptional(p, 'PlotXCorrFigures', 0);

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% If trials are all NaNs, all zeroes, or empty -- this will skip function

if ~any(traces)
    fig1 = []; MotorCorrData = []; cc = []; lag_in_s = [];
    return
end

%% Apply options from MotorCorrOps

% Optional smooth traces with a moving window
if MotorCorrOps.smTime > 0
    smoothval = ceil(MotorCorrOps.smTime/(1/MotorCorrOps.fs)); %window size
    traces = single(smoothdata(traces, 2,'movmean',smoothval));
end

% Convert raw traces to zscore or df/f
if MotorCorrOps.zscore
    if size(traces,1) > 1
        traces = zscore(traces,'',2);
        if any(any(isnan(traces))); error('NaNs found -- fix so that zscoring will work'); end
    else
        %For 1xnFrames traces, remove nans (e.g. whisker/pupil traces)
        traces_withoutnans = traces(~isnan(traces));
        traces_withoutnans_zscored = zscore(traces_withoutnans,'',2);
        traces_zscored = traces;
        traces_zscored(~isnan(traces)) = traces_withoutnans_zscored;
        traces = traces_zscored;
    end
    
else %df/f
    mean_F = mean(traces,2,'omitnan');
    traces = (traces - mean_F)./mean_F; %(total-mean)/mean
end

%% Check which behavior data is present in the block

TrialData = ExtractedData.TrialData(tableRow); %Traces are already resampled and cropped to be the same # of frames in ExtractedData

[hasLoco, hasXY, hasZcorr, hasResiduals, hasWhisker, hasPupil] = deal(false);
columns = fields(TrialData);
if any(ismember(columns, 'Loco')); hasLoco = ~isempty(TrialData.Loco); end
if any(ismember(columns, 'Full_xoff')); hasXY = ~isempty(TrialData.Full_xoff); end
if any(ismember(columns, 'Full_zcorr')); hasZcorr = ~isempty(TrialData.Full_zcorr); end
if hasLoco && hasXY; hasResiduals = true; end
if any(ismember(columns, 'Full_Whisker')); hasWhisker = ~isempty(TrialData.Full_Whisker); end
if any(ismember(columns, 'Full_Pupil')); hasPupil = ~isempty(TrialData.Full_Pupil); end

timestamp = TrialData.Timestamp;

if strcmp(ops.TraceToAnalyze, 'Loco')
    if ~hasLoco; error('Missing loco data'); end
elseif strcmp(ops.TraceToAnalyze, 'XY')
    if ~hasXY; error('Missing x/y data'); end
elseif strcmp(ops.TraceToAnalyze, 'zcorr')
    if ~hasZcorr; error('Missing zcorr data'); end
elseif strcmp(ops.TraceToAnalyze, 'Residuals')
    if ~hasResiduals; error('Missing loco or x/y data'); end
elseif strcmp(ops.TraceToAnalyze, 'Whisker')
    if ~hasWhisker; error('Missing whisker data'); end
elseif strcmp(ops.TraceToAnalyze, 'Pupil')
    if ~hasPupil; error('Missing pupil data'); end
end

%% Extract and plot locomotor activity & convert to Z-score

if ops.PlotFigure; fig1 = figure; hold on; else; fig1 = []; end
nPanels = sum([1,hasLoco,hasXY,hasZcorr,hasResiduals,hasWhisker,hasPupil]);
PanelCount = 1;
    
%Plot one example trace
if ops.PlotFigure
    subplot(nPanels,1,PanelCount); hold on; PanelCount = PanelCount + 1;
    title('Trace')
    ylabel('Z-score')
    plot(timestamp, traces(1,:), 'LineWidth', 0.25, 'Color', 'k')
    xlim([timestamp(1) timestamp(end)])
end

if hasLoco %PLOT LOCO
    
    %Time-corrected loco trace matched to Bruker data
    loco_speed = TrialData.Full_Loco;
    loco_speed_zscored = zscore(loco_speed,'',2);

    %Plot locomotor activity
    if ops.PlotFigure
        subplot(nPanels,1,PanelCount); hold on; PanelCount = PanelCount + 1;
        title('Locomotor activity')
        ylabel('Speed (cm/s)') 
        plotMax = max(loco_speed) + 1;
        area(timestamp,(loco_speed > 0)*plotMax, 'EdgeColor', 'none', 'Facecolor', [210/255, 248/255, 210/255])
        area(timestamp,(loco_speed == 0)*plotMax, 'EdgeColor', 'none', 'Facecolor', [238/255, 144/255, 144/255])
        plot(timestamp, loco_speed, 'LineWidth', 0.25, 'Color', 'k')
        hline(0.7) %Noise floor of the wheel
        try ylim([0 plotMax]); catch; end
        xlim([timestamp(1) timestamp(end)])
    end
end

if hasXY %PLOT XY
    
    xoff = TrialData.Full_xoff;
    yoff = TrialData.Full_yoff;
    XY = sqrt(xoff.^2 + yoff.^2); %proxy for motion in the FOV: compute euclidean distance from x/y shift
    XY_zscored = zscore(XY,'',2);
    
    if ops.PlotFigure
        subplot(nPanels,1,PanelCount); hold on; PanelCount = PanelCount + 1;
        title('X/Y movement')
        ylabel('pixels') 
        plot(timestamp, xoff, 'LineWidth', 0.25)
        plot(timestamp, yoff, 'Linewidth', 0.25)
        plot(timestamp, XY, 'Linewidth', 0.25, 'Color', 'k')
        xlim([timestamp(1) timestamp(end)])
        xlabel('Seconds')
        legend({'xoff', 'yoff', 'combined'})
    end
end


if hasLoco && hasXY %PLOT LOCO+XY RESIDUALS
    
    %Perform linear regression between loco and XY trace and compute residuals
    p = polyfit(XY_zscored, loco_speed_zscored, 1);
    LocoFit = polyval(p, XY_zscored);
    LocoResiduals = loco_speed_zscored - LocoFit;
    LocoResiduals_zscored = zscore(LocoResiduals,'',2);
    
    if ops.PlotFigure
        subplot(nPanels,1,PanelCount); hold on; PanelCount = PanelCount + 1;
        title('Residuals of Loco + X/Y movement')
        ylabel('Z-scores') 
        plot(timestamp, XY_zscored, 'LineWidth', 0.25, 'Color', [210/255, 248/255, 210/255])
        plot(timestamp, loco_speed_zscored, 'LineWidth', 2, 'Color', [238/255, 144/255, 144/255])
        plot(timestamp, LocoResiduals, 'LineWidth', 0.25, 'Color', 'k')
        xlim([timestamp(1) timestamp(end)])
        xlabel('Seconds')
        legend({'Z-scored XY', 'Z-scored Loco', 'Residuals'})
    end
end

if hasZcorr %PLOT zcorr
    
    zcorr = TrialData.Full_zcorr;
    
    %Determine best z plane and amount of drift
    max_drift = 10; %Number of planes that the imaging plane can drift in +/- Z
    [~, best_z] = max(zcorr,[],1);
    imaging_plane = mode(best_z);
    A = imaging_plane - max_drift;
    B = imaging_plane + max_drift;
    if A < 1
        A = 1;
    end
    %It is necessary to crop otherwise superfluous high correlations with far away planes can happen
    zcorr_cropped = zcorr(A:B,:);
    [~, best_z_cropped] = max(zcorr_cropped,[],1);
    imaging_plane_cropped = mode(best_z_cropped);

    %Determine amount of zdrift (this is what we will correlate with brain activity)
    zdrift = imaging_plane_cropped - best_z_cropped;
    zdrift_zscored = zscore(zdrift,'',2);

    if ops.PlotFigure
        subplot(nPanels,1,PanelCount); hold on; PanelCount = PanelCount + 1;
        title('Z movement')
        ylabel('Z slices') 
        plot(timestamp, zdrift, 'LineWidth', 0.25)
        xlim([timestamp(1) timestamp(end)])
        ylim([-max_drift*2, max_drift*2])
        xlabel('Seconds')
        legend('Z drift')

        subplot(5,1,5); hold on
        imagesc(zcorr_cropped)
        plot(best_z_cropped, '-w', 'LineWidth', 0.05)
        xlabel('Frames')
        ylabel('Z slices')
        xlim([1 size(zcorr_cropped,2)])
        ylim([1 size(zcorr_cropped,1)])
    end
end

if hasWhisker %PLOT XY
    
    whisker = TrialData.Full_Whisker;
    whisker_withoutnans = whisker(~isnan(whisker));
    whisker_withoutnans_zscored = zscore(whisker_withoutnans,'',2);
    whisker_zscored = whisker;
    whisker_zscored(~isnan(whisker)) = whisker_withoutnans_zscored;
    
    if ops.PlotFigure
        subplot(nPanels,1,PanelCount); hold on; PanelCount = PanelCount + 1;
        title('Whisker')
        ylabel('z-scored whisker') 
        plot(timestamp, whisker_zscored, 'LineWidth', 0.25)
        xlim([timestamp(1) timestamp(end)])
        xlabel('Seconds')
    end
end

if hasPupil %PLOT XY
    
    pupil = TrialData.Full_Pupil;
    pupil_withoutnans = pupil(~isnan(pupil));
    pupil_withoutnans_zscored = zscore(pupil_withoutnans,'',2);
    pupil_zscored = pupil;
    pupil_zscored(~isnan(pupil)) = pupil_withoutnans_zscored;
    
    if ops.PlotFigure
        subplot(nPanels,1,PanelCount); hold on; PanelCount = PanelCount + 1;
        title('Pupil')
        ylabel('z-scored pupil') 
        plot(timestamp, pupil_zscored, 'LineWidth', 0.25)
        xlim([timestamp(1) timestamp(end)])
        xlabel('Seconds')
    end
end

if ops.PlotFigure; tak_suptitle(regexprep(TrialData.Block,'_',' ','emptymatch')); end

%% Use tak_compute_xcorr to loop through each cell and save cross-correlations

%Set trace to correlate with neural data
switch ops.TraceToAnalyze
    case 'Loco'
        traceToAnalyze = loco_speed_zscored;
        
    case 'XY'
        traceToAnalyze = XY_zscored;
        
    case 'zcorr'
        traceToAnalyze = zdrift_zscored;
        
    case 'Residuals'
        traceToAnalyze = LocoResiduals_zscored;
        
    case 'Whisker'
        traceToAnalyze = whisker_zscored;
        
    case 'Pupil'
        traceToAnalyze = pupil_zscored;
        
    otherwise
        error('Unexpected TraceToAnalyze')        
end

[cc, lag_in_s, MotorCorrData, ~] = tak_compute_xcorr(traces, traceToAnalyze, MotorCorrOps.maxlag_in_s, MotorCorrOps.fs, 'nShuffles', MotorCorrOps.shuffles, 'pvalue', MotorCorrOps.p_value, 'PlotFigures', ops.PlotXCorrFigures);

end %end function