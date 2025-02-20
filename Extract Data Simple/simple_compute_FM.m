function [FMData, fig1] = simple_compute_FM(isResponsive, PeakData, R, StimInfo, plot_figures)
% Analyze responses to FM sweep
%
% Argument(s): 
%   isResponsive - array of 0s and 1s indicating whether the
%   cell is significantly responsive to that stim combination
%   PeakData - from simple_check_if_responsive
%   R - enter all NaNs size of isResponsive to ignore, otherwise provide R from ReliabilityData or XCorrData
%   StimInfo - from simple_prepare_variables
%   plot_figures - 0 or 1 to plot figures
%
% Returns:
%   FMData table:
%     - BestSpeed = best FM speed according to highest amplitude response
%     - fit_width = bandwidth of speed preference
%     - fit_width_log = bandwidth of speed preference (computed in equal increments)
%     - MI = modulation index (preference for down vs. up sweeps)
%     - MI_p = p value for test against shuffled distribution
%     - MI_p_bs = p value for comparing bootstrapped distribution to 0
%     - MI_reliability = same as MI but using R as the value
%     - MI_binary = same as MI but using 0 and 1 as the value
%     - BestAbsSpeed = best FM speed regardless of direction according to highest amplitude response
%     - SMI = speed mod. index (preference for fast vs. slow sweeps)
%     - SMI_p = p value for test against shuffled distribution
%     - SMI_p_bs = p value for comparing bootstrapped distribution to 0
%     - SMI_reliability = same as SMI but using R as the value
%     - SMI_binary = same as SMI but using 0 and 1 as the value
%     - Slope = slope of responses to absolute speed speed
%     - Slope_p = p value for test against shuffled distribution
%     - Slope_p_bs = p value for comparing bootstrapped distribution to 0
%     - SSI = speed selectivity index (preference for sweep speed)
%     - SSI_p = p value for test against shuffled distribution
%     - SSI_p_bs = p value for comparing bootstrapped distribution to 0
%     - SSI_reliability = same as SSI but using R as the value
%     - SSI_binary = same as SSI but using 0 and 1 as the value
%     - FRA_area = % of RF == 1
%     - Raw_BestSpeed = computed ignoring isRF (best values not restricted to resposive & reliable)
%     - Raw_fit_width_log = " "
%   fig1 - Gaussian fit + Neurometric curve
% 
% Version History:
%  - V1 = current version made for simple_extract_data
%
% TODO: Find best starting values for fitting
% Search 'TODO'
%% Initial parameters

nShuffles = 1000;
alpha = 0.05;
speed_boundary = 50; %oct/s to delineate fast vs. slow speeds

%Figures will return as [] if plot_figures is 0
fig1 = [];

%FMData will return as NaN if there are no isResponsive stim conditions
FMData = table;
FMData.BestSpeed = single(nan);
FMData.fit_width = single(nan);
FMData.fit_width_log = single(nan);
FMData.MI = single(nan);
FMData.MI_p = single(nan);
FMData.MI_p_bs = single(nan);
FMData.MI_reliability = single(nan);
FMData.MI_binary = single(nan);
FMData.BestAbsSpeed = single(nan);
FMData.BestAbsSpeed_fit = single(nan);
FMData.abs_fit_width = single(nan);
FMData.SMI = single(nan);
FMData.SMI_p = single(nan);
FMData.SMI_p_bs = single(nan);
FMData.SMI_reliability = single(nan);
FMData.SMI_binary = single(nan);
FMData.Slope = single(nan);
FMData.Slope_p = single(nan);
FMData.Slope_p_bs = single(nan);
FMData.SSI = single(nan);
FMData.SSI_reliability = single(nan);
FMData.SSI_binary = single(nan);
FMData.FRA_area = single(nan);
FMData.Raw_BestSpeed = single(nan);
FMData.Raw_fit_width_log = single(nan);
FMData.Response = single(nan(1, size(StimInfo.V1_Label,1)));
FMData.Abs_Response = single(nan(1, size(StimInfo.V1_Label,1)/2));

%Check RF for responsive stim conditions and make sure there is >1 stim condition:
if ~any(isResponsive) || isequal(size(isResponsive),[1,1])
    return
end

%% compute FRA area using isResponsive

%Add NaNs to IsRF where there were no trial to analyze (e.g. in Running vs. NotRunning)
type = PeakData.ResponseType;
IsRF_withNans = single(isResponsive);
IsRF_withNans(strcmp(type, "undetermined")) = NaN;

%Compute FRA on non-NaN values only
IsRF_removeNans = IsRF_withNans;
IsRF_removeNans(isnan(IsRF_removeNans)) = [];

FMData.FRA_area = sum(IsRF_removeNans)/numel(IsRF_removeNans);

%% Reshape responses
speeds = StimInfo.Combined(:,1).*StimInfo.Combined(:,2);

%Which PeakData to use for Response value?
%--> Use AUC, but make suppressed AUCs negative for proper inhibitory/excitatory RF detection
%--> If ALL significant responses are suppressed, flip the signs so that we can measure suppressed RFs
RF_type = determine_RF_response_type(type, isResponsive);

%Generate RF from AUC, assuming excitatory or mixed RF for now
response = nan(size(PeakData.Peak_AUC));
for i = 1:length(response)
    response_type = PeakData.ResponseType(i);
    if strcmp(response_type, 'suppressed')
        response(i) = -PeakData.Trough_AUC(i);
    else
        response(i) = PeakData.Peak_AUC(i);
    end
end

%If all significant responses are suppressed, flip the sign of RF
response_orig = response; %Save for storing in FMData without sign flip
if strcmp(RF_type, 'inhibitory')
    response = -response;
end

%Reorganize data from -sweeps to +sweeps
if ~issorted(speeds(speeds > 0),'ascend') || ~issorted(speeds(speeds < 0), 'descend')
    error('Speeds are in unexpected format for reordering')
end
reorder_idx = [flipud(find(speeds < 0)); find(speeds > 0)];
speeds = speeds(reorder_idx);
type = type(reorder_idx);
isResponsive = isResponsive(reorder_idx);
response = response(reorder_idx);
response_orig = response_orig(reorder_idx);
R = R(reorder_idx);
Binary = IsRF_withNans(reorder_idx);

%Get responses for absolute speeds (average together up and down sweeps)
abs_speeds = abs(speeds);
unique_speeds = unique(abs_speeds);

abs_response = nan(size(unique_speeds));
abs_response_orig = nan(size(unique_speeds)); %Save for storing in FMData without sign flip
abs_R = nan(size(unique_speeds));
abs_Binary = nan(size(unique_speeds));
for i = 1:length(abs_response)
    abs_response(i) = mean(response(abs_speeds == unique_speeds(i)), 'omitnan');
    abs_response_orig(i) = mean(response_orig(abs_speeds == unique_speeds(i)), 'omitnan');
    abs_R(i) = mean(R(abs_speeds == unique_speeds(i)), 'omitnan');
    abs_Binary(i) = any(Binary(abs_speeds == unique_speeds(i)));
end

%Store response in FMData before removing values (for plotting)
FMData.Response = single(response_orig');
FMData.Abs_Response = single(abs_response_orig'); %low to high

%% d prime sensitivity index (from Romero, Polley et al Cerebral Cortex 2019)
%This is the only part of the code where we keep the responses in a matrix

%Get RF mapping from StimInfo
RF_Map = StimInfo.RF_Map;
nanmat = nan(max(RF_Map));
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

%Assign RF values to RF map
[RF, thresholded_RF] = deal(nanmat);
RF(RF_ind) = response_orig;
temp_response_thresholded = response_orig;
temp_response_thresholded(~isResponsive) = nan;
thresholded_RF(RF_ind) = temp_response_thresholded;
[M,I] = max(thresholded_RF(:));
[best_row, best_col] = ind2sub(size(RF),I);

%% Eliminate stim that aren't represented (e.g. due to loco or changes in stim protocol) or have NaN response
missing_ind = strcmp(type, 'undetermined');
nan_ind = isnan(response);
remove_ind = (missing_ind + nan_ind) > 0;
isResponsive(remove_ind) = [];
speeds(remove_ind) = [];
abs_speeds(remove_ind) = [];
response(remove_ind) = [];
%type(remove_ind) = [];
R(remove_ind) = [];
Binary(remove_ind) = [];

%Do the same for abs response
unique_speeds(isnan(abs_response)) = [];
abs_R(isnan(abs_response)) = [];
abs_Binary(isnan(abs_response)) = [];
abs_response(isnan(abs_response)) = [];

%% Find preferred speed [Two ways]
% Method #1: Use speed with greatest response

%Raw values
if ~all(isnan(response))
    [~, max_ind] = max(response);
    
    FMData.Raw_BestSpeed = single(speeds(max_ind));
    
    %For fitting
    raw_best_ind = zeros(size(speeds));
    raw_best_ind(max_ind) = 1;
else
    raw_best_ind = [];
end

%Thresholded values
response_thresholded = response;
response_thresholded(~isResponsive) = nan;

if ~all(isnan(response_thresholded))
    [~, max_ind] = max(response_thresholded);
    
    FMData.BestSpeed = single(speeds(max_ind));
    
    %For fitting
    best_ind = zeros(size(speeds));
    best_ind(max_ind) = 1;
else
    best_ind = []; %For fitting
end

%% Method #2: Fit a gaussian across all speeds 

FitTable = tak_fit_gaussian(speeds, response, best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

%Fit table will be empty if there were not at least 3 points to fit a gaussian to
if ~isempty(FitTable) 
    FMData.fit_width = single(FitTable.Backtransformed_Width);
    FMData.fit_width_log = single(FitTable.Width);
end

%Compute width with raw values
Raw_FitTable = tak_fit_gaussian(speeds, response, raw_best_ind, 'LogTransform', 1, 'PlotFigure', 0);
Raw_FitTable = Raw_FitTable(Raw_FitTable.BestFit == 1,:); %Only keep best fit

%Fit table will be empty if there were not at least 3 points to fit a gaussian to
if ~isempty(Raw_FitTable) 
    FMData.Raw_fit_width_log = single(Raw_FitTable.Width);
end

%% Up vs. Down Sweep Modulation Index
% Resp_Up - Resp_Down / Resp_Up + Resp_Down
% Asymmetry index betwen -1 and 1 indicating whether cell responds more to up or down sweeps
% If value is closer to -1, cell prefers down sweeps, if value is closer to  +1, cell prefers up sweeps

up_ind = speeds > 0;
down_ind = speeds < 0;
FMData.MI = tak_compute_asymmetry_index(response, up_ind, down_ind);
FMData.MI_reliability = tak_compute_asymmetry_index(R, up_ind, down_ind);
FMData.MI_binary = tak_compute_asymmetry_index(Binary, up_ind, down_ind);

%% Test Modulation Index (Two ways)

if ~isnan(FMData.MI) %MI will be NaN if there weren't both up and down values
    
    MI_test = struct;

    %% FIRST WAY: Create shuffled distribution by shuffling stim IDs
    %Generate list of shuffled IDs
    N1 = size(speeds,1); % number of stim

    shuffled_responses = nan(nShuffles,N1);
    for u = 1:nShuffles
        permIDs = randperm(N1);

        %Make sure the original order isn't included
        while isequal(permIDs,1:N1)
            permIDs = randperm(N1);
        end

        %Threshold negative values to be 0
        tempResponse = response(permIDs);
        tempResponse(tempResponse < 0) = 0;
        
        shuffled_responses(u,:) = tempResponse;
    end

    %Compute MI on all shuffled IDs
    MI_Shuffled = tak_compute_asymmetry_index(shuffled_responses, up_ind, down_ind);

    %Perform two-tailed z test
    [MI_test.Shuffle.Test, MI_test.Shuffle.P, MI_test.Shuffle.CI, MI_test.Shuffle.R] = ztest(FMData.MI, mean(MI_Shuffled), std(MI_Shuffled), 'Alpha', alpha);
    
    %Store p value
    FMData.MI_p = single(MI_test.Shuffle.P);
    
    %% SECOND WAY: Compare bootstrapped MI distribution to 0

    %Calculate the average of r subsamples nShuffle number of times with replacement
    N1 = size(speeds,1); % number of stim
    
    if N1 > 2
        r = N1-1; %Subsample (keep as much of the data as possible to make sure we have both + and - sweeps)

        shuffle_inds = sort(repmat((1:nShuffles)',[r,1]));
        bootstrap_ind = randsample(1:N1,r*nShuffles,'true')';

        MI_bootstrap = nan(nShuffles,1);
        for u = 1:nShuffles
            %Define sample
            sample_ind = bootstrap_ind(shuffle_inds == u);
            sample_speeds = speeds(sample_ind);

            %Make sure we don't have all downs or all ups
            if all(sample_speeds < 0) || all(sample_speeds > 0)
                hasUpAndDown = 0;
                while ~hasUpAndDown
                    sample_ind = randsample(1:N1,r,'true')';
                    sample_speeds = speeds(sample_ind);
                    if ~(all(sample_speeds < 0) || all(sample_speeds > 0))
                        hasUpAndDown = 1;
                    end
                end
            end

            %Compute MI
            sample_response = response(sample_ind);
            sample_up_ind = sample_speeds > 0;
            sample_down_ind = sample_speeds < 0;
            MI_bootstrap(u) = tak_compute_asymmetry_index(sample_response, sample_up_ind, sample_down_ind);

            if isnan(MI_bootstrap(u))
                if sample_Resp_Up == 0 && sample_Resp_Down == 0
                    %If value is NaN because numerator and denominator are 0, set MI to 0
                    MI_bootstrap(u) = 0;
                else
                    error('Correct NaNs')
                end
            end
        end

        %Perform Z-test asking whether our bootstrapped MI distribution is significantly different from 0
        if FMData.MI > 0; tail = 'left'; elseif FMData.MI < 0; tail = 'right'; else; tail = 'both'; end %Determine whether to do left or right-tailed test
        [MI_test.Bootstrap.Test, MI_test.Bootstrap.P, MI_test.Bootstrap.CI, MI_test.Bootstrap.R] = ztest(0, mean(MI_bootstrap), std(MI_bootstrap), 'Alpha', alpha, 'tail', tail);

        %Store p value
        FMData.MI_p_bs = single(MI_test.Bootstrap.P);
    else
        %Cannot compute
        MI_test.Bootstrap.Test = NaN;
        MI_test.Bootstrap.P = NaN;
        MI_test.Bootstrap.CI = [NaN, NaN];
        MI_test.Bootstrap.R = NaN;
        FMData.MI_p_bs = NaN;
    end
end

%% Speed selectivity (Best Abs Speed and SSI)

%Find best absolute speed based on max response
%Can't threshold response here because we averaged together up/down sweeps, so we don't have responsive & reliability info
[~, ind] = max(abs_response);
FMData.BestAbsSpeed = unique_speeds(ind);

%Compute speed selectivity index (kind of like sparseness)
%How selective is the preferred speed? Speed Selectivity Index:
%1 - (avg response at all other speeds/preferred speed)
%Result will be between 0 and 1: speed selective = 1, not selective = 0
[~, FMData.SSI] = tak_compute_selectivity_index(abs_response, ind, 1, 'Method', 'Sparseness', 'PlotFigure', 0);
[~, FMData.SSI_reliability] = tak_compute_selectivity_index(abs_R, ind, 1, 'Method', 'Sparseness', 'PlotFigure', 0);
[~, FMData.SSI_binary] = tak_compute_selectivity_index(abs_Binary, ind, 1, 'Method', 'Sparseness', 'PlotFigure', 0);

%Fit best absolute speed using a gaussian
AbsFitTable = tak_fit_gaussian(unique_speeds, abs_response, [], 'LogTransform', 1, 'PlotFigure', 0);
AbsFitTable = AbsFitTable(AbsFitTable.BestFit == 1,:); %Only keep best fit

%Fit table will be empty if there were not at least 3 points to fit a gaussian to
if ~isempty(AbsFitTable) 
    FMData.BestAbsSpeed_fit = single(AbsFitTable.Backtransformed_Best_X);
    FMData.abs_fit_width = single(AbsFitTable.Backtransformed_Width);
end

%% Speed Modulation Index (SMI)
% Resp_Fast - Resp_Slow / Resp_Fast + Resp_Slow
% Asymmetry index betwen -1 and 1 indicating whether cell responds more to  slow or fast sweeps
% If value is closer to -1, cell prefers slow sweeps, if value is closer to  +1, cell prefers fast sweeps

fast_ind = unique_speeds >= speed_boundary;
slow_ind = unique_speeds < speed_boundary; 
FMData.SMI = tak_compute_asymmetry_index(abs_response, fast_ind, slow_ind);
FMData.SMI_reliability = tak_compute_asymmetry_index(abs_R, fast_ind, slow_ind);
FMData.SMI_binary = tak_compute_asymmetry_index(abs_Binary, fast_ind, slow_ind);

%% Test Speed Modulation Index (Two ways)

if ~isnan(FMData.SMI) && length(unique_speeds) >= 3 %SMI will be nan if there weren't both fast and slow speeds
    %Need at least 3 speeds to compute SMI bootstrap
    
    SMI_test = struct;
    
    %FIRST WAY: Create shuffled distribution by shuffling stim IDs
    %Generate list of shuffled IDs
    N1 = size(unique_speeds,1); % number of stim

    shuffled_responses = nan(nShuffles,N1);
    for u = 1:nShuffles
        permIDs = randperm(N1);

        %Make sure the original order isn't included
        while isequal(permIDs,1:N1)
            permIDs = randperm(N1);
        end
        
        %Threshold negative values to be 0
        tempResponse = abs_response(permIDs);
        tempResponse(tempResponse < 0) = 0;
    
        shuffled_responses(u,:) = tempResponse;
    end

    %Compute SMI on all shuffled IDs
    SMI_Shuffled = tak_compute_asymmetry_index(shuffled_responses, fast_ind, slow_ind);

    % perform two-tailed z-test
    [SMI_test.Shuffle.Test, SMI_test.Shuffle.P, SMI_test.Shuffle.CI, SMI_test.Shuffle.R] = ztest(FMData.SMI, mean(SMI_Shuffled), std(SMI_Shuffled), 'Alpha', alpha);

    %Store p value
    FMData.SMI_p = single(SMI_test.Shuffle.P);
    
    %% SECOND WAY: Compare bootstrapped MI distribution to 0

    %Calculate the average of r subsamples nShuffle number of times with replacement
    N1 = size(unique_speeds,1); % number of stim
    r = N1-1; %Subsample (keep as much of the data as possible to make sure we have both + and - sweeps)

    shuffle_inds = sort(repmat((1:nShuffles)',[r,1]));
    bootstrap_ind = randsample(1:N1,r*nShuffles,'true')';

    SMI_bootstrap = nan(nShuffles,1);
    for u = 1:nShuffles
        %Define sample
        sample_ind = bootstrap_ind(shuffle_inds == u);
        sample_speeds = unique_speeds(sample_ind);

        %Make sure we don't have all fast or all slow
        if all(sample_speeds < speed_boundary) || all(sample_speeds >= speed_boundary)
            hasFastAndSlow = 0;
            while ~hasFastAndSlow
                sample_ind = randsample(1:N1,r,'true')';
                sample_speeds = unique_speeds(sample_ind);
                if ~(all(sample_speeds < speed_boundary) || all(sample_speeds >= speed_boundary))
                    hasFastAndSlow = 1;
                end
            end
        end

        %Compute SMI
        sample_response = abs_response(sample_ind);
        sample_fast_ind = sample_speeds >= speed_boundary;
        sample_slow_ind = sample_speeds < speed_boundary;
        
        SMI_bootstrap(u) = tak_compute_asymmetry_index(sample_response, sample_fast_ind, sample_slow_ind);

        if isnan(SMI_bootstrap(u))
            if sample_Resp_Fast == 0 && sample_Resp_Slow == 0
                %If value is NaN because numerator and denominator are 0, set MI to 0
                SMI_bootstrap(u) = 0;
            else
                error('Correct NaNs')
            end
        end
    end

    %Perform Z-test asking whether our bootstrapped SMI distribution is significantly different from 0
    if FMData.SMI > 0; tail = 'left'; elseif FMData.SMI < 0; tail = 'right'; else; tail = 'both'; end %Determine whether to do left or right-tailed test
    [SMI_test.Bootstrap.Test, SMI_test.Bootstrap.P, SMI_test.Bootstrap.CI, SMI_test.Bootstrap.R] = ztest(0, mean(SMI_bootstrap), std(SMI_bootstrap), 'Alpha', alpha, 'tail', tail);

    %Store p value
    FMData.SMI_p_bs = single(SMI_test.Bootstrap.P);
end

%% Slope of the responses

if length(unique_speeds) >= 2
    x = 1:length(unique_speeds);
    y = abs_response;
    coefficients = polyfit(x, y, 1); % Get coefficients of a line fit through the data.
    xFit = linspace(min(x), max(x), 1000); % Create a new x axis with exactly 1000 points (or whatever you want).
    yFit = polyval(coefficients , xFit); % Get the estimated yFit value for each of those 1000 new x locations.
    FMData.Slope = single(coefficients(1));
end

%% Test Slope (Two ways)

if ~isnan(FMData.Slope) %Slope will be nan if there weren't at least 2 points
    
    Slope_test = struct;
    
    %FIRST WAY: Create shuffled distribution by shuffling stim IDs
    %Generate list of shuffled IDs
    N1 = size(unique_speeds,1); % number of stim
    slope_shuffled = nan(nShuffles,1);
    for u = 1:nShuffles
        permIDs = randperm(N1);

        %Make sure the original order isn't included
        while isequal(permIDs,1:N1)
            permIDs = randperm(N1);
        end

        y = response(permIDs);

        %Compute slope
        coefficients = polyfit(x, y, 1); % Get coefficients of a line fit through the data.
        slope_shuffled(u) = coefficients(1);
    end

    %Perform two-tailed z test
    [Slope_test.Shuffle.Test, Slope_test.Shuffle.P, Slope_test.Shuffle.CI, Slope_test.Shuffle.R] = ztest(FMData.Slope, mean(slope_shuffled), std(slope_shuffled), 'Alpha', alpha);

    %Store p value
    FMData.Slope_p = single(Slope_test.Shuffle.P);

    %% SECOND WAY: Compare bootstrapped slope distribution to 0

    N1 = size(unique_speeds,1); % number of stim
    r = 2;%2 is minimum number of points needed to calculate slope

    %Calculate the average of r subsamples nShuffle number of times WITHOUT replacement
    slope_bootstrap = nan(nShuffles,1);
    for u = 1:nShuffles
        %Define sample
        sample_ind = sort(randsample(1:N1,r,'false')); %sort randsample so that speeds will be in order
        %sample_speeds = unique_speeds(sample_ind);
        sample_responses = abs_response(sample_ind);

        %Compute slope
        x_bootstrap = 1:length(sample_responses);
        y_bootstrap = sample_responses;
        coefficients = polyfit(x_bootstrap, y_bootstrap, 1); % Get coefficients of a line fit through the data.
        slope_bootstrap(u) = coefficients(1);
    end

    %Catch issues
    if any(isnan(slope_bootstrap))
        error('Figure out why there are NaNs')
    end

    %Perform Z-test asking whether our bootstrapped slope distribution is significantly different from 0
    if FMData.Slope > 0; tail = 'left'; elseif FMData.Slope < 0; tail = 'right'; else; tail = 'both'; end %Determine whether to do left or right-tailed test
    [Slope_test.Bootstrap.Test, Slope_test.Bootstrap.P, Slope_test.Bootstrap.CI, Slope_test.Bootstrap.R] = ztest(0, mean(slope_bootstrap), std(slope_bootstrap), 'Alpha', alpha, 'tail', tail);    

    %Store p value
    FMData.Slope_p_bs = single(Slope_test.Bootstrap.P);

end

%% Plot figure

if plot_figures
    fig1 = figure;
    
    %Plot fit
    subplot(3,3,1); hold on
    if isempty(FitTable)
        plot(speeds, response, 'Linewidth', 2);
        xlim([speeds(1) speeds(end)])
        set(gca, 'XTick', speeds)
        set(gca, 'XTickLabel', num2str(speeds))
    else
        tak_plot_gaussian(FitTable)
    end
    xlabel('Octaves per second')
    ylabel('Response');
    title(['Best Speed: ' num2str(FMData.BestSpeed)])
        
    %Plot MI tests    
    if ~isnan(FMData.MI)
        subplot(3,3,4); hold on
        histogram(MI_Shuffled, 'facecolor', 'g', 'facealpha', 0.2)
        vline(FMData.MI)
        ylabel('Shuffles')
        title(['MI (up vs. down): ' num2str(FMData.MI)])
        if MI_test.Shuffle.Test; textcolour = ' {\color{green}p=} '; else; textcolour = ' {\color{black}p=} '; end
        subtitle(strcat('Z=', sprintf('%.2f', MI_test.Shuffle.R), textcolour, sprintf('%.2f', MI_test.Shuffle.P)))
    
        subplot(3,3,7); hold on
        histogram(MI_bootstrap, 'facecolor', 'b', 'facealpha', 0.2)
        vline(FMData.MI)
        vline(0, 'k')
        ylabel('Bootstraps')
        xlabel('Modulation Index')
        if MI_test.Bootstrap.Test; textcolour = ' {\color{green}p=} '; else; textcolour = ' {\color{black}p=} '; end
        subtitle(strcat('Z=', sprintf('%.2f', MI_test.Bootstrap.R), textcolour, sprintf('%.2f', MI_test.Bootstrap.P)))
    else
        subplot(3,3,4);
        title('MI could not be computed')
    end
    
    %Plot abs speed, fit
    subplot(3,3,2); hold on
    
    if isempty(AbsFitTable)
        plot(unique_speeds, abs_response, 'Linewidth', 2)
        xlim([unique_speeds(1) unique_speeds(end)])
        set(gca, 'XTick', unique_speeds)
        set(gca, 'XTickLabel', num2str(unique_speeds))
    else
        tak_plot_gaussian(AbsFitTable)
    end
    xlabel('Absolute Speeds')
    ylabel('Response')
    title(['Best Abs Speed = ' num2str(FMData.BestAbsSpeed) '  Fit = ' num2str(round(FMData.BestAbsSpeed_fit,2))])
    
    %Plot abs speed, slope
    subplot(3,3,3); hold on
    plot(abs_response, 'Linewidth', 2)
    if ~isnan(FMData.Slope)
        plot(xFit, yFit, 'r-'); % Plot fitted line.
    end
    
    %Plot points for each response
    for i = 1:length(unique_speeds)
        points = response(abs_speeds == unique_speeds(i));
        scatter([i, i], points, 25, 'MarkerEdgeColor', 'k');
    end
    xlim([0.5 length(unique_speeds) + 0.5])
    set(gca, 'XTick', 1:length(unique_speeds))
    set(gca, 'XTickLabel', num2str(unique_speeds))
    xlabel('Absolute Speeds')
    ylabel('AUC')
    title(['Slope = ' num2str(round(FMData.Slope,2)) '  SSI = ' num2str(round(FMData.SSI,2))])
    
    %Plot SMI tests    
    if ~isnan(FMData.SMI)
        subplot(3,3,5); hold on
        histogram(SMI_Shuffled, 'facecolor', 'g', 'facealpha', 0.2)
        vline(FMData.SMI)
        ylabel('Shuffles')
        title(['SMI (fast vs. slow): ' num2str(FMData.SMI)])
        if isnan(SMI_test.Shuffle.Test)
            textcolour = ' {\color{black}p=} ';
        elseif SMI_test.Shuffle.Test
            textcolour = ' {\color{green}p=} ';
        else
            textcolour = ' {\color{black}p=} ';
        end
        subtitle(strcat('Z=', sprintf('%.2f', SMI_test.Shuffle.R), textcolour, sprintf('%.2f', SMI_test.Shuffle.P)));

        subplot(3,3,8); hold on
        histogram(SMI_bootstrap, 'facecolor', 'b', 'facealpha', 0.2)
        vline(FMData.SMI)
        vline(0, 'k')
        ylabel('Bootstraps')
        xlabel('Modulation Index')
        if isnan(SMI_test.Bootstrap.Test)
            textcolour = ' {\color{black}p=} ';
        elseif SMI_test.Bootstrap.Test
            textcolour = ' {\color{green}p=} ';
        else
            textcolour = ' {\color{black}p=} ';
        end
        subtitle(strcat('Z=', sprintf('%.2f', SMI_test.Bootstrap.R), textcolour, sprintf('%.2f', SMI_test.Bootstrap.P)))
    else
        subplot(3,3,5);
        title('SMI could not be computed')
    end
    
    %Plot slope tests
    if ~isnan(FMData.Slope)
        subplot(3,3,6); hold on
        histogram(slope_shuffled, 'facecolor', 'g', 'facealpha', 0.2)
        vline(FMData.Slope)
        ylabel('Shuffles')
        title(['Slope: ' num2str(round(FMData.Slope,2))])
        if Slope_test.Shuffle.Test; textcolour = ' {\color{green}p=} '; else; textcolour = ' {\color{black}p=} '; end
        subtitle(strcat('Z=', sprintf('%.2f', Slope_test.Shuffle.R), textcolour, sprintf('%.2f', Slope_test.Shuffle.P)))

        subplot(3,3,9); hold on
        histogram(slope_bootstrap, 'facecolor', 'b', 'facealpha', 0.2)
        vline(FMData.Slope)
        vline(0, 'k')
        ylabel('Bootstraps')
        xlabel('Slope')
        if Slope_test.Bootstrap.Test; textcolour = ' {\color{green}p=} '; else; textcolour = ' {\color{black}p=} '; end
        subtitle(strcat('Z=', sprintf('%.2f', Slope_test.Bootstrap.R), textcolour, sprintf('%.2f', Slope_test.Bootstrap.P)))
    else
        subplot(3,3,6);
        title('Slope could not be computed')
    end
end
