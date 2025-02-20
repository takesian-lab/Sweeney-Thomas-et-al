function [SAMFreqData, fig1] = simple_compute_SAMfreq(IsResponsive, PeakData, StimInfo, plot_figures)
% Analyze responses to SAM data
%
% Argument(s): 
%   IsResponsive - 0s and 1s indicating whether the cell is significantly responsive to that stim combination
%   PeakData - from simple_check_if_resopnsive
%   StimInfo - from simple_prepare_variables
%   plot_figures - 0 or 1 to plot figures
%
% Returns:
%   SAMFreqData table:
%     - BF = freq (kHz) with highest amplitude response
%     - BF_fit_halfwidth = bandwidth of freq preference (in octaves)
%     - BestDepth = SAM depth (%) with highest amplitude response (corresponds to best rate)
%     - BestDepth_fit_halfwidth = bandwidth of depth preference
%     - BestDepth_fit_halfwidth_log = bandwidth of depth preference (computed in equal increments)
%     - DSI = Depth Selectivity Index (adapted from Mesik et al ISI)
%     - dprime = sensitivity index (adapted from Romero et al)
%     - Raw_dprime  =  " " (best value not restricted to responsive and reliable)
%     - Binary_dprime = " " computed on zeros and ones
%     - Raw_BF = computed ignoring isRF (best values not restricted to resposive & reliable)
%     - Raw_BestDepth = computed ignoring isRF (best values not restricted to resposive & reliable)
%     - Raw_BF_fit_halfwidth = " "
%     - Raw_BestDepth_fit_halfwidth = " "
%     - Raw_BestDepth_fit_halfwidth_log = " " computed on scale with equal increments
%     - FRA_type = simple or complex (multipeaked)
%     - FRA_area = % of RF deemed to be responsive

% Version History:
%  - V1 = current version made for simple_extract_data
%
% TODO: 
% Search 'TODO'
%% Initial parameters

%Get RF mapping from StimInfo
RF_Map = StimInfo.RF_Map;
nanmat = nan(max(RF_Map));
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

%Reconstruct RF maps
[IsRF, RF] = deal(nanmat);
IsRF(RF_ind) = IsResponsive;

%Generate RF from AUC, assuming excitatory or mixed RF for now
temp_RF = nan(size(IsResponsive));
for i = 1:length(temp_RF)
    response_type = PeakData.ResponseType(i);
    if strcmp(response_type, 'suppressed')
        temp_RF(i) = -PeakData.Trough_AUC(i);
    else
        temp_RF(i) = PeakData.Peak_AUC(i);
    end
end

%Assign RF values to RF map
RF(RF_ind) = temp_RF;
RF_orig = RF; %Save for plotting

%If all significant responses are suppressed, flip the sign of RF
temp_RF_orig = temp_RF; %Save for storing in SAMFreqData without sign flip
RF_type = determine_RF_response_type(PeakData.ResponseType, IsResponsive);
if strcmp(RF_type, 'inhibitory')
    temp_RF = -temp_RF;
    RF = -RF;
end

%Figures will return as [] if plot_figures is 0
fig1 = [];

%Variables will return as NaN if there are no isResponsive stim conditions
SAMFreqData = table;
SAMFreqData.BF = single(nan);
SAMFreqData.BF_fit_halfwidth = single(nan);
SAMFreqData.BestDepth = single(nan);
SAMFreqData.BestDepth_fit_halfwidth = single(nan);
SAMFreqData.BestDepth_fit_halfwidth_log = single(nan);
SAMFreqData.DSI = single(nan);
SAMFreqData.dprime = single(nan);
SAMFreqData.Raw_dprime = single(nan);
SAMFreqData.Binary_dprime = single(nan);
SAMFreqData.Raw_BF = single(nan);
SAMFreqData.Raw_BestDepth = single(nan);
SAMFreqData.Raw_BF_fit_halfwidth = single(nan);
SAMFreqData.Raw_BestDepth_fit_halfwidth = single(nan);
SAMFreqData.Raw_BestDepth_fit_halfwidth_log = single(nan);
SAMFreqData.FRA_type = "none"; %Double quotes = string
SAMFreqData.FRA_area = single(nan);

%Store response in SAMFreqData (for plotting)
SAMFreqData.Response = single(temp_RF_orig');

%Check RF for responsive stim conditions and make sure there is >1 stim condition:
if ~any(IsResponsive) || isequal(size(IsResponsive),[1,1])
    return
end

%% FRA - Frequency Receptive Area
%make frequency receptive area (FRA) contour
%determine type, single peak if continuous, dual peak, or complex (multi-peak)   

type = strings(size(nanmat));
type(RF_ind) = PeakData.ResponseType;

act_RF = ismember(type, 'activated');
pro_RF = ismember(type, 'prolonged');
inh_RF = ismember(type, 'suppressed');
exc_RF = act_RF + pro_RF; %Combine activated and prolonged

%determine those stimuli that are excitatory or inhibitory and responsive/reliable 
excR_RF = (exc_RF + IsRF) == 2; 
inhR_RF = (inh_RF + IsRF) == 2;

if strcmp(RF_type, 'inhibitory')
    [FRA,L] = bwboundaries(inhR_RF, 8, 'noholes');
else
    [FRA,L] = bwboundaries(excR_RF, 8, 'noholes');
end

if length(FRA) == 1 %Number of objects found by bwboundaries
    FRA_type = 'single peak';
else
    FRA_type = 'complex peak';
end

%% Find best rate and depth [Two ways]
%Method #1: Use stimulus with the greatest response

freq = StimInfo.V1; %use for labelling graphs
freq_log = log2(freq); %use log scale for fitting 
depths = StimInfo.V2*100; %Turn into true percents

%Raw values
[M,I] = max(RF(:));
[raw_best_row, raw_best_col] = ind2sub(size(RF),I);
Raw_BF = freq(raw_best_col);
Raw_BestDepth = depths(raw_best_row);

%Thresholded values
thresholded_RF = RF;
thresholded_RF(~IsRF) = nan;
[M,I] = max(thresholded_RF(:));
[best_row, best_col] = ind2sub(size(RF),I);
BF = freq(best_col);
BestDepth = depths(best_row);

%% compute FRA area using IsRF

%Add NaNs to IsRF where there were no trial to analyze (e.g. in Running vs. NotRunning)
IsRF_withNans = IsRF;
IsRF_withNans(isnan(RF)) = NaN;

%Compute FRA on non-NaN values only
IsRF_removeNans = IsRF_withNans(:);
IsRF_removeNans(isnan(IsRF_removeNans)) = [];

FRA_area = sum(IsRF_removeNans)/numel(IsRF_removeNans);

%% d prime sensitivity index (from Romero, Polley et al Cerebral Cortex 2019)

dprime = tak_compute_RF_dprime(RF, best_row, best_col);
Raw_dprime = tak_compute_RF_dprime(RF, raw_best_row, raw_best_col);
Binary_dprime = tak_compute_RF_dprime(~isnan(thresholded_RF), best_row, best_col);

%% Plot RF and FRA

if plot_figures
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    
    subplot(2,2,1); %Plot RF for reference
    RFmask = ones(size(nanmat))*0.3;
    RFmask(IsRF == 1) = 1;
    imagesc(RF_orig, 'AlphaData', RFmask);
    colorbar
    xlabel('Freq (kHz)')
    ylabel('Depth (%)')
    set(gca,'XTick',1:length(freq))
    set(gca,'YTick',1:length(depths))
    set(gca,'XTickLabel',num2str(freq))
    set(gca,'YTickLabel',num2str(depths))
    title('SAM receptive field')
    subtitle(['dPrime =  ' num2str(round(dprime,1))]);
    text(best_col - 0.25, best_row, 'Best', 'Color', 'w') %-0.25 centers the text
    
    colormap(gca, bluewhitered(256))
    c = colorbar;
    c.Title.String = 'Z';
    
    subplot(2,2,3) %plot FRA contour
    imagesc(RF_orig, 'AlphaData', RFmask);
    colorbar; hold on
    for k = 1:length(FRA)
        boundary = FRA{k};
        if size(boundary,1) <= 2
            scatter(boundary(:,2), boundary(:,1), 'w', 'filled')
        else
            plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth',2);
        end
    end
    xlabel('Freq (kHz)')
    ylabel('Depth (%)')
    set(gca,'XTick',1:length(freq))
    set(gca,'YTick',1:length(depths))
    set(gca,'XTickLabel',num2str(freq))
    set(gca,'YTickLabel',num2str(depths))
    title('Frequency Receptive Area');
    subtitle(strcat({'Type =  '}, FRA_type, {' '}, RF_type));
    
    colormap(gca, bluewhitered(256))
    c = colorbar;
    c.Title.String = 'Z';
end

%% Compute depth selectivity (DSI)
%DSI = Depth Selectivity Index (1 = responds to 0% depth, 0 = responds to 100% depth)

[~, DSI] = tak_compute_selectivity_index(RF, best_row, best_col, 'PlotFigure', 0);

%% Fit #1 - Freq tuning
% determine freq tuning using peak of gaussian fit across freqs at best depth

best_ind = zeros(size(freq));
best_ind(best_col) = 1;
FitTable = tak_fit_gaussian(freq, RF(best_row,:), best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    fit_halfwidth = FitTable.Width; %Keep in octaves
else
    fit_halfwidth = nan;
end

if plot_figures
    subplot(2,2,2); hold on %plot gaussian fits
    if isempty(FitTable)
        plot(freq_log, RF(best_row,:), 'Linewidth', 2);
        set(gca, 'XTick', freq_log)
        set(gca, 'XTickLabel', num2str(freq))
        xlim([freq_log(1) freq_log(end)])
    else
        tak_plot_gaussian(FitTable)
    end
    title(['Best Freq at ' num2str(BestDepth) '% = ' num2str(BF) 'kHz']);
    xlabel('Freq (kHz) log scale');
    ylabel('Response');    
end

% Compute width with raw values
best_ind = zeros(size(freq));
best_ind(raw_best_col) = 1;
FitTable = tak_fit_gaussian(freq, RF(raw_best_row,:), best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Raw_fit_halfwidth = FitTable.Width; %Keep in octaves
else
    Raw_fit_halfwidth = nan;
end

%% Fit #2 - Depth tuning
% determine depth tuning using peak of gaussian fit across depths at best rate

best_ind = zeros(size(depths));
best_ind(best_row) = 1;
FitTable = tak_fit_gaussian(depths, RF(:,best_col), best_ind, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Depth_halfwidth = FitTable.Width;
else
    Depth_halfwidth = nan;
end

%Do it again but this time set x axis to equal increments
FitTable = tak_fit_gaussian(depths, RF(:,best_col), best_ind, 'LogTransform', 2, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Depth_halfwidth_log = FitTable.Width;
else
    Depth_halfwidth_log = nan;
end

if plot_figures
    subplot(2,2,4); hold on %plot gaussian fits
    if isempty(FitTable)
        plot(depths, RF(:,best_col), 'Linewidth', 2);
        set(gca, 'XTick', sort(depths,'ascend'))
        set(gca, 'XTickLabel', num2str(sort(depths,'ascend')))
        xlim([depths(end) depths(1)])
    else
        tak_plot_gaussian(FitTable)
    end
    title(['Best Depth at ' num2str(BF) 'kHz = ' num2str(BestDepth) '%']);
    subtitle(['DSI = ' num2str(round(DSI,1))])
    xlabel('Depth (%)');
    ylabel('Response');    
end

%Compute width with raw values
best_ind = zeros(size(depths));
best_ind(raw_best_row) = 1;
FitTable = tak_fit_gaussian(depths, RF(:,raw_best_col), best_ind, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Raw_Depth_halfwidth = FitTable.Width;
else
    Raw_Depth_halfwidth = nan;
end

%Compute again with log scale
FitTable = tak_fit_gaussian(depths, RF(:,raw_best_col), best_ind, 'LogTransform', 2, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Raw_Depth_halfwidth_log = FitTable.Width;
else
    Raw_Depth_halfwidth_log = nan;
end

%% Store Data
SAMFreqData.BF = single(BF);
SAMFreqData.BF_fit_halfwidth = single(fit_halfwidth);
SAMFreqData.BestDepth = single(BestDepth/100);
SAMFreqData.BestDepth_fit_halfwidth = single(Depth_halfwidth/100);
SAMFreqData.BestDepth_fit_halfwidth_log = single(Depth_halfwidth_log);
SAMFreqData.DSI = single(DSI);
SAMFreqData.dprime = single(dprime);
SAMFreqData.Raw_dprime = single(Raw_dprime);
SAMFreqData.Binary_dprime = single(Binary_dprime);
SAMFreqData.Raw_BF = Raw_BF;
SAMFreqData.Raw_BestDepth = Raw_BestDepth;
SAMFreqData.Raw_BF_fit_halfwidth = single(Raw_fit_halfwidth);
SAMFreqData.Raw_BestDepth_fit_halfwidth = single(Raw_Depth_halfwidth);
SAMFreqData.Raw_BestDepth_fit_halfwidth_log = single(Raw_Depth_halfwidth_log);
SAMFreqData.FRA_type = string(FRA_type);
SAMFreqData.FRA_area = single(FRA_area);

end %end function
