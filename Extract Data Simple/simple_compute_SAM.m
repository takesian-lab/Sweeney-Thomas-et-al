function [SAMData, fig1] = simple_compute_SAM(IsResponsive, PeakData, StimInfo, plot_figures)
% Analyze responses to SAM data
%
% Argument(s): 
%   IsResponsive - 0s and 1s indicating whether the cell is significantly responsive to that stim combination
%   PeakData - from simple_check_if_resopnsive
%   StimInfo - from simple_prepare_variables
%   plot_figures - 0 or 1 to plot figures
%
% Returns:
%   SAMData table:
%     - BestRate = SAM rate (Hz) with highest amplitude response
%     - BestRate_fit_halfwidth = bandwidth of rate preference (computed in log scale, backtransformed)
%     - BestRate_fit_halfwidth_log = bandwidth of rate preference (computed in log scale)
%     - BestDepth = SAM depth (%) with highest amplitude response (corresponds to best rate)
%     - BestDepth_fit_halfwidth = bandwidth of depth preference
%     - BestDepth_fit_halfwidth_log = bandwidth of depth preference (computed in equal increments)
%     - RSI = Rate Selectivity Index (adapted from Mesik et al ISI)
%     - DSI = Depth Selectivity Index (adapted from Mesik et al ISI)
%     - dprime = sensitivity index (adapted from Romero et al)
%     - Raw_dprime  =  " " (best value not restricted to responsive and reliable)
%     - Binary_dprime = " " computed on zeros and ones
%     - Raw_BestRate = computed ignoring isRF (best values not restricted to resposive & reliable)
%     - Raw_BestDepth = computed ignoring isRF (best values not restricted to resposive & reliable)
%     - Raw_BestRate_fit_halfwidth_log = " "
%     - Raw_BestDepth_fit_halfwidth = " "
%     - Raw_BestDepth_fit_halfwidth_log = " "
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
temp_RF_orig = temp_RF; %Save for storing in SAMData without sign flip
RF_type = determine_RF_response_type(PeakData.ResponseType, IsResponsive);
if strcmp(RF_type, 'inhibitory')
    temp_RF = -temp_RF;
    RF = -RF;
end

%Figures will return as [] if plot_figures is 0
fig1 = [];

%Variables will return as NaN if there are no isResponsive stim conditions
SAMData = table;
SAMData.BestRate = single(nan);
SAMData.BestRate_fit_halfwidth = single(nan);
SAMData.BestRate_fit_halfwidth_log = single(nan);
SAMData.BestDepth = single(nan);
SAMData.BestDepth_fit_halfwidth = single(nan);
SAMData.BestDepth_fit_halfwidth_log = single(nan);
SAMData.RSI = single(nan);
SAMData.DSI = single(nan);
SAMData.dprime = single(nan);
SAMData.Raw_dprime = single(nan);
SAMData.Binary_dprime = single(nan);
SAMData.Raw_BestRate = single(nan);
SAMData.Raw_BestDepth = single(nan);
SAMData.Raw_BestRate_fit_halfwidth_log = single(nan);
SAMData.Raw_BestDepth_fit_halfwidth = single(nan);
SAMData.Raw_BestDepth_fit_halfwidth_log = single(nan);
SAMData.FRA_type = "none"; %Double quotes = string
SAMData.FRA_area = single(nan);

%Store response in SAMData (for plotting)
SAMData.Response = single(temp_RF_orig');

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

rates = StimInfo.V1; %use for labelling graphs
depths = StimInfo.V2; %Turn into true percents

%Raw values
[M,I] = max(RF(:));
[raw_best_row, raw_best_col] = ind2sub(size(RF),I);
Raw_BestRate = rates(raw_best_col);
Raw_BestDepth = depths(raw_best_row);

%Thresholded values
thresholded_RF = RF;
thresholded_RF(~IsRF) = nan;
[M,I] = max(thresholded_RF(:));
[best_row, best_col] = ind2sub(size(RF),I);
BestRate = rates(best_col);
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
    xlabel('Rate (Hz)')
    ylabel('Depth (%)')
    set(gca,'XTick',1:length(rates))
    set(gca,'YTick',1:length(depths))
    set(gca,'XTickLabel',num2str(rates))
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
    xlabel('Rate (Hz)')
    ylabel('Depth (%)')
    set(gca,'XTick',1:length(rates))
    set(gca,'YTick',1:length(depths))
    set(gca,'XTickLabel',num2str(rates))
    set(gca,'YTickLabel',num2str(depths))
    title('Frequency Receptive Area');
    subtitle(strcat({'Type =  '}, FRA_type, {' '}, RF_type));
    
    colormap(gca, bluewhitered(256))
    c = colorbar;
    c.Title.String = 'Z';
end

%% Compute rate and depth selectivity (RSI and DSI)
%DSI = Depth Selectivity Index (1 = responds to 0% depth, 0 = responds to 100% depth)
%RSI = Rate Selctivity Index (1 = responds to 0Hz, 0 = responds to 128Hz)

[RSI, DSI] = tak_compute_selectivity_index(RF, best_row, best_col, 'PlotFigure', 0);

%% Fit #1 - Rate tuning
% determine rate tuning using peak of gaussian fit across rates at best depth

best_ind = zeros(size(rates));
best_ind(best_col) = 1;
FitTable = tak_fit_gaussian(rates, RF(best_row,:), best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Rate_halfwidth = FitTable.Backtransformed_Width;
    Rate_halfwidth_log = FitTable.Width;
else
    Rate_halfwidth = nan;
    Rate_halfwidth_log = nan;
end

if plot_figures
    subplot(2,2,2); hold on %plot gaussian fits
    
    if isempty(FitTable)
        plot(rates, RF(best_row,:), 'Linewidth', 2);
        xlim([rates(1) rates(end)])
        set(gca, 'XTick', rates)
        set(gca, 'XTickLabel', num2str(rates))
    else
        tak_plot_gaussian(FitTable)
    end
 
    title(['Best Rate at ' num2str(BestDepth) '% = ' num2str(BestRate) 'Hz']);
    subtitle(['RSI = ' num2str(round(RSI,1))])
    xlabel('Rate (Hz) log scale');
    ylabel('Response');
end

% Compute width with raw values
best_ind = zeros(size(rates));
best_ind(raw_best_col) = 1;
FitTable = tak_fit_gaussian(rates, RF(raw_best_row,:), best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Raw_Rate_halfwidth_log = FitTable.Width;
else
    Raw_Rate_halfwidth_log = nan;
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
        xlim([min(depths) max(depths)])
        set(gca, 'XTick', sort(depths,'ascend'))
        set(gca, 'XTickLabel', num2str(sort(depths, 'ascend')))
        vline(BestDepth_fit)
    else
        tak_plot_gaussian(FitTable)
    end
 
    title(['Best Depth at ' num2str(BestRate) 'Hz = ' num2str(BestDepth) '%']);
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

%Do it again but this time set x axis to equal increments
FitTable = tak_fit_gaussian(depths, RF(:,raw_best_col), best_ind, 'LogTransform', 2, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Raw_Depth_halfwidth_log = FitTable.Width;
else
    Raw_Depth_halfwidth_log = nan;
end

%% Store Data
SAMData.BestRate = single(BestRate);
SAMData.BestRate_fit_halfwidth = single(Rate_halfwidth);
SAMData.BestRate_fit_halfwidth_log = single(Rate_halfwidth_log);
SAMData.BestDepth = single(BestDepth);
SAMData.BestDepth_fit_halfwidth = single(Depth_halfwidth);
SAMData.BestDepth_fit_halfwidth_log = single(Depth_halfwidth_log);
SAMData.RSI = single(RSI);
SAMData.DSI = single(DSI);
SAMData.dprime = single(dprime);
SAMData.Raw_dprime = single(Raw_dprime);
SAMData.Binary_dprime = single(Binary_dprime);
SAMData.Raw_BestRate = single(Raw_BestRate);
SAMData.Raw_BestDepth = single(Raw_BestDepth);
SAMData.Raw_BestRate_fit_halfwidth_log = single(Raw_Rate_halfwidth_log);
SAMData.Raw_BestDepth_fit_halfwidth = single(Raw_Depth_halfwidth);
SAMData.Raw_BestDepth_fit_halfwidth_log = single(Raw_Depth_halfwidth_log);
SAMData.FRA_type = string(FRA_type);
SAMData.FRA_area = single(FRA_area);

end %end function
