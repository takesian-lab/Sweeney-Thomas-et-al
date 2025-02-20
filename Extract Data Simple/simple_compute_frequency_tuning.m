function [TuningData, fig1] = simple_compute_frequency_tuning(IsResponsive, PeakData, StimInfo, plot_figures)
% Compute tuning properties of a single cell and plot figures
%
% Argument(s): 
%   IsResponsive - 0s and 1s indicating whether the cell is significantly responsive to that stim combination
%   PeakData - from simple_check_if_resopnsive
%   StimInfo - from simple_prepare_variables
%   plot_figures
%
% Returns:
%   TuningData:
%     - BF          = best frequency, freq with the strongest response (kHz)
%     - BF_I        = intensity with the strongest response (dB), equivalent to intestity at BF
%     - BWBF_I      = width of the gaussian tuning curve at BF at half the maximum (octaves) computed using all raw values
%     - BWBF_BinaryMask = BWBF_I computed with RF binary mask applied (~RF set to 0)
%     - CF          = frequency with the strongest response at threshold (kHz)
%     - CF_I        = threshold, intensity corresponding to CF (dB)
%     - BW20        = width of the gaussian tuning curve 20dB above threshold at half the maximum (octaves)
%     - BWInt       = (measured at BF) width of the gaussian intensity tuning curve at half the maximum (dB)
%     - BWInt_BinaryMask = BWInt computed with RF binary mask applied (~RF set to 0)
%     - ISI         = Intensity Selectivity Index (from Mesik et al)
%     - dprime      = sensitivity index (from Romero et al) 
%     - Raw_dprime  =  " " (best value not restricted to responsive and reliable)
%     - Binary_dprime = " " computed on zeros and ones
%     - Raw_BF      = compute BF ignoring isRF (best values not restricted to responsive & reliable)
%     - Raw_BF_I    = compute BF_I  "      "      
%     - Raw_BWInt   = compute BWInt "      "     
%     - Raw_ISI     = compute ISI    "      " 
%     - FRA_type    = simple or complex (multipeaked)
%     - FRA_area    = % of RF deemed to be reliable and responsive
% 
% Version history:
% - V1 = current version based off of compute_frequency_tuning for extract_data pipeline
%
% TODO: 
% Search 'TODO'

%% Setup

%Get RF mapping from StimInfo
RF_Map = StimInfo.RF_Map;
nanmat = nan(max(RF_Map));
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

%Reconstruct RF maps
[IsRF, RF] = deal(nanmat);
IsRF(RF_ind) = IsResponsive;

%Which PeakData to use for Response value?
%--> Use AUC, but make suppressed AUCs negative for proper inhibitory/excitatory RF detection
%--> If ALL significant responses are suppressed, flip the signs so that we can measure suppressed RFs

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
temp_RF_orig = temp_RF; %Save for storing in TuningData without sign flip
RF_type = determine_RF_response_type(PeakData.ResponseType, IsResponsive);
if strcmp(RF_type, 'inhibitory')
    temp_RF = -temp_RF;
    RF = -RF;
end

%Variables will return as NaN if fits cannot be made
fig1 = [];
TuningData                  = table;
TuningData.BF               = single(nan);
TuningData.BF_I             = single(nan);
TuningData.BWBF_I           = single(nan);
TuningData.BWBF_BinaryMask  = single(nan);
TuningData.CF               = single(nan);
TuningData.CF_I             = single(nan);
TuningData.BW20             = single(nan);
TuningData.BWInt            = single(nan);
TuningData.BWInt_BinaryMask = single(nan);
TuningData.ISI              = single(nan);
TuningData.dPrime           = single(nan);
TuningData.Raw_dPrime       = single(nan);
TuningData.Binary_dPrime    = single(nan);
TuningData.Raw_BF           = single(nan);
TuningData.Raw_BF_I         = single(nan);
TuningData.Raw_BWBF_I       = single(nan);
TuningData.Raw_BWInt        = single(nan);
TuningData.Raw_ISI          = single(nan);
TuningData.FRA_type         = "none"; %Double quotes = string
TuningData.FRA_area         = single(nan);

%Store response in TuningData (for plotting)
TuningData.Response = single(temp_RF_orig');

%convert freqs to octaves
freqs = StimInfo.V1;
ints = StimInfo.V2;
octaves = log2(freqs);

%Check RF for responsive stim conditions and make sure there is >1 stim condition:
if ~any(any(IsRF)) || isequal(size(IsRF),[1,1])
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

%% Find threshold and BF

%find CF_I (threshold)
if ismember(RF_type, 'inhibitory')
    [row,col] = find(inhR_RF == 1);
else
    [row,col] = find(excR_RF == 1); % only excitatory and responsive
end    

[CF_row] = max(row);
CF_I = ints(CF_row);

%find BF
thresholded_RF = RF;
thresholded_RF(~IsRF) = nan;
[M,I] = max(thresholded_RF(:));
[BF_row, BF_col] = ind2sub(size(RF),I);
BF_I = ints(BF_row);
BF = freqs(BF_col);

%find raw BF (without using RF mask)
[M2,I2] = max(RF(:));
[raw_BF_row, raw_BF_col] = ind2sub(size(RF),I2);
Raw_BF_I = ints(raw_BF_row);
Raw_BF = freqs(raw_BF_col);

%% compute FRA area using IsRF

%Add NaNs to IsRF where there were no trial to analyze (e.g. in Running vs. NotRunning)
IsRF_withNans = IsRF;
IsRF_withNans(isnan(RF)) = NaN;

%Compute FRA on non-NaN values only
IsRF_removeNans = IsRF_withNans(:);
IsRF_removeNans(isnan(IsRF_removeNans)) = [];

FRA_area = sum(IsRF_removeNans)/numel(IsRF_removeNans);

%% d prime sensitivity index (from Romero, Polley et al Cerebral Cortex 2019)

dprime = tak_compute_RF_dprime(RF, BF_row, BF_col);
Raw_dprime = tak_compute_RF_dprime(RF, raw_BF_row, raw_BF_col);
Binary_dprime = tak_compute_RF_dprime(~isnan(thresholded_RF), BF_row, BF_col);

%% Plot RF and FRA

if plot_figures
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    
    subplot(2,3,1); %Plot RF for reference
    RFmask = ones(size(nanmat))*0.1;
    RFmask(IsRF == 1) = 1;
    imagesc(RF_orig, 'AlphaData', ones(size(nanmat))); %Plot without mask first
    colorbar
    xlabel('Frequency (kHz)')
    ylabel('Intensity (dB)')
    set(gca,'XTick',1:length(freqs))
    set(gca,'YTick',1:length(ints))
    set(gca,'XTickLabel',num2str(freqs))
    set(gca,'YTickLabel',num2str(ints))
    title('Receptive field')
    subtitle(['dPrime =  ' num2str(round(dprime,1))]);
    text(BF_col - 0.25, BF_row, 'BF', 'Color', 'w') %-0.25 centers the text
    
    colormap(gca, bluewhitered(256)); %TODO: use bluewhitered_TAKlab
    c = colorbar;
    c.Title.String = 'Z';

    subplot(2,3,4) %plot FRA contour
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
    xlabel('Frequency (kHz)')
    ylabel('Intensity (dB)')
    set(gca,'XTick',1:length(freqs))
    set(gca,'YTick',1:length(ints))
    set(gca,'XTickLabel',num2str(freqs))
    set(gca,'YTickLabel',num2str(ints))
    title('Frequency Receptive Area');
    subtitle(strcat({'Type =  '}, string(FRA_type), {' '}, string(RF_type)));
    
    colormap(gca, bluewhitered(256)); %TODO: use bluewhitered_TAKlab
    c = colorbar;
    c.Title.String = 'Z';
end

%% Compute intensity tuning (ISI)

[~, ISI] = tak_compute_selectivity_index(RF, BF_row, BF_col, 'PlotFigure', 0);
[~, Raw_ISI] = tak_compute_selectivity_index(RF, raw_BF_row, raw_BF_col, 'PlotFigure', 0);

%% Fit #1 - Intensity Tuning
% determine intensity tuning using peak of gaussian fit across intensities at BF

best_ind = zeros(size(ints));
best_ind(BF_row) = 1;
FitTable = tak_fit_gaussian(ints, RF(:,BF_col), best_ind, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

%Currently tak_fit_gaussian constrains the peak to the max value, so FitTable BestInt and ints(BF_row) should always be the same
if ~isempty(FitTable)
    BestInt = FitTable.Best_X;
    Int_halfwidth = FitTable.Width;
else
    BestInt = ints(BF_row);
    Int_halfwidth = nan;
end

if plot_figures
    subplot(2,3,2); hold on %plot gaussian fits    
    if isempty(FitTable)
        plot(ints, RF(:,BF_col), 'Linewidth', 2);
        xlim([ints(end) ints(1)])
        vline(BestInt)
    else
        tak_plot_gaussian(FitTable)
    end
    title(['Best Intensity at ' num2str(BF) 'kHz = ' num2str(BestInt) 'dB']);
    subtitle(['ISI = ' num2str(round(ISI,1))])
    xlabel('Intensity (dB)');
    ylabel('Response');    
end

%Fit width at best Int with binary mask
best_ind = zeros(size(ints));
best_ind(BF_row) = 1;
thresholdedRF0 = thresholded_RF;
thresholdedRF0(isnan(thresholded_RF)) = 0;
FitTable = tak_fit_gaussian(ints, thresholdedRF0(:,BF_col), best_ind, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Int_halfwidth_binary_mask = FitTable.Width;
else
    Int_halfwidth_binary_mask = nan;
end

%Fit width at best Int with raw data
best_ind = zeros(size(ints));
best_ind(raw_BF_row) = 1;
FitTable = tak_fit_gaussian(ints, RF(:,raw_BF_col), best_ind, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Raw_BWInt = FitTable.Width;
else
    Raw_BWInt = nan;
end

%% Fit #2 - CF
%determine CF using peak of gaussian fit at threshold intensity

%Threshold CF row by responsive/reliable to find best_ind
CFfitrow = RF(CF_row,:)';
IsRF_CFfitrow = IsRF(CF_row,:)';
CFfitrow(~IsRF_CFfitrow) = nan;
[~, ind] = max(CFfitrow);
best_ind = zeros(size(octaves));
best_ind(ind) = 1;

FitTable = tak_fit_gaussian(freqs, RF(CF_row,:)', best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

%Currently tak_fit_gaussian constrains the peak to the max value, so FitTable CF and freqs(ind) should always be the same
if ~isempty(FitTable)
    CF_kHz = FitTable.Backtransformed_Best_X;
else
    CF_kHz = freqs(ind);
end

if plot_figures
    subplot(2,3,3); hold on %plot gaussian fits
    if isempty(FitTable)
        plot(octaves, RF(CF_row,:), 'Linewidth', 2);
        set(gca, 'XTick', octaves)
        set(gca, 'XTickLabel', num2str(freqs))
        xlim([octaves(1) octaves(end)])
    else
        tak_plot_gaussian(FitTable)
    end
    title(['CF at Threshold ' num2str(ints(CF_row)) 'dB = ' num2str(round(CF_kHz,1)) 'kHz']);
    xlabel('Frequency (kHz)');
    ylabel('Response');
end

%% Fit #3 - BW20 
%find bandwidth 20dB above threshold
if CF_row - 2 > 0 %Check if BW20 exists
    
    %Try to threshold BW20 row by responsive/reliable
    BW20row = RF(CF_row-2,:)';
    IsRF_BW20row = IsRF(CF_row-2,:)';
    BW20row(~IsRF_BW20row) = nan;
    
    %It is possible none are responive/reliable
    if all(isnan(BW20row))
        best_ind = [];
    else
        %But if any are, set that as the best peak ind
        [~, ind] = max(BW20row);
        best_ind = zeros(size(octaves));
        best_ind(ind) = 1;
    end
    
    FitTable = tak_fit_gaussian(freqs, RF(CF_row-2,:)', best_ind, 'LogTransform', 1, 'PlotFigure', 0);
    FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

    if ~isempty(FitTable)
        BW_20_halfwidth = FitTable.Width; %Keep width in octaves, do not use backtransformed width
    else
        BW_20_halfwidth = nan;
    end

    if plot_figures
        subplot(2,3,5); hold on %plot gaussian fits
        if isempty(FitTable)
            plot(octaves, RF(CF_row-2,:), 'Linewidth', 2);
            set(gca, 'XTick', octaves)
            set(gca, 'XTickLabel', num2str(freqs))
        else
            tak_plot_gaussian(FitTable)
        end
        title(['BW20 at ' num2str(ints(CF_row-2)) 'dB = ' num2str(round(BW_20_halfwidth,1)) ' octaves']);
        xlabel('Octaves');
        ylabel('Response');
    end
else
    BW_20_halfwidth = nan;
end

%% Fit #4 - bandwidth at BF
%find bandwidth at BF

best_ind = zeros(size(octaves));
best_ind(BF_col) = 1;
FitTable = tak_fit_gaussian(freqs, RF(BF_row,:)', best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    BW_BF_halfwidth = FitTable.Width;
else
    BW_BF_halfwidth = nan;
end

if plot_figures
    subplot(2,3,6); hold on %plot gaussian fits    
    if isempty(FitTable)
        plot(octaves, RF(BF_row,:), 'Linewidth', 2);
        set(gca, 'XTick', octaves)
        set(gca, 'XTickLabel', num2str(freqs))
        xlim([octaves(1) octaves(end)])
    else
        tak_plot_gaussian(FitTable)
    end
    title(['Bandwidth at BF ' num2str(round(BF,1)) 'kHz = ' num2str(round(BW_BF_halfwidth,1)) ' octaves']);
    xlabel('Octaves');
    ylabel('Response');
end

%Fit bandwidth at BF with binary mask
best_ind = zeros(size(octaves));
best_ind(BF_col) = 1;
thresholdedRF0 = thresholded_RF;
thresholdedRF0(isnan(thresholded_RF)) = 0;
FitTable = tak_fit_gaussian(freqs, thresholdedRF0(BF_row,:)', best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    BW_BF_binary_mask = FitTable.Width;
else
    BW_BF_binary_mask = nan;
end

%Fit bandwidth at BF with raw data
best_ind = zeros(size(octaves));
best_ind(raw_BF_col) = 1;
FitTable = tak_fit_gaussian(freqs, RF(raw_BF_row,:)', best_ind, 'LogTransform', 1, 'PlotFigure', 0);
FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit

if ~isempty(FitTable)
    Raw_BWBF_I = FitTable.Width;
else
    Raw_BWBF_I = nan;
end
%% Store Data

TuningData.BF                = single(BF);
TuningData.BF_I              = single(BF_I);
TuningData.BWBF_I            = single(BW_BF_halfwidth);
TuningData.BWBF_BinaryMask   = single(BW_BF_binary_mask);
TuningData.CF                = single(CF_kHz);
TuningData.CF_I              = single(CF_I);
TuningData.BW20              = single(BW_20_halfwidth);
TuningData.BWInt             = single(Int_halfwidth);
TuningData.BWInt_BinaryMask  = single(Int_halfwidth_binary_mask);
TuningData.ISI               = single(ISI);
TuningData.dPrime            = single(dprime);
TuningData.Raw_dPrime        = single(Raw_dprime);
TuningData.Binary_dPrime     = single(Binary_dprime);
TuningData.Raw_BF            = single(Raw_BF);
TuningData.Raw_BF_I          = single(Raw_BF_I);
TuningData.Raw_BWBF_I        = single(Raw_BWBF_I);
TuningData.Raw_BWInt         = single(Raw_BWInt);
TuningData.Raw_ISI           = single(Raw_ISI);
TuningData.FRA_type          = string(FRA_type);
TuningData.FRA_area          = single(FRA_area);

end %end function