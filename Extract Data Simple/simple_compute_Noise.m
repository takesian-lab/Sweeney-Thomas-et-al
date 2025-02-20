function [NoiseData, fig1] = simple_compute_Noise(IsResponsive, PeakData, StimInfo, plot_figures)
% Compute intensity tuning for noise data with more than one dB
%
% Argument(s): 
%   IsResponsive - 0s and 1s indicating whether the cell is significantly responsive to that stim combination
%   PeakData - from simple_check_if_resopnsive
%   StimInfo - from simple_prepare_variables
%   plot_figures - 0 or 1 to plot figures
%
% Returns:
%   NoiseData table:
%     - BestInt = Intensity (dB) with highest amplitude response
%     - Raw_BestInt = computed ignoring isRF (best values not restricted to resposive & reliable)
%     - ISI = Intensity Selectivity Index (adapted from Mesik et al ISI)
%     - ISI = Intensity Selectivity Index (best values not restricted to resposive & reliable)
%     - FRA_area = % of RF deemed to be responsive

% Version History:
%  - V1 = current version made for simple_extract_data
%
% TODO: 
% Search 'TODO'
%% Skip if data is  not in correct format
% We are expecting data that only varies in the intensity dimension

if ~strcmp(StimInfo.Parameters{1}, 'Intensity') || (numel(StimInfo.V2) > 1)
    return;
end

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
NoiseData = table;
NoiseData.BestInt = single(nan);
NoiseData.Raw_BestInt = single(nan);
NoiseData.ISI = single(nan);
NoiseData.Raw_ISI = single(nan);
NoiseData.FRA_area = single(nan);

%Store response in NoiseData (for plotting)
NoiseData.Response = single(temp_RF_orig');

%Check RF for responsive stim conditions and make sure there is >1 stim condition:
if ~any(IsResponsive) || isequal(size(IsResponsive),[1,1])
    return
end

%% Find best intensity 

%Reorder ints from quietest to loudest
ints = StimInfo.V1;
[ints, stimorder] = sort(ints,'ascend');
RF = RF(stimorder);
IsRF = IsRF(stimorder);

%Raw values
[M,I] = max(RF);
NoiseData.Raw_BestInt = ints(I);

%Thresholded values
thresholded_RF = RF;
thresholded_RF(~IsRF) = nan;
[M,I] = max(thresholded_RF);
NoiseData.BestInt = ints(I);

%% compute FRA area using IsRF

%Add NaNs to IsRF where there were no trial to analyze (e.g. in Running vs. NotRunning)
IsRF_withNans = IsRF;
IsRF_withNans(isnan(RF)) = NaN;

%Compute FRA on non-NaN values only
IsRF_removeNans = IsRF_withNans(:);
IsRF_removeNans(isnan(IsRF_removeNans)) = [];

NoiseData.FRA_area = sum(IsRF_removeNans)/numel(IsRF_removeNans);

%% Compute intensity selectivity  (ISI)
%ISI = Intensity Selectivity Index (1 = responds to lowest dB, 0 = responds to highest dB)

[NoiseData.ISI, ~] = tak_compute_selectivity_index(RF, 1, I, 'PlotFigure', 0);

%Raw ISI
[~, ind] = max(RF);
[NoiseData.Raw_ISI, ~] = tak_compute_selectivity_index(RF, 1, ind, 'PlotFigure', 0);

end %end function
