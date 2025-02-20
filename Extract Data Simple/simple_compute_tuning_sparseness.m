function [SparsenessData, fig1] = simple_compute_tuning_sparseness(PeakData, IsResponsive, StimInfo, plot_figures, varargin)
% Calculate lifetime sparseness (for 2D matrices only)
%
% Isaacson papers: Lin et al 2019 PNAS Arousal Regulates Frequency
% Tuning... and Kato, Asinof, Isaacson 2017 Neuron Network-Level Control...
%
% Lifetime sparseness (Rolls and Tovee, 1995; Willmore and Tolhurst, 2001), which was calculated as:
% (1 - {[sum of Rj/N]^2 / [sum of Rj^2/N]}) / (1 - 1/N)
% where Rj was the response peak amplitude of the cell to tone j, and N was the total number of tones.
% 1 – Sp provides a measure of how much the response probability of a neuron was distributed equally among all tones
% (non-selective: 1 - Sp = 1) versus attributable entirely to one tone (highly selective: 1- Sp = 0).
% Since calculation of 1 - Sp does not rely on the TRF shape, it can also be used for quantifying the selectivity of
% inhibitory responses, which often do not have clear V-shape.
%
% Argument(s): 
%   PeakData - From simple_extract_data
%   IsResponsive (X x Y matrix) - matrix of 0s and 1s corresponding to if cell is responsive to that stimuli
%   StimInfo - From simple_extract_data
%   plot_figures - 0 or 1 to plot figure
%
% Returns:
% SparsenessData
% 
% Version history:
% - V1
% - V2 = current version
%
% TODO: 
% Search 'TODO'

%% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% Use custom RF values instead of PeakData.Peak_AUC/Trough_AUC
% **NEEDS TO BE A 1D ARRAY -> Make sure to plot figure to confirm ordering is correct**
addOptional(p, 'RF', []); 

% Set the RF floor to 0 every time (1 *RECOMMENDED*) or only set negative values to 0 (0)
addOptional(p, 'MakeFloor0', 1);

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% Establish variable to use for sparseness
%Isaacson uses peak. Lets try our Peak_3fr variable

%Get RF mapping from StimInfo
RF_Map = StimInfo.RF_Map;
nanmat = nan(max(RF_Map));
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

%Fill RF
RF = nanmat;

if isempty(ops.RF)
    %Generate RF from peaks and troughs
    temp_RF = nan(size(RF_Map,1),1);
    for i = 1:length(temp_RF)
        response_type = PeakData.ResponseType(i);
        if strcmp(response_type, 'suppressed')
            temp_RF(i) = -PeakData.Trough_AUC(i);
        else
            temp_RF(i) = PeakData.Peak_AUC(i);
        end
    end
    RF(RF_ind) = temp_RF;
    
    %If all significant responses are suppressed, flip the sign of RF
    RF_type = determine_RF_response_type(PeakData.ResponseType, IsResponsive);
    if strcmp(RF_type, 'inhibitory')
        RF = -RF;
    end
else
    RF(RF_ind) = ops.RF; %Use user-defined values for RF
end

%% Set inf to nan

RF(isinf(RF)) = nan;

%% Remove negatives (If you don't do this, sp can be > 1)

if ops.MakeFloor0
    %Make floor 0 (any negative values will mess up computation)
    if any(any(RF < 0))
        RF = RF + abs(min(min(RF)));
    elseif all(all(RF > 0))
        RF = RF - min(min(RF));
    end
else
    %Only set negative values to zero
    RF(RF < 0) = 0;
end

%% Calculate sparseness of full matrix not including NaNs
R = RF(~isnan(RF));
N = numel(R); %N = number of stimuli
E1 = sum((R/N))^2;
E2 = sum(R.^2/N);
if E2 == 0 %Can't divide by zero
    sp = 0;
else
    sp = (1 - E1/E2)/(1 - 1/N);
end

%% Calculate V1 sparseness by V2
% e.g. freq. sparseness at each dB intensity

sp_by_V2 = nan(size(RF,1),1);
for i = 1:size(RF,1)
    R = RF(i,:);
    R = R(~isnan(R));
    if isempty(R)
        continue;
    end
    N = length(R); %N = number of tones
    E1 = sum((R/N))^2;
    E2 = sum(R.^2/N);
    if E2 == 0 %Can't divide by zero
        sp_by_V2(i) = 0;
    else
        sp_by_V2(i) = (1 - E1/E2)/(1 - 1/N);
    end
end

%% Calculate V2 sparseness by V1
% e.g. int. sparseness at each freq

sp_by_V1 = nan(size(RF,2),1);
for i = 1:size(RF,2)
    R = RF(:,i);
    R = R(~isnan(R));
    if isempty(R)
        continue;
    end
    N = length(R); %N = number of tones
    E1 = sum((R/N))^2;
    E2 = sum(R.^2/N);
    if E2 == 0 %Can't divide by zero
        sp_by_V1(i) = 0;
    else
        sp_by_V1(i) = (1 - E1/E2)/(1 - 1/N);
    end
end

%% Catch errors

if any(sp_by_V1 > 1) || any(sp_by_V2 > 1)
    error('Found sparseness values greater than 1')
end

%% Figure

if plot_figures
    
    V1 = StimInfo.V1;
    V2 = StimInfo.V2;
    V1_label = StimInfo.Parameters{1};
    V2_label = StimInfo.Parameters{2};
    
    fig1 = figure; hold on
    subplot(2,2,3)
    imagesc(RF); hold on
    ylabel(V2_label)
    xlabel(V1_label)
    set(gca,'YTick',1:length(V2))
    set(gca,'YTickLabel', V2)
    set(gca,'XTick',1:length(V1))
    set(gca,'XTickLabel', V1)
    title(['Overall sp. = ' num2str(round(sp,1))])
    colormap(gca, bluewhitered(256))
    c = colorbar;
    c.Title.String = 'Z';
    
    subplot(2,2,4)
    plot(sp_by_V2); hold on
    scatter(1:length(sp_by_V2),sp_by_V2)
    camroll(-90)
    xlabel(V2_label)
    ylabel([V1_label ' sparseness by ' V2_label])
    ylim([0 1])
    xlim([0.5 length(V2)+0.5])
    set(gca,'XTick',1:length(V2))
    set(gca,'XTickLabel', V2)
    
    subplot(2,2,1)
    plot(sp_by_V1); hold on
    scatter(1:length(sp_by_V1),sp_by_V1)
    xlabel(V1_label)
    ylabel([V2_label ' sparseness by ' V1_label])
    ylim([0 1])
    xlim([0.5 length(V1)+0.5])
    set(gca,'XTick',1:length(V1))
    set(gca,'XTickLabel', V1)
else
    fig1 = [];
end

%% Make table

SparsenessData = table;
SparsenessData.Sparseness(1) = single(sp);
SparsenessData.V1_Sparseness{1} = single(sp_by_V2);
SparsenessData.V2_Sparseness{1} = single(sp_by_V1);
SparsenessData.V1_Mean_Sparseness(1) = single(mean(sp_by_V2,'omitnan'));
SparsenessData.V2_Mean_Sparseness(1) = single(mean(sp_by_V1,'omitnan'));

%% USE THIS CODE TO TEST FUNCTION (run in separate script)
% 
% plotFigure = 1;
% units = 'A.U.';
% freqs = [4, 5.7, 8, 11.3, 16, 22.6, 32, 45.2];
% ints = [80 70 60 50 40 30 20 10];
% cellNumber = 0;
% 
% RF = zeros(8); %Divide by zero error: nans
% RF = ones(8); %Sparseness = 0
% RF = -ones(8); %Same as above
% RF = eye(8); %Sparseness = 1 at all freqs and ints
% RF = flipud(eye(8)); %Same as above
% RF = -eye(8); %Same as above
% RF = rand(8);
% RF = -rand(8);
% RF = zeros(8); RF(randi(64,10,1)) = 1;
% 
% [sp, sp_by_int, sp_by_freq, h] = compute_tuning_sparseness_v2(RF, plotFigure, units, freqs, ints, cellNumber);

end