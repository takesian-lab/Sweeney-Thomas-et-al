%helper_script_to_test_ISI

%QUESTIONS:
%How many cells have a BF value lower than the greatest value?
%What do the results look like if we set all - values to 0
%What do the results look like if we always set the min value to zero
%What do the results look like if we make the norm range 0 to 1?

%What about 0 floor for sparseness??

%Load RF ExtractedData file first

%% Subset data to only include responsive cells
Ops = struct;
Ops.Loco                = 'All'; %All, Running, or NotRunning
Ops.ResponseType        = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
Ops.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
Ops.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
Ops.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
Ops.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
Ops.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
Ops.SuppressOutput      = 1; %0 to print Ops to command line, 1 to suppress

SubsetData = simple_subset_ExtractedData(ExtractedData, Ops);

StimInfo = SubsetData.StimInfo;
StimAnalysis = SubsetData.StimAnalysis;
PeakData = SubsetData.CellDataByStim;

%% Recompute ISI to compare different ways of normalizing

%Get RF mapping from StimInfo
RF_Map = StimInfo.RF_Map;
nanmat = nan(max(RF_Map));
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

%Normalization options
ISI_nan = nan(size(StimAnalysis,1),1); %Do not compute
ISI_fixneg  = nan(size(StimAnalysis,1),1); %Only set negative values to 0
ISI_floor0  = nan(size(StimAnalysis,1),1); %Always set floor to 0

for c = 1:size(StimAnalysis,1)
    
    %Reconstruct RF maps
    RF = nanmat;
    RF(RF_ind) = StimAnalysis.Response(c,:);
    
    %Figure out BF ind
    BF_col = find(StimInfo.V1 == StimAnalysis.BF(c));
    BF_row = find(StimInfo.V2 == StimAnalysis.BF_I(c));

    [~, ISI_nan(c)] = tak_compute_selectivity_index(RF, BF_row, BF_col, 'MakeFloor0', 0, 'NormalizeToBF', 0, 'PlotFigure', 0);
    [~, ISI_fixneg(c)] = tak_compute_selectivity_index(RF, BF_row, BF_col, 'MakeFloor0', 0, 'NormalizeToBF', 1, 'PlotFigure', 0);
    [~, ISI_floor0(c)] = tak_compute_selectivity_index(RF, BF_row, BF_col, 'MakeFloor0', 1, 'NormalizeToBF', 1, 'PlotFigure', 0);
end

%% Results

%Number of cells affected by this
percentcells = sum(isnan(ISI_nan))/length(ISI_nan); %34 out of 632
nanind = isnan(ISI_nan);

%Distributions
figure;
subplot(2,2,1); hold on
histogram(ISI_fixneg)
%histogram(ISI_fixneg(nanind == 1))
title('Old method - make negatives 0')

subplot(2,2,3); hold on
histogram(ISI_floor0)
%histogram(ISI_floor0(nanind == 1))
title('New method - make floor 0')

subplot(2,2,[2 4]); hold on
scatter(ISI_fixneg, ISI_floor0, 'jitter', 'off')
scatter(ISI_fixneg(nanind == 1), ISI_floor0(nanind == 1), 'r', 'jitter', 'off')
xlabel('Old method')
ylabel('New method')

sgtitle('ISI')

%Takeaways: with NormalizeToBF, all cells that have stronger responses at
%highest intensity than BF will have ISI 0 because we set highest response == BF

%With old method of correcting zeros, we had a lot of ISI = 1 because we
%were probably setting lots of the RF to zero. Whereas with the floor
%method, we only make 1 value zero and the rest get shifted accordingly, so
%we end up with less 0/BF

%% Recompute sparseness to compare different ways of normalizing

%Normalization options
Sp_fixneg  = nan(size(PeakData,2),1); %Only set negative values to 0
Sp_floor0  = nan(size(PeakData,2),1); %Always set floor to 0

for c = 1:size(PeakData,2)
    [SparsenessData, ~] = simple_compute_tuning_sparseness(PeakData(c).PeakData, PeakData(c).RF, StimInfo, 0, 'MakeFloor0', 0);
    Sp_fixneg(c) = SparsenessData.Sparseness;
    
    [SparsenessData, ~] = simple_compute_tuning_sparseness(PeakData(c).PeakData, PeakData(c).RF, StimInfo, 0, 'MakeFloor0', 1);
    Sp_floor0(c) = SparsenessData.Sparseness;
end

%% Results

%Distributions
figure;
subplot(2,2,1); hold on
histogram(Sp_fixneg)
title('Old method - make negatives 0')

subplot(2,2,3); hold on
histogram(Sp_floor0)
title('New method - make floor 0')

subplot(2,2,[2 4]); hold on
scatter(Sp_fixneg, Sp_floor0, 'jitter', 'off')
xlabel('Old method')
ylabel('New method')

sgtitle('Sparseness')