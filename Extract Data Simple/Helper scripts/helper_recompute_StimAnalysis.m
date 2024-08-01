%% Helper function: Recompute StimAnalysis for ExtractedData
% This script recomputes only the sparseness and StimAnalysis measures for
% a previously saved ExtractedData file, without changing the identity of
% sound-responsive/reliable cells.

% Update September 19th: Add reliability to every stim condition (without
% controls) to make Reliability Tuning maps

%% Load ExtractedData first
clearvars -except ExtractedData

%Options
save_path = pwd;
save_new_ExtractedData = 1;
save_summary_figures = 0;

%% Recompute StimAnalysis

Ops = ExtractedData.Ops;
StimInfo = ExtractedData.StimInfo;
Loco = ExtractedData.Ops.loco_filter;
CellList = ExtractedData.CellList;
CellDataByStim = ExtractedData.CellDataByStim;

oldBlock = '';

for c = 1:size(CellList,1)
    
    blockName = char(CellList.Block(c));
    if ~strcmp(oldBlock, blockName)
        oldBlock = blockName;
        disp(oldBlock)
    end
    cellNumber = CellList.CellNumber(c);
    cc = CellList.CellOrder(c);
    disp(cellNumber)
    
    %% Add reliability to every stim condition
    if Ops.use_reliability && Ops.reliability_type == "Pearsons"
        for L = 1:length(Loco)
            for v = 1:size(StimInfo.Combined,1)
                ReliabilityOps = Ops; %Duplicate Ops for use in simple_reliability function
                ReliabilityOps.ComputeReliabilityControls = 0; %DO NOT RECOMPUTE CONTROLS
                
                %Compute reliability by comparing trials to either blank trials or shifted traces
                trials = CellDataByStim.(Loco{L})(c).StimTraces{v};
                blanks = CellDataByStim.(Loco{L})(c).BlankTraces; %If empty, script will automatically use shifted traces
                ReliabilityData = simple_reliability(trials, blanks, StimInfo, ReliabilityOps);
                
                %Add R and Pearsons to ReliabilityData
                CellDataByStim.(Loco{L})(c).ReliabilityData(v).R = ReliabilityData.R;
                CellDataByStim.(Loco{L})(c).ReliabilityData(v).Pearsons = ReliabilityData.Pearsons;
            end
        end
    end
    
    %% Compute sparseness
    SparsenessData = struct;
    BinarySparsenessData = struct;
    
    for L = 1:length(Loco)
        switch Ops.stim_protocol
            case {2 3 5 6} %These are the stim that have been tested so far, other stim could also work here as long as they have > 1 stim parameter
                %Compute tuning sparseness
                [SparsenessData.(Loco{L}), ~] = simple_compute_tuning_sparseness(CellDataByStim.(Loco{L})(c).PeakData, CellDataByStim.(Loco{L})(c).RF, StimInfo, Ops.plot_figures);

                %Compute again using just zeros and ones for RF
                TempRF = single(CellDataByStim.(Loco{L})(c).RF);
                TempRF(CellDataByStim.(Loco{L})(c).nTrials == 0) = NaN; %Correct RF to include NaNs where there were no trials
                [BinarySparsenessData.(Loco{L}), ~] = simple_compute_tuning_sparseness(CellDataByStim.(Loco{L})(c).PeakData, CellDataByStim.(Loco{L})(c).RF, StimInfo, Ops.plot_figures, 'RF', TempRF);
                BinarySparsenessData.(Loco{L}).Properties.VariableNames =  ["Binary_Sparseness","Binary_V1_Sparseness","Binary_V2_Sparseness","Binary_V1_Mean_Sparseness","Binary_V2_Mean_Sparseness"];
        end
    end
    
    %% Perform stim-specific receptive field analyses
    %Compute tuning, combine with sparseness data, and store results in StimAnalysis
    
    for L = 1:length(Loco)
        
        switch Ops.stim_protocol

            case 2 %RF
                [TuningData, ~] = simple_compute_frequency_tuning(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), TuningData];

            case 3 %FM
                R = nan(size(CellDataByStim.(Loco{L})(c).RF));
                if Ops.use_reliability && Ops.reliability_type == "Pearsons"
                    R = [CellDataByStim.(Loco{L})(c).ReliabilityData.R]';
                end
                [FMData, ~] = simple_compute_FM(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, R, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), FMData];

            case 5 %SAM
                [SAMData, ~] = simple_compute_SAM(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), SAMData];

            case 6 %SAMfreq
                [SAMFreqData, ~] = simple_compute_SAMfreq(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), SAMFreqData];
                
            case 10
                %Compute intensity tuning for noise stim with more than one dB
               [StimAnalysis.(Loco{L})(c,:), ~] = simple_compute_Noise(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                
            case 12 %SPONTANEOUS
                %Ca_detection pulls out spontaneous calcium transients and calculates their amplitude, frequency, and duration
                [~, StimAnalysis.(Loco{L})(c,:), ~] = simple_Ca_detection_spontaneous(TrialData(b).Timestamp, TrialData(b).F7, TrialData(b).Full_Loco, StimInfo.fs, cc, Loco{L}, Ops.plot_figures);

            case 13 %MARYSE BEHAVIOR STIM
                [StimAnalysis.(Loco{L})(c,:), ~] = simple_compute_tuning_maryse_behavior(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
        end
    end
    
end

ExtractedData.CellDataByStim = CellDataByStim;
ExtractedData.StimAnalysis = StimAnalysis;
ExtractedData.Summary = make_simple_extracted_data_summary(ExtractedData); %Put main results into stats table
   
%% Save and plot figures

if save_new_ExtractedData
    cd(save_path)
    disp('Saving data...')
    
    filename = ExtractedData.Filename;
    save([filename(1:end-4) '_newStimAnalysis.mat'], 'ExtractedData', '-v7.3');
    for L = 1:length(Loco)
        writetable(ExtractedData.Summary.(Loco{L}),[filename(1:end-4) '_newStimAnalysis_' Loco{L} '_Summary.csv'],'Delimiter',',');
    end
end

%% Plot summary figures

if save_summary_figures
    cd(save_path)
    disp('Saving figures...')
    
    %These functions use simple_subset_ExtractedData to filter plotted cells by the folowing parameters.
    %You can adjust them here:
    PlotOps = struct;
    PlotOps.Loco                = 'All'; %All, Running, or NotRunning
    PlotOps.ResponseType        = 'activated'; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
    PlotOps.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
    PlotOps.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
    PlotOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1)
    PlotOps.Group               = ''; %'' -> no filtering, or add group name here to only keep that group
    PlotOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
    PlotOps.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    PlotOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    PlotOps.plotCorrelations	= 0; %0 = don't plot dependent variable correlations, 1 = plot
    PlotOps.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
    PlotOps.FibPhotChannel      = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
    PlotOps.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress

    plot_simple_ExtractedData(ExtractedData, PlotOps);               %Rasters and pie charts for activated vs. prolonged. vs. suppressed cells
    plot_simple_ExtractedData_Reliability(ExtractedData, PlotOps);   %Scatter plots with reliability measures (function will skip if no reliability computed)
    plot_simple_ExtractedData_PeakData(ExtractedData, PlotOps);      %Peak data histograms and rasters
    plot_simple_ExtractedData_Controls(ExtractedData, PlotOps);      %Plot results vs. sex, age, auditory field, GCaMP type, etc

    %Stim-specific plots (function will skip if not the right data type)
    plot_simple_ExtractedData_RF(ExtractedData, PlotOps);
    plot_simple_ExtractedData_FM(ExtractedData, PlotOps);
    plot_simple_ExtractedData_SAM(ExtractedData, PlotOps);
    plot_simple_ExtractedData_SAMfreq(ExtractedData, PlotOps);
    plot_simple_ExtractedData_MaryseBehavior(ExtractedData, PlotOps);
    plot_simple_ExtractedData_AllBehavior(ExtractedData, PlotOps);

    %Stim-specific plots split by activated, prolonged, suppressed 
    plot_simple_ExtractedData_RF_APS(ExtractedData, PlotOps);
    plot_simple_ExtractedData_FM_APS(ExtractedData, PlotOps);
    plot_simple_ExtractedData_SAM_APS(ExtractedData, PlotOps);
    plot_simple_ExtractedData_SAMfreq_APS(ExtractedData, PlotOps);

    %Save all open figures as FIG, PNG, and EPS
    [~, ~, ~] = mkdir(save_path,'SummaryFigures_newStimAnalysis');
    cd([save_path '/SummaryFigures_newStimAnalysis'])
    save_plots(pwd, ['SummaryFig_' ExtractedData.StimInfo.StimType '_' ExtractedData.Date])
end

%%

disp('All done!')