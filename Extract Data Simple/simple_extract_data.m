% simple_extract_data
%
% Combine data across blocks and prepare variables for statistical analysis
% After data has been extracted, use visualize_extracted_data to plot figures for individual cells
% You can also run visualize_responsive_cells on individual blocks to preview sound responses
%
% Version history:
% - V1 = continued from extract_data_v11, archived Sept 2022
% - V2 = current version with GUI
%
% Authors: Takesian Lab 2022 anne_takesian@meei.harvard.edu

%% Launch GUI and set up environment

clear
close all

% Pop open GUI to change extraction parameters
% This script will PAUSE until the GUI window is closed (either via 'Run' or top-right 'Exit')
% Questions? ask 2P code Slack or Latane (latanebullock@gmail.com)
app = simple_extract_data_GUI();
while isvalid(app); pause(0.1); end
clearvars app

%If the user opened the GUI and exited without clicking 'Run...', abort execution. The script stops here
if ~exist('Ops', 'var')
    warning("No 'Ops' variable found. Aborting script. Please click 'Run...' in the GUI to execute the script with the variables in the GUI.")
    return
end

%Print Ops to command window
disp(Ops)

%Transform Ops if necessary
Ops.stim_protocol = tak_stimName2stimID(Ops.stim_protocol); %convert from stimulus string to protocol ID
Loco = Ops.loco_filter;

%For TESTING purposes
c_end = inf; %set an upper limit on the number of cells. inf -> no limit, 5 -> compute 5

%% Load data

%Load info sheet and select files corresponding to stim_protocol
Info = tak_read_info_table(Ops.info_file_path, 'StimProtocol', Ops.stim_protocol, 'StripFields', false);

%Create data structure and load blocks
disp('Loading blocks to include...')
[BlockInfo, Blocks, CellList] = fillSimpleDataTableFromInfo(Info, Ops.blocks_path, 'StimProtocol', Ops.stim_protocol, 'ReturnBlocks', 1, 'ReturnCellList', 1);

%Add match info
if Ops.include_matches
    CellList = add_matches_to_CellList(CellList, Ops.matching_path);
end

%Concatenate blocks:
[Blocks, BlockInfo, CellList, Archive] = simple_Concat_Blocks(Ops, Blocks, BlockInfo, CellList);

%Pull out all single-trial and stim data from blocks
[TrialData, StimInfo] = makeTrialList(BlockInfo, Blocks, Ops.save_full_trace);

%Add loco percent and running speed to BlockInfo
BlockInfo = add_loco_to_BlockInfo(BlockInfo, TrialData);

%Remove unbalanced stim (stim that not all mice are presented with) from trial list
%StimInfo will include a list of removed stim
if Ops.remove_unbalanced_stim
    [TrialData, StimInfo] = removeUnbalancedStim(TrialData, StimInfo);
    
    %Maryse and Carolyn special option to use 0% SAM as noise data
    %[TrialData, StimInfo] = removeUnbalancedStim(TrialData, StimInfo, 'KeepOnlyStim', [0 0]);
end

%If Markpoints/Activation analysis, split blocks by targeted ensemble
[BlockInfo, CellList, TrialData] = simple_Split_Blocks_by_MarkpointsExperiment(BlockInfo, CellList, TrialData, StimInfo, 'Analysis', Ops.markpoints_analysis);

%Clear blocks to save memory
clear('Blocks');

%% Optional: Check whether indvidiual trials are sound-responsive according to Z-score threshold and AUC_level

if Ops.check_single_trial_responses
    TrialData = simple_check_single_trial_responses(TrialData, StimInfo, Ops);
end

%% Optional: Behavior analysis

% if Ops.stim_protocol == 13
%     BehaviorData = extract_maryse_behavior(BlockInfo, TrialData, 0); %Ops.plot_extracted_data);
%     
%     %Save all open figures as FIG, PNG, and EPS
% %     if Ops.save_figures
% %         [~, ~, ~] = mkdir(Ops.save_path,'SummaryFigures');
% %         cd([Ops.save_path '/SummaryFigures'])
% %         save_plots(pwd, ['SummaryFig_' StimInfo.StimType '_' ExtractedData.Date])
% %         close all
% %     end
% end

if Ops.stim_protocol == 7 % Frequency Discrimination/associative learning
    [BehaviorData] = extract_FreqDisc_behavior(TrialData);
end

%% Loop through all cells/channels in each block

disp('Extracting data...')
oldBlock = '';
CellData = struct;
CellDataByStim = struct;
Traces = struct;
StimAnalysis = struct;

for c = 1:size(CellList,1)
    
    if c > c_end
        continue;
    end
    
    figure_names = {}; %Temporarily store figure titles
    figures = {}; %Temporarily store figures
    blockName = char(CellList.Block(c));
    b = find(strcmp({TrialData.Block}, blockName)); %Block index
    if ~strcmp(oldBlock, blockName)
        oldBlock = blockName;
        disp(oldBlock)
        
        %Make folders to save figures (one folder per block)
        if Ops.plot_figures && Ops.save_figures
            cd(Ops.save_path)
            [~, ~, ~] = mkdir(Ops.save_path,CellList.MouseID(c));
            [~, ~, ~] = mkdir(CellList.MouseID(c),blockName);
            for L = 1:length(Loco)
                [~, ~, ~] = mkdir(strcat(CellList.MouseID(c), '\', blockName),Loco{L});
            end
        end
    end
    cellNumber = CellList.CellNumber(c);
    cc = CellList.CellOrder(c);
    disp(cellNumber)
    
    %For Fibphot data, this skips blocks without functional data
    if isnan(cellNumber)
        for L = 1:length(Loco)
            %If cellNumber is the last cell in the loop, some tables won't have the last row(s) filled
            %Make sure they get filled here:
            CellData.(Loco{L}).RF(c) = 0;
            CellDataByStim.(Loco{L})(c).BlankTraces = [];
            StimAnalysis.(Loco{L})(c) = NaN;
        end
        continue;
    end
    
    %Get trial information from TrialData
    isBlank = TrialData(b).IsBlank;
    stim_v1 = TrialData(b).Stim_V1(~isBlank);
    stim_v2 = TrialData(b).Stim_V2(~isBlank);
    trials = squeeze(TrialData(b).Trials(cc,~isBlank,:)); %Final variable should be nTrials x nFrames
    blank_trials = squeeze(TrialData(b).Trials(cc,isBlank,:)); %Final variable should be nTrials x nFrames
    stim_isRunning = TrialData(b).IsRunning(~isBlank); %This will match the number of stim trials
    blank_isRunning = TrialData(b).IsRunning(isBlank); %This will match the number of blank trials
    
    %Reshape matrices if they only contain one trial
    if size(trials,2) == 1; trials = trials'; end
    if size(blank_trials,2) == 1; blank_trials = blank_trials'; end
    
    %% Check whether cell is responsive to each stimulus condition
    
    %Store average trace for each stim condition (filtered by locomotor activity)
    for L = 1:length(Loco)
        CellDataByStim.(Loco{L})(c).BlankTraces = [];
        CellDataByStim.(Loco{L})(c).StimTraces =  cell(size(StimInfo.Combined,1),1);
        CellDataByStim.(Loco{L})(c).StimTracesAveraged = single(nan(size(StimInfo.Combined,1),StimInfo.nFrames));
        CellDataByStim.(Loco{L})(c).nTrials = cast(zeros(size(StimInfo.Combined,1),1), 'int16');  %Caution, largest number you can represent with int16 is 32767 (intmax('int16'))
        
        for v = 1:size(StimInfo.Combined,1)
            v1 = StimInfo.Combined(v,1);
            v2 = StimInfo.Combined(v,2);
            
            %Find rows corresponding to current stim
            stim_rows = intersect(find(stim_v1 == v1), find(stim_v2 == v2));
            blank_rows = 1:size(blank_trials,1);
            
            %Filter by locomotor activity
            if strcmp(Loco{L}, 'Running')
                stim_rows = intersect(stim_rows, find(stim_isRunning == 1)); %Must compare to 0 and 1 because 9 is used for missing loco data
                blank_rows = intersect(blank_rows, find(blank_isRunning == 1));
            elseif strcmp(Loco{L}, 'NotRunning')
                stim_rows = intersect(stim_rows, find(stim_isRunning == 0)); %Must compare to 0 and 1 because 9 is used for missing loco data
                blank_rows = intersect(blank_rows, find(blank_isRunning == 0));
            end
            
            %Store StimTraces
            if ~isempty(stim_rows) %Some stim combinations might not exist in this block
                CellDataByStim.(Loco{L})(c).BlankTraces = single(blank_trials(blank_rows,:));
                CellDataByStim.(Loco{L})(c).nTrials(v,1) = size(trials(stim_rows,:),1);
                CellDataByStim.(Loco{L})(c).StimTraces{v} = single(trials(stim_rows,:)); %To be used with reliability
                CellDataByStim.(Loco{L})(c).StimTracesAveraged(v,:) = single(mean(trials(stim_rows,:),1,'omitnan'));
            end
        end
        
        %Determine if cell is significantly responsive based on average response to all trials
        figure_names{L} = {strcat('ResponseByStim_', Loco{L})};
        [CellDataByStim.(Loco{L})(c).PeakData, figures{L}] = simple_check_if_responsive(CellDataByStim.(Loco{L})(c).StimTracesAveraged, StimInfo.nBaselineFrames, StimInfo.fs, Ops.Z_level, Ops.AUC_level, Ops.smTime, Ops.plot_figures, 'MinimumPeakDuration', Ops.minimum_peak_duration, 'FigureDimensions', [length(StimInfo.V2),length(StimInfo.V1)], 'Subtitles', StimInfo.Combined_Label);
        %Store IsResponsive as RF value
        CellDataByStim.(Loco{L})(c).RF = CellDataByStim.(Loco{L})(c).PeakData.IsResponsive;
    end
    
    %% Check whether cell is reliable for each stim condition
    
    if Ops.use_reliability && Ops.reliability_type == "Pearsons"
        for L = 1:length(Loco)
            nResponsiveStim = sum(CellDataByStim.(Loco{L})(c).PeakData.IsResponsive == 1);
            for v = 1:size(StimInfo.Combined,1)
                %Compute reliability for all stim conditions but only compute
                %control on those that pass Z-score threshold (as this is the rate-limiting step)
                ReliabilityOps = Ops; %Duplicate Ops for use in simple_reliability function
                ReliabilityOps.ComputeReliabilityControls = CellDataByStim.(Loco{L})(c).PeakData.IsResponsive(v);
                
                %Adjust p-value for number of comparisons
                if Ops.adjust_pvalue
                    ReliabilityOps.reliability_p_value = Ops.reliability_p_value/nResponsiveStim;
                end
                
                
                %Compute reliability by comparing trials to either blank trials or shifted traces
                trials = CellDataByStim.(Loco{L})(c).StimTraces{v};
                blanks = CellDataByStim.(Loco{L})(c).BlankTraces; %If empty, script will automatically use shifted traces
                ReliabilityData = simple_reliability(trials, blanks, StimInfo, ReliabilityOps);
                
                %Concatenate results for each stim combination
                if v == 1
                    CellDataByStim.(Loco{L})(c).ReliabilityData = ReliabilityData;
                else
                    CellDataByStim.(Loco{L})(c).ReliabilityData = concatenateStruct(CellDataByStim.(Loco{L})(c).ReliabilityData, ReliabilityData);
                end
            end
            
            %Record whether cell is both reliable and responsive in RF
            IsResponsive = [CellDataByStim.(Loco{L})(c).PeakData.IsResponsive] == 1; %Convert NaNs to zeros
            IsReliable = [CellDataByStim.(Loco{L})(c).ReliabilityData.IsReliable]' == 1; %Convert NaNs to zeros
            CellDataByStim.(Loco{L})(c).RF = (IsResponsive + IsReliable) == 2; %Find overlap
            
            %Plot reliability figures
            if Ops.plot_figures
                figure_names{end+1} = strcat('Reliability_Traces_', Loco{L});
                figure_names{end+1} = strcat('Reliability_Blanks_', Loco{L});
                figure_names{end+1} = strcat('Reliability_Histograms_', Loco{L});
                [figures{end+1}, figures{end+2}, figures{end+3}] = plot_simple_reliability(CellDataByStim.(Loco{L})(c), StimInfo, Ops);
            end
        end
    end
    
    %% Compute reliability using xcorr
    
    if Ops.use_reliability && Ops.reliability_type == "XCorr"
        for L = 1:length(Loco)
            nResponsiveStim = sum(CellDataByStim.(Loco{L})(c).PeakData.IsResponsive == 1);
            for v = 1:size(StimInfo.Combined,1)
                %Only compute xcorr on stim conditions that pass Z-score threshold (use low Z_level for this reason)
                ReliabilityOps = Ops; %Duplicate Ops for use in simple_reliability function
                ReliabilityOps.ComputeReliabilityControls = CellDataByStim.(Loco{L})(c).PeakData.IsResponsive(v);
                
                %Adjust p-value for number of comparisons
                if ReliabilityOps.adjust_pvalue
                    ReliabilityOps.reliability_p_value = Ops.reliability_p_value/nResponsiveStim; 
                end
                
                %Compute reliability by comparing trials to either blank trials or shifted traces
                trials = CellDataByStim.(Loco{L})(c).StimTraces{v};
                blanks = CellDataByStim.(Loco{L})(c).BlankTraces;
                XcorrData = simple_xcorr(trials, blanks, StimInfo, v, ReliabilityOps);
                
                %Concatenate results for each stim combination
                if v == 1
                    CellDataByStim.(Loco{L})(c).XcorrData = XcorrData;
                else
                    CellDataByStim.(Loco{L})(c).XcorrData = [CellDataByStim.(Loco{L})(c).XcorrData; XcorrData];
                end
            end
            
            %Record whether cell is both reliable and responsive in RF
            IsResponsive = [CellDataByStim.(Loco{L})(c).PeakData.IsResponsive] == 1; %Convert NaNs to zeros
            IsReliable = CellDataByStim.(Loco{L})(c).XcorrData.cc_zero_z == 1; %Convert NaNs to zeros
            CellDataByStim.(Loco{L})(c).RF = (IsResponsive + IsReliable) == 2; %Find overlap
            
            %Plot xcorr figures
            if Ops.plot_figures
                figure_names{end+1} = strcat('Xcorr_Traces_', Loco{L});
                figure_names{end+1} = strcat('Xcorr_Blanks_', Loco{L});
                figure_names{end+1} = strcat('Xcorr_Max_Histograms_', Loco{L});
                figure_names{end+1} = strcat('Xcorr_Zero_Histograms_', Loco{L});
                [figures{end+1}, figures{end+2}, figures{end+3}, figures{end+4}] = plot_simple_xcorr(CellDataByStim.(Loco{L})(c), StimInfo, Ops);
            end
        end
    end
    
    %% Find neighbors of reliable/responsive stimuli within receptive field that are significantly correlated
    
    if Ops.use_neighbors
        for L = 1:length(Loco)
            %Reliability neighbors
            figure_names{end+1} = strcat('ReliabilityNeighbors_RF_', Loco{L});
            figure_names{end+1} = strcat('ReliabilityNeighbors_Traces_', Loco{L});
            
            if Ops.reliability_type == "Pearsons"
                %Pearson's R reliability
                [New_RF, NeighborData, figures{end+1}, figures{end+2}] = simple_reliability_neighbors(CellDataByStim.(Loco{L})(c), StimInfo, Ops.plot_figures, Ops);
            elseif Ops.reliability_type == "XCorr"
                %Xcorr reliability
                [New_RF, NeighborData, figures{end+1}, figures{end+2}] = simple_xcorr_neighbors(CellDataByStim.(Loco{L})(c), StimInfo, Ops.plot_figures, Ops);
            end
            
            %Store New_RF values
            CellDataByStim.(Loco{L})(c).RF = New_RF;
        end
    end
    
    %% Use the cell's receptive field responses for final trace
    
    for L = 1:length(Loco)
        T = table;
        
        %Store traces in array prefilled with nans
        if ~isfield(Traces, Loco{L})
            Traces.(Loco{L}) = single(nan(size(CellList,1),StimInfo.nFrames));
        end
        
        %If none of the stim conditions end up being significant, use all stim for average
        %Otherwise, only use significant conditions
        if ~any(CellDataByStim.(Loco{L})(c).RF)
            final_trace = mean(CellDataByStim.(Loco{L})(c).StimTracesAveraged,1,'omitnan');
            T.RF = 0; %If reliability was used, this will reflect reliability
            T.RF_Type = "none";
        else
            final_trace = mean(CellDataByStim.(Loco{L})(c).StimTracesAveraged(CellDataByStim.(Loco{L})(c).RF == 1,:),1,'omitnan');
            T.RF = 1; %If reliability was used, this will reflect reliability
            %Determine RF_type (excitatory, inhibitory, or mixed) using only significant conditions
            T.RF_Type = determine_RF_response_type(CellDataByStim.(Loco{L})(c).PeakData.ResponseType, CellDataByStim.(Loco{L})(c).RF);
        end
        
        %Compute response properties one more time for final trace.
        figure_names{end+1} = strcat('FinalTrace_', Loco{L});
        [PeakData, figures{end+1}] = simple_check_if_responsive(final_trace, StimInfo.nBaselineFrames, StimInfo.fs, Ops.Z_level, Ops.AUC_level, Ops.smTime, Ops.plot_figures, 'MinimumPeakDuration', Ops.minimum_peak_duration);
        
        %If no peak detected in final trace, force categorization by setting Z score, AUC level, and minimum peak duration to 0
        if PeakData.IsResponsive == 0
            [PeakData, figures{end}] = simple_check_if_responsive(final_trace, StimInfo.nBaselineFrames, StimInfo.fs, 0, 0, Ops.smTime, Ops.plot_figures, 'MinimumPeakDuration', 0);
        end
        
        %Remove IsResponsive column (this has no meaning since we forced categorization)
        PeakData.IsResponsive = [];
        
        %Make ResponseType none if cell was not found to be responsive and reliable
        %Here, ResponseType refers to that of the final average trace
        if T.RF == 0
            PeakData.ResponseType = "none";
        end
        
        % Add summary reliability measures to CellData
        ReliabilitySummary = simple_compute_ReliabilitySummary(CellDataByStim.(Loco{L})(c));
        
        %Store new data
        CellData.(Loco{L})(c,:) = [T, PeakData, ReliabilitySummary];
        
        %Store trace
        if ~isempty(final_trace)
            Traces.(Loco{L})(c,:) = single(final_trace);
        end
    end
    
    %% Compute sparseness
    SparsenessData = struct;
    BinarySparsenessData = struct;
    
    for L = 1:length(Loco)
        switch Ops.stim_protocol
            case {2 3 5 6 8 16 18} %These are the stim that have been tested so far, other stim could also work here as long as they have > 1 stim parameter

                %Compute tuning sparseness
                figure_names{end+1} = strcat('Sparseness_', Loco{L});
                [SparsenessData.(Loco{L}), figures{end+1}] = simple_compute_tuning_sparseness(CellDataByStim.(Loco{L})(c).PeakData, CellDataByStim.(Loco{L})(c).RF, StimInfo, Ops.plot_figures);
                
                %Compute again using just zeros and ones for RF
                figure_names{end+1} = strcat('BinarySparseness_', Loco{L});
                TempRF = single(CellDataByStim.(Loco{L})(c).RF);
                TempRF(CellDataByStim.(Loco{L})(c).nTrials == 0) = NaN; %Correct RF to include NaNs where there were no trials
                [BinarySparsenessData.(Loco{L}), figures{end+1}] = simple_compute_tuning_sparseness(CellDataByStim.(Loco{L})(c).PeakData, CellDataByStim.(Loco{L})(c).RF, StimInfo, Ops.plot_figures, 'RF', TempRF);
                BinarySparsenessData.(Loco{L}).Properties.VariableNames =  ["Binary_Sparseness","Binary_V1_Sparseness","Binary_V2_Sparseness","Binary_V1_Mean_Sparseness","Binary_V2_Mean_Sparseness"];
        end
    end
    
    %% Perform stim-specific receptive field analyses
    %Compute tuning, combine with sparseness data, and store results in StimAnalysis
    
    for L = 1:length(Loco)
        figure_names{end+1} = strcat(StimInfo.StimType, Loco{L});
        
        switch Ops.stim_protocol
            
            case 2 %RF
                [TuningData, figures{end+1}] = simple_compute_frequency_tuning(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), TuningData];
                
            case 3 %FM
                R = nan(size(CellDataByStim.(Loco{L})(c).RF));
                if Ops.use_reliability && Ops.reliability_type == "Pearsons"
                    R = [CellDataByStim.(Loco{L})(c).ReliabilityData.R]';
                end
                [FMData, figures{end+1}] = simple_compute_FM(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, R, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), FMData];
                
            case 5 %SAM
                [SAMData, figures{end+1}] = simple_compute_SAM(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), SAMData];
                
            case 6 %SAMfreq
                [SAMFreqData, figures{end+1}] = simple_compute_SAMfreq(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), SAMFreqData];
                
            case 7 %FreqDisc Behavior
                [FreqDiscAnalysis.(Loco{L})(c,:)]  = simple_compute_FreqDisc_behavior(Loco{L},c, TrialData, StimInfo,CellList);
                StimAnalysis = struct;
                
            case 8 %ABI Levels
               [TuningData, figures{end+1}] =  simple_compute_frequency_tuning_VA(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
%                 StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), TuningData];
                
            case 10
                %Compute intensity tuning for noise stim with more than one dB
               [StimAnalysis.(Loco{L})(c,:), figures{end+1}] = simple_compute_Noise(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                
            case 12 %SPONTANEOUS
                %Ca_detection pulls out spontaneous calcium transients and calculates their amplitude, frequency, and duration
                [~, StimAnalysis.(Loco{L})(c,:), figures{end+1}] = simple_Ca_detection_spontaneous(TrialData(b).Timestamp, TrialData(b).F7, TrialData(b).Full_Loco, StimInfo.fs, cc, Loco{L}, Ops.plot_figures);
                
            case 13 %MARYSE BEHAVIOR STIM
                [StimAnalysis.(Loco{L})(c,:), figures{end+1}] = simple_compute_tuning_maryse_behavior(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);

            case 18 %vocalizations
                [VocalData, figures{end+1}] = simple_compute_vocalization(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), VocalData];
            
            case 16 %ABI-SAM
                [ABISAMData, figures{end+1}] = simple_compute_SAM(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), ABISAMData];
            
            case 19 % Victor Behavior to do TODO 

            case 1018 %vocalizations_noise
                [VocalData, figures{end+1}] = simple_compute_vocalization(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, Ops.plot_figures);
                StimAnalysis.(Loco{L})(c,:) = [SparsenessData.(Loco{L}), BinarySparsenessData.(Loco{L}), VocalData];

        end
    end
    
    %% Plot receptive field heat maps using peak data
    
    if Ops.plot_figures
        for L = 1:length(Loco)
            figure_names{end+1} = strcat('RF_', Loco{L});
            [figures{end+1}] = simple_plot_PeakData(StimInfo, CellDataByStim.(Loco{L})(c).PeakData, 'RF', CellDataByStim.(Loco{L})(c).RF);
        end
    end
    
    %% Save figures
    
    if Ops.plot_figures && Ops.save_figures
        
        %Skip cells that were deemed non-responsive
        isResponsive = false;
        for L = 1:length(Loco)
            isResponsive = isResponsive || CellData.(Loco{L}).RF(c);
        end
        
        if Ops.save_only_active_cells && ~isResponsive
            close all
            continue
        end
        
        %Figure out how many figures there were not including replicates for Loco
        nFigures = size(figures,2)/length(Loco);
        figure_numbers = sort(repmat(1:nFigures, 1, size(figures,2)/nFigures));
        
        %Label each figure with block name and cell number and save
        for f = 1:size(figures,2)
            %Skip if figure is empty or has been closed
            if isempty(figures{1,f}) || ~isvalid(figures{1,f})
                continue
            end
            
            set(0, 'currentfigure', figures{1,f});
            sgtitle(strcat(regexprep(blockName(1,10:end),'_',' ','emptymatch'), ' Cell ', num2str(cellNumber))) %Replace _ with a space
            
            h_name = strcat(CellList.MouseID(c), '_CellRow', num2str(c), '_ID', num2str(cellNumber), '_Fig', num2str(figure_numbers(f)), '_', figure_names{f});
            
            if contains(figure_names{f}, '_All')
                save_folder = strcat(CellList.MouseID(c), '\', blockName,'\All\');
            elseif contains(figure_names{f}, '_Running')
                save_folder = strcat(CellList.MouseID(c), '\', blockName,'\Running\');
            elseif contains(figure_names{f}, '_NotRunning')
                save_folder = strcat(CellList.MouseID(c), '\', blockName,'\NotRunning\');
            end
            
            saveas(figures{1,f}, strcat(save_folder, h_name), 'fig');
            saveas(figures{1,f}, strcat(save_folder, h_name), 'jpg');
        end
        
        close all
    end
end

%% Save ExtractedData

ExtractedData = struct;
ExtractedData.Date = datestr(now,'yyyymmdd-HHMMSS');
ExtractedData.Ops = Ops;
ExtractedData.Ops.Info = Info;
ExtractedData.StimInfo = StimInfo;
ExtractedData.BlockInfo = BlockInfo;
ExtractedData.TrialData = TrialData;
ExtractedData.CellList = CellList;
ExtractedData.CellData = CellData;
ExtractedData.CellDataByStim = CellDataByStim;
ExtractedData.StimAnalysis = StimAnalysis;
ExtractedData.FinalTraces = Traces;
ExtractedData.Summary = make_simple_extracted_data_summary(ExtractedData); %Put main results into stats table
if exist('BehaviorData', 'var'); ExtractedData.BehaviorData = BehaviorData; end
if exist('FreqDiscAnalysis', 'var'); ExtractedData.FreqDiscAnalysis = FreqDiscAnalysis; end
if exist('Archive', 'var'); ExtractedData.Archive = Archive; end
filename = strcat('ExtractedData_', StimInfo.StimType, '_', ExtractedData.Date);
ExtractedData.Filename = filename;

if Ops.save_extracted_data
    disp('Saving ExtractedData...')
    cd(Ops.save_path)
    save([filename '.mat'], 'ExtractedData', '-v7.3');
    for L = 1:length(Loco)
        writetable(ExtractedData.Summary.(Loco{L}),[filename '_' Loco{L} '_Summary.csv'],'Delimiter',',');
    end
end

%% Plot summary figures

if Ops.preview_blocks
    visualize_extracted_block(ExtractedData);               %Plot sound-responsive cell averages for each block (Loco = All by default)
    if ExtractedData.StimInfo.StimProtocol == 17
        visualize_extracted_ephys_block(ExtractedData)
    end
end

if Ops.plot_extracted_data
    
    %These functions use simple_subset_ExtractedData to filter plotted cells by the folowing parameters.
    %You can adjust them here:
    PlotOps = struct;
    PlotOps.Loco                = 'All'; %All, Running, or NotRunning
    PlotOps.ResponseType        = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
    PlotOps.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
    PlotOps.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
    PlotOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1), put ~ in front to keep everything excep that FOV
    PlotOps.Group               = ''; %'' -> no filtering, or add group name here to only keep that group, put ~ in front to keep everything excep that Group
    PlotOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
    PlotOps.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    PlotOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    PlotOps.plotCorrelations	= 0; %0 = don't plot dependent variable correlations, 1 = plot
    PlotOps.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
    PlotOps.plotShadedErrorBars = 0; %0 = don't plot, 1 = plot shaded error bars (some functions only)
    PlotOps.FibPhotChannel      = 1; %1 = blue, 2 = green, 3 = red, 4 = all of them
    PlotOps.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress
    PlotOps.ByMouse             = 1; %Sort by unique mice. Currently only utilized in FreqDisc (stimProtocol =7)

    plot_simple_ExtractedData(ExtractedData, PlotOps);               %Rasters and pie charts for activated vs. prolonged. vs. suppressed cells
    plot_simple_ExtractedData_Reliability(ExtractedData, PlotOps);   %Scatter plots with reliability measures (function will skip if no reliability computed)
    plot_simple_ExtractedData_PeakData(ExtractedData, PlotOps);      %Peak data histograms and rasters
    plot_simple_ExtractedData_Controls(ExtractedData, PlotOps);      %Plot results vs. sex, age, auditory field, GCaMP type, etc
    
    %Stim-specific plots (function will skip if not the right data type)
    plot_simple_ExtractedData_RF(ExtractedData, PlotOps);
    plot_simple_ExtractedData_FM(ExtractedData, PlotOps);
    plot_simple_ExtractedData_SAM(ExtractedData, PlotOps);
    plot_simple_ExtractedData_SAMfreq(ExtractedData, PlotOps);
    plot_simple_ExtractedData_Noise(ExtractedData, PlotOps);
    plot_simple_ExtractedData_MaryseBehavior(ExtractedData, PlotOps);
    plot_simple_ExtractedData_AllBehavior(ExtractedData, PlotOps);
    plot_simple_ExtractedData_FreqDisc(ExtractedData, PlotOps);
    plot_simple_ExtractedData_ABI(ExtractedData, PlotOps);
    plot_simple_ExtractedData_ABISAM(ExtractedData, PlotOps);
    
    %Stim-specific plots split by activated, prolonged, suppressed
    plot_simple_ExtractedData_RF_APS(ExtractedData, PlotOps);
    plot_simple_ExtractedData_FM_APS(ExtractedData, PlotOps);
    plot_simple_ExtractedData_SAM_APS(ExtractedData, PlotOps);
    plot_simple_ExtractedData_SAMfreq_APS(ExtractedData, PlotOps);
    plot_simple_ExtractedData_ABISAM_APS(ExtractedData, PlotOps);
    
end

%Save all open figures as FIG, PNG, and EPS
if Ops.save_figures && (Ops.preview_blocks || Ops.plot_extracted_data)
    [~, ~, ~] = mkdir(Ops.save_path,'SummaryFigures');
    cd([Ops.save_path '/SummaryFigures'])
    save_plots(pwd, ['SummaryFig_' ExtractedData.StimInfo.StimType '_' ExtractedData.Date])
end

disp('Done extracting data :-)')
