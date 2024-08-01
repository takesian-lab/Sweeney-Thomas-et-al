%Helper script to analyze_single_trial_reliability
%We are interested in getting the peak response amplitudes for activated/prolonged cells
%and the trough responses for the suppressed cells

%TODO: Add xcorr vs. pearsons reliability mode

%THIS SCRIPT IS CURRENTLY HARDCODED FOR NOISEITI WITH 20 TRIALS AND 20 BLANKS
%Load NoiseITI ExtractedData

save_folder = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Figures\Figure 1\1D';
plot_figures = 1;
save_figures = 0;
RF_filter_type = 2; %1 = reponsive, 2 = responsive & reliable
write_table = 0;
recompute_R = 1;
recompute_type = 1; %1 for Pearsons, 2 for XCorr

T = 20; %Max number of trials

TrialData = ExtractedData.TrialData;
CellDataByStim = ExtractedData.CellDataByStim.All;
StimInfo = ExtractedData.StimInfo;
Ops = ExtractedData.Ops;

if recompute_R
    Summary = ExtractedData.Summary.All;
    
    if recompute_type == 1 %Pearsons
        Summary.R_shifted = nan(height(Summary),1);
        Summary.R_shifted_Z = nan(height(Summary),1);
        Summary.R_shifted_blank = nan(height(Summary),1);
        Summary.R_first5 = nan(height(Summary),1);
        Summary.R_first5_Z = nan(height(Summary),1);
        Summary.R_first5_blank = nan(height(Summary),1);
        Summary.R_first5_blank_Z = nan(height(Summary),1);
        Summary.R_second5 = nan(height(Summary),1);
        Summary.R_second5_Z = nan(height(Summary),1);
        Summary.R_second5_blank = nan(height(Summary),1);
        Summary.R_second5_blank_Z = nan(height(Summary),1);
        Summary.R_third5 = nan(height(Summary),1);
        Summary.R_third5_Z = nan(height(Summary),1);
        Summary.R_third5_blank = nan(height(Summary),1);
        Summary.R_third5_blank_Z = nan(height(Summary),1);
        Summary.R_last5 = nan(height(Summary),1);
        Summary.R_last5_Z = nan(height(Summary),1);
        Summary.R_last5_blank = nan(height(Summary),1);
        Summary.R_last5_blank_Z = nan(height(Summary),1);
        Summary.R_last15 = nan(height(Summary),1);
        Summary.R_last15_Z = nan(height(Summary),1);
        Summary.R_last15_blank = nan(height(Summary),1);
        Summary.R_last15_blank_Z = nan(height(Summary),1);
    else %Xcorr
        %TODO: Add xcorr vs. pearsons reliability mode

        %Summary.
    end
        
else
    %LOAD SUMMARY - This functionality not in place yet: As long as you
    %aren't using the Shifted reliability mode, this script doesn't take
    %too long **(Update: I think the Shifted mode should be fast now)**
end

cd(save_folder)

table_initialized = 0;
for c = 1:height(Summary)
    
    %Only include responsive cells
    if RF_filter_type == 1
        if CellDataByStim(c).PeakData.IsResponsive == 0 
            continue;
        end
    elseif RF_filter_type == 2
        if Summary.RF(c) == 0
            continue;
        end
    end
    
    disp(c);
    
    Block = Summary.Block(c);
    Group = Summary.Group(c);
    block_idx = strcmp(Block, {TrialData.Block});
    cellOrder = Summary.CellOrder(c);
    cellNumber = Summary.CellNumber(c);
    
    IsBlank = TrialData(block_idx).IsBlank;
    Trial_ResponseType = TrialData(block_idx).Trial_ResponseType(cellOrder,:);
    Trial_IsResponsive = TrialData(block_idx).Trial_IsResponsive(cellOrder,:);
    Trial_PeakAUC = TrialData(block_idx).Peak_AUC(cellOrder,:);
    Trial_TroughAUC = TrialData(block_idx).Trough_AUC(cellOrder,:); %store as positive
    Trial_Peak = TrialData(block_idx).Peak_3fr(cellOrder,:);
    nTrials = length(IsBlank);
    
    %Store order of trials with blanks separated
    Trial_Order = nan(size(IsBlank));
    Trial_Order(~IsBlank) = 1:sum(~IsBlank);
    Trial_Order(IsBlank) = 1:sum(IsBlank);
    
    %Combine AUCs based on response type
    AUC = Trial_PeakAUC;
    for t = 1:length(AUC)
        if strcmp(Trial_ResponseType(t), 'suppressed')
            AUC(t) = -Trial_TroughAUC(t); %Make negative here
        end
    end

    %Normalize to first trial (separately for blank and not blank)
    firstNotBlankInd = find(~IsBlank,1,'first');
    firstBlankInd = find(IsBlank,1,'first');
    
    Peak_NormToFirstTrial = nan(size(AUC));
    AUC_NormToFirstTrial = nan(size(AUC));
    PeakAUC_NormToFirstTrial = nan(size(AUC));
    TroughAUC_NormToFirstTrial = nan(size(AUC));
    
    %Raise Peak and AUC floor to 1 before normalizing so that there are no negatives (can't raise to 0 or else we would be dividing by 0s)
    if min(Trial_Peak) < 0
        Trial_Peak_pos = Trial_Peak + abs(min(Trial_Peak)) + 1;
    else
        Trial_Peak_pos = Trial_Peak;
    end
    
    if min(AUC) < 0
        AUC_pos = AUC + abs(min(AUC)) + 1;
    else
        AUC_pos = AUC;
    end
    
    %Not Blank
    Peak_NormToFirstTrial(~IsBlank) = Trial_Peak_pos(~IsBlank)./Trial_Peak_pos(firstNotBlankInd);
    AUC_NormToFirstTrial(~IsBlank) = AUC_pos(~IsBlank)./AUC_pos(firstNotBlankInd);
    PeakAUC_NormToFirstTrial(~IsBlank) = Trial_PeakAUC(~IsBlank)./Trial_PeakAUC(firstNotBlankInd);
    TroughAUC_NormToFirstTrial(~IsBlank) = -(Trial_TroughAUC(~IsBlank)./Trial_TroughAUC(firstNotBlankInd));
    
    %Blank
    Peak_NormToFirstTrial(IsBlank) = Trial_Peak_pos(IsBlank)./Trial_Peak_pos(firstBlankInd);
    AUC_NormToFirstTrial(IsBlank) = AUC_pos(IsBlank)./AUC_pos(firstBlankInd);
    PeakAUC_NormToFirstTrial(IsBlank) = Trial_PeakAUC(IsBlank)./Trial_PeakAUC(firstBlankInd);
    TroughAUC_NormToFirstTrial(IsBlank) = -(Trial_TroughAUC(IsBlank)./Trial_TroughAUC(firstBlankInd));
    
    %% Compute reliability
    
    Trials = squeeze(TrialData(block_idx).Trials(cellOrder,:,:));
    BlankTrials = Trials(IsBlank,:);
    StimTrials = Trials(~IsBlank,:);
    
    %Only compute reliability on cells that are responsive
    if recompute_R && CellDataByStim(c).PeakData.IsResponsive == 1 && recompute_type == 1 %Pearsons

        %Store original R based on shifted traces
        Summary.R_shifted(c) = CellDataByStim(c).ReliabilityData.R;
        Summary.R_shifted_Z(c) = CellDataByStim(c).ReliabilityData.IsReliable;
        
        %Compute R value on all 20 blank trials
        [ReliabilityData] = simple_reliability(BlankTrials, BlankTrials, StimInfo, ExtractedData.Ops);
        Summary.R_shifted_blank(c) = ReliabilityData.R;
        
        %Reliability based on blanks for first 5 trials
        [ReliabilityData] = simple_reliability(StimTrials(1:5,:), BlankTrials, StimInfo, ExtractedData.Ops);
        Summary.R_first5(c) = ReliabilityData.R;
        Summary.R_first5_Z(c) = ReliabilityData.IsReliable;
        [ReliabilityDataBlank] = simple_reliability(BlankTrials(1:5,:), BlankTrials(6:end,:), StimInfo, ExtractedData.Ops);
        Summary.R_first5_blank(c) = ReliabilityDataBlank.R;
        Summary.R_first5_blank_Z(c) = ReliabilityDataBlank.IsReliable;

        %Reliability based on blanks for second 5 trials
        [ReliabilityData_Second5] = simple_reliability(StimTrials(6:10,:), BlankTrials, StimInfo, ExtractedData.Ops);
        Summary.R_second5(c) = ReliabilityData_Second5.R;
        Summary.R_second5_Z(c) = ReliabilityData_Second5.IsReliable;
        remaining_blankidx = 1:size(BlankTrials,1);
        for i = 6:10
            remaining_blankidx(remaining_blankidx == i) = [];
        end
        [ReliabilityDataBlank_Second5] = simple_reliability(BlankTrials(6:10,:), BlankTrials(remaining_blankidx,:), StimInfo, ExtractedData.Ops);
        Summary.R_second5_blank(c) = ReliabilityDataBlank_Second5.R;
        Summary.R_second5_blank_Z(c) = ReliabilityDataBlank_Second5.IsReliable;
        
        %Reliability based on blanks for third 5 trials
        [ReliabilityData_Third5] = simple_reliability(StimTrials(11:15,:), BlankTrials, StimInfo, ExtractedData.Ops);
        Summary.R_third5(c) = ReliabilityData_Third5.R;
        Summary.R_third5_Z(c) = ReliabilityData_Third5.IsReliable;
        remaining_blankidx = 1:size(BlankTrials,1);
        for i = 11:15
            remaining_blankidx(remaining_blankidx == i) = [];
        end
        [ReliabilityDataBlank_Third5] = simple_reliability(BlankTrials(11:15,:), BlankTrials(remaining_blankidx,:), StimInfo, ExtractedData.Ops);
        Summary.R_third5_blank(c) = ReliabilityDataBlank_Third5.R;
        Summary.R_third5_blank_Z(c) = ReliabilityDataBlank_Third5.IsReliable;

        %Reliability for last 5 trials
        [ReliabilityData_Last5] = simple_reliability(StimTrials(end-4:end,:), BlankTrials, StimInfo, ExtractedData.Ops);
        Summary.R_last5(c) = ReliabilityData_Last5.R;
        Summary.R_last5_Z(c) = ReliabilityData_Last5.IsReliable;
        [ReliabilityDataBlank_Last5] = simple_reliability(BlankTrials(end-4:end,:), BlankTrials(1:end-5,:), StimInfo, ExtractedData.Ops);
        Summary.R_last5_blank(c) = ReliabilityDataBlank_Last5.R;
        Summary.R_last5_blank_Z(c) = ReliabilityDataBlank_Last5.IsReliable;
        
        %Reliability for last 15 trials
        [ReliabilityData_Last15] = simple_reliability(StimTrials(6:end,:), BlankTrials, StimInfo, ExtractedData.Ops);
        Summary.R_last15(c) = ReliabilityData_Last15.R;
        Summary.R_last15_Z(c) = ReliabilityData_Last15.IsReliable;
        [ReliabilityDataBlank_Last15] = simple_reliability(BlankTrials(6:end,:), BlankTrials, StimInfo, ExtractedData.Ops);
        Summary.R_last15_blank(c) = ReliabilityDataBlank_Last15.R;
        Summary.R_last15_blank_Z(c) = ReliabilityDataBlank_Last15.IsReliable;
    elseif recompute_R && CellDataByStim(c).PeakData.IsResponsive == 1 && recompute_type == 2 %XCorr
        
        %TODO: Add xcorr vs. pearsons reliability mode
    end
    
    %% Plot figure 
    if plot_figures
        h = figure;
       
        ylim0 = [min(min(Trials)) max(max(Trials))];
        ylim1 = [-1 max(Trial_Peak)]; 
        ylim2 = [min(Peak_NormToFirstTrial) max(Peak_NormToFirstTrial)];
        
        
        for f = 1:2
            if f == 1 %Stim
                trialsToPlot = StimTrials;
                peakToPlot = Trial_Peak(~IsBlank);
                norm_peakToPlot = Peak_NormToFirstTrial(~IsBlank);
                title1 = 'Stim trials';
                title2 = 'Normalized to first trial';
                %titleR1 = ['First 5 stim: R = ' num2str(round(Summary.R_first5(c),2)) ' Z = ' num2str(round(Summary.R_first5_Z(c),2))];
                %titleR2 = ['Last 5 stim: R = ' num2str(round(Summary.R_last5(c),2)) ' Z = ' num2str(round(Summary.R_last5_Z(c),2))];
                titleR1 = ['All stim: R = ' num2str(round(Summary.R_shifted(c),2)) ' Z = ' num2str(round(Summary.R_shifted_Z(c),2))];
                titleR2 = ['Last 15 stim: R = ' num2str(round(Summary.R_last15(c),2)) ' Z = ' num2str(round(Summary.R_last15_Z(c),2))];
                
            elseif f == 2 %Blank
                trialsToPlot = BlankTrials;
                peakToPlot = Trial_Peak(IsBlank);
                norm_peakToPlot = Peak_NormToFirstTrial(IsBlank);
                title1 = 'Blank trials';
                title2 = 'Normalized to first trial';
                %titleR1 = ['First 5 blanks: R = ' num2str(round(Summary.R_first5_blank(c),2)) ' Z = ' num2str(round(Summary.R_first5_blank_Z(c),2))];
                %titleR2 = ['Last 5 blanks: R = ' num2str(round(Summary.R_last5_blank(c),2)) ' Z = ' num2str(round(Summary.R_last5_blank_Z(c),2))];
                titleR1 = ['All blanks: R = ' num2str(Summary.R_shifted_blank(c),2) ' Z = N/A'];
                titleR2 = ['Last 15 blanks: R = ' num2str(round(Summary.R_last15_blank(c),2)) ' Z = ' num2str(round(Summary.R_last15_blank_Z(c),2))];
            end
                        
            %All responses
            for t = 1:size(trialsToPlot,1)
                subplot(4, T, T*(f-1) + t); hold on
                plot(smooth(trialsToPlot(t,:),5), 'Linewidth', 1.5)
                xlim([1 size(trialsToPlot,2)])
                ylim(ylim0)
                vline(StimInfo.nBaselineFrames) %nBaselineFrames
                if f == 1
                    title(['Trial ' num2str(t)])
                end
            end
  
            %Response by trial
            subplot(4, T, T*f + [21:25]); hold on
            scatter(1:length(peakToPlot), peakToPlot)
            plot(1:length(peakToPlot), peakToPlot)
            xlim([0.5 T+0.5])
            ylim(ylim1)
            title(title1)
            if f == 2; xlabel('Trial'); end

            %Normalized response by trial
            subplot(4, T, T*f + [26:30]); hold on
            scatter(1:length(peakToPlot), norm_peakToPlot)
            plot(1:length(peakToPlot), norm_peakToPlot)
            xlim([0.5 T+0.5])
            ylim(ylim2)
            hline(1)
            title(title2)
            if f == 2; xlabel('Trial'); end
            
            %Plot all traces
            subplot(4, T, T*f + [31:33]); hold on
            for t = 1:T
                plot(smooth(trialsToPlot(t,:),5))
            end
            plot(smooth(mean(trialsToPlot,1,'omitnan'),5), 'k', 'Linewidth', 2)
            xlim([1 size(trialsToPlot,2)])
            ylim(ylim0)
            vline(StimInfo.nBaselineFrames) %nBaselineFrames
            title(titleR1)
    
            subplot(4, T, T*f + [34 35]); hold on
            Pearsons1 = corrcoef(trialsToPlot(1:20,:)','Rows','pairwise'); 
            imagesc(Pearsons1);
            xlim([0.5 size(Pearsons1,1)+0.5])
            ylim([0.5 size(Pearsons1,1)+0.5])
            caxis([0 1])
            colormap(bluewhitered_TAKlab)
            set(gca,'XTick', 1:2:size(Pearsons1,1))
            set(gca,'YTick', 1:2:size(Pearsons1,1))
            
            %Plot last 15 traces
            subplot(4, T, T*f + [36:38]); hold on
            for t = (size(trialsToPlot,1)-14):size(trialsToPlot,1)
                plot(smooth(trialsToPlot(t,:),5))
            end
            plot(smooth(mean(trialsToPlot(end-14:end,:),1,'omitnan'),5), 'k', 'Linewidth', 2)
            xlim([1 size(trialsToPlot,2)])
            ylim(ylim0)
            vline(StimInfo.nBaselineFrames) %nBaselineFrames
            title(titleR2)
            
            subplot(4, T, T*f + [39 40]); hold on
            Pearsons2 = corrcoef(trialsToPlot(end-14:end,:)','Rows','pairwise'); 
            imagesc(Pearsons2);
            set(gca,'XTick', 1:2:size(Pearsons2,1))
            set(gca,'YTick', 1:2:size(Pearsons2,1))
            xlim([0.5 size(Pearsons2,1)+0.5])
            ylim([0.5 size(Pearsons2,1)+0.5])
            caxis([0 1])
            colormap(bluewhitered_TAKlab)
            if f == 2; xlabel('Trials'); end
        end
        
        sgtitle([regexprep(Block,'_',' ','emptymatch') ' Cell ' num2str(cellNumber)])
        
        if save_figures
            fullscreen(h);
            FigName = strcat(Group, '_CellRow', num2str(c));
            savefig(h, strcat(FigName, '.fig'));
            saveas(h,strcat(FigName, '.png'))
            close(h)
        end
    end
    
    %% Make Results Table
    Temp = repmat(Summary(c,:), nTrials, 1);
    Temp.Trial = (1:nTrials)';
    Temp.Trial_Separated_ByBlanks = Trial_Order;
    Temp.Trial_IsBlank = IsBlank;
    Temp.Trial_ResponseType = Trial_ResponseType';
    Temp.Trial_IsResponsive = Trial_IsResponsive';
    Temp.Trial_Peak_AUC = Trial_PeakAUC';
    Temp.Trial_Trough_AUC = Trial_TroughAUC';
    Temp.Trial_AUC = AUC';
    Temp.Trial_Peak_AUC_NormToFirstTrial = PeakAUC_NormToFirstTrial';
    Temp.Trial_Trough_AUC_NormToFirstTrial = TroughAUC_NormToFirstTrial';
    Temp.Trial_AUC_NormToFirstTrial = AUC_NormToFirstTrial';
    Temp.Trial_Peak = Trial_Peak';
    Temp.Trial_Peak_NormToFirstTrial = Peak_NormToFirstTrial';
    
    if table_initialized == 0
        Results = Temp;
        table_initialized = 1;
    else
        Results = [Results; Temp];
    end
end

%% Write table

if write_table
    writetable(Results,'Noise_Single_Trial_Stats.csv','Delimiter',',');
    save('Results.mat', 'Results');
    writetable(Summary,'Noise_Summary.csv','Delimiter',',');
    save('Summary.mat', 'Summary');
end

%% To update ExtractedData with newly computed R value
%Use last 5 trials to avoid any effect of adaptation
return;

newR = Summary.R_last15_Z == 1;
ExtractedData.Summary.All.RF = newR;
ExtractedData.CellData.All.RF = newR;
for c = 1:length(newR)
    ExtractedData.CellDataByStim.All(c).RF = newR(c);
    ExtractedData.Summary.All.ResponseType(c) = ExtractedData.CellDataByStim.All(c).PeakData.ResponseType;
    ExtractedData.CellData.All.ResponseType(c) =  ExtractedData.CellDataByStim.All(c).PeakData.ResponseType;
    %[PeakData,~] = simple_check_if_responsive(ExtractedData.FinalTraces.All(c,:), StimInfo.nBaselineFrames, StimInfo.fs, Ops.Z_level, Ops.AUC_level, Ops.smTime, Ops.plot_figures);
end

save([ExtractedData.Filename(1:end-4) '_updatedReliabilityAllOnly.mat'],'ExtractedData', '-v7.3')