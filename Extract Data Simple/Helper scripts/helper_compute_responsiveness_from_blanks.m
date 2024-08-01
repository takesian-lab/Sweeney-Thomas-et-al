%helper_compute_responsiveness_from_blanks

%% Q1: How spontaneously active are our cell populations NDNF vs PYR vs VIP?
% What is the percentage of blank trials (not averaged) on which we find a
% peak or trough? and do those peaks/troughs differ?

%% Q2: How likely is our analysis to detect realiable and responsive peaks from blanks?
%For every cell, subsample 5 blank trials and compute whether the average of those 5 trials is responsive (has a peak)
%and reliable (compared to the rest of the blank trials). Do this 1000 times subsampling different blank trials.
%The average of the 1000 yes/no responses is the “percent likelihood that we would have found a responsive and reliable peak from the blank trials”
%Three measures:
%1. Percent likelihood of detecting a peak (or trough)
%2. Percent likelihood of detecting a reliable peak (or trough)
%3. Out of all peaks detected, percent likelihood they would be reliable 

%%

%load extracted data
clearvars -except ExtractedData

nBootstraps = 1; %n bootstraps

Ops = ExtractedData.Ops;
StimInfo = ExtractedData.StimInfo;

%Preallocate for everything here (but I haven't pulled out data for stim
%peaks because how do we deal with sounds that cell is not responsive to?
[PercentResponsive, PercentReliableOfTotal, PercentReliableOfResponsive,... %Q2
    PercentBlankPeaks_A, PercentBlankPeaks_P, PercentBlankPeaks_S, PercentBlankPeaks_AP,... %Q1
    PercentStimPeaks_A, PercentStimPeaks_P, PercentStimPeaks_S, PercentStimPeaks_AP,...
    MeanBlankPeaks_A, STDBlankPeaks_A, MeanStimPeaks_A, STDStimPeaks_A,...
    MeanBlankPeaks_P, STDBlankPeaks_P, MeanStimPeaks_P, STDStimPeaks_P,...
    MeanBlankPeaks_S, STDBlankPeaks_S, MeanStimPeaks_S, STDStimPeaks_S] = ...
    deal(nan(size(ExtractedData.CellDataByStim.All,2),1));

for c = 1:size(ExtractedData.CellDataByStim.All,2)
    disp(c);
    
    nTrials = 5; %mode(ExtractedData.CellDataByStim.All(c).nTrials);
    BlankTraces = ExtractedData.CellDataByStim.All(c).BlankTraces;

    %Get peak response properties from blanks Q1
    BlankPeakData = simple_check_if_responsive(BlankTraces,StimInfo.nBaselineFrames,StimInfo.fs,Ops.Z_level,Ops.AUC_level, 0, 0);
    
    PercentBlankPeaks_A(c) = sum(BlankPeakData.IsResponsive(strcmp(BlankPeakData.ResponseType,'activated')))/size(BlankPeakData,1);
    PercentBlankPeaks_P(c) = sum(BlankPeakData.IsResponsive(strcmp(BlankPeakData.ResponseType,'prolonged')))/size(BlankPeakData,1);
    PercentBlankPeaks_S(c) = sum(BlankPeakData.IsResponsive(strcmp(BlankPeakData.ResponseType,'suppressed')))/size(BlankPeakData,1);
    PercentBlankPeaks_AP(c) = (sum(BlankPeakData.IsResponsive(strcmp(BlankPeakData.ResponseType,'activated')))...
        + sum(BlankPeakData.IsResponsive(strcmp(BlankPeakData.ResponseType,'prolonged'))))/size(BlankPeakData,1);
    
    MeanBlankPeaks_A(c) = mean(BlankPeakData.Peak_AUC(strcmp(BlankPeakData.ResponseType,'activated')),'omitnan');
    STDBlankPeaks_A(c) = std(BlankPeakData.Peak_AUC(strcmp(BlankPeakData.ResponseType,'activated')),'omitnan');
    MeanBlankPeaks_P(c) = mean(BlankPeakData.Peak_AUC(strcmp(BlankPeakData.ResponseType,'prolonged')),'omitnan');
    STDBlankPeaks_P(c) = std(BlankPeakData.Peak_AUC(strcmp(BlankPeakData.ResponseType,'prolonged')),'omitnan');
    MeanBlankPeaks_S(c) = mean(BlankPeakData.Trough_AUC(strcmp(BlankPeakData.ResponseType,'suppressed')),'omitnan');
    STDBlankPeaks_S(c) = std(BlankPeakData.Trough_AUC(strcmp(BlankPeakData.ResponseType,'suppressed')),'omitnan');
    
    if isempty(BlankTraces)
        continue;
    end
    
    %Figure out how likely it would have been for us to call this a sound response Q2
    if size(BlankTraces,1) <= nTrials
        %If we didn't have more blanks than trials (e.g. NoiseITI)....
        Traces = BlankTraces;
        TracesAveraged = mean(Traces,1,'omitnan');
        [PeakData, ~] = simple_check_if_responsive(TracesAveraged, StimInfo.nBaselineFrames, StimInfo.fs, Ops.Z_level, Ops.AUC_level, Ops.smTime, Ops.plot_figures);
        PercentResponsive(c) = PeakData.IsResponsive;
        
        %Compute reliability for significant peaks only
        %This will use the circle-shifted method since Traces and BlankTraces are the same size
        if PeakData.IsResponsive
            [ReliabilityData] = simple_reliability(Traces, BlankTraces, StimInfo, Ops);
            PercentReliableOfTotal(c) = ReliabilityData.IsReliable/PeakData.IsResponsive;
        end

    else
        %If we did... subsample nTrials nBootstraps times from BlankTraces to see what % of time we would have found a significant peak
        IsResponsive = nan(nBootstraps,1);
        IsReliable = nan(nBootstraps,1);
        for r = 1:nBootstraps
            ind = randperm(size(BlankTraces,1), nTrials); %Choose n blank trials at random without replacement
            Traces = BlankTraces(ind,:);
            TracesAveraged = mean(Traces,1,'omitnan');
            [PeakData, ~] = simple_check_if_responsive(TracesAveraged, StimInfo.nBaselineFrames, StimInfo.fs, Ops.Z_level, Ops.AUC_level, Ops.smTime, Ops.plot_figures);
            IsResponsive(r) = PeakData.IsResponsive;

            %Compute reliability for significant peaks only
            %Use the rest of the blanks as a fill-in for blank trials
            if PeakData.IsResponsive
                blankind = setdiff(1:size(BlankTraces,1), ind);
                [ReliabilityData] = simple_reliability(Traces, BlankTraces(blankind,:), StimInfo, Ops);
                IsReliable(r) = ReliabilityData.IsReliable;
            end
        end

        PercentResponsive(c) = sum(IsResponsive)/nBootstraps;
        PercentReliableOfTotal(c) = sum(IsReliable,'omitnan')/nBootstraps;
        PercentReliableOfResponsive(c) = sum(IsReliable,'omitnan')/sum(IsResponsive);
    end
end

return % MAKE USER INTENTIONALLY PLOT FIGURES AND SAVE

%% Plot figure Q1

Summary = ExtractedData.Summary.All;

%Groups = unique(Summary.Group);
Groups = ["PYR", "VIP", "NDNF"]; %Reorder

BlankPercentData_A = nan(2,length(Groups));
BlankPercentData_P = nan(2,length(Groups));
BlankPercentData_S = nan(2,length(Groups));
BlankPercentData_AP = nan(2,length(Groups));
BlankMeanData_A = nan(2,length(Groups));
BlankMeanData_P = nan(2,length(Groups));
BlankMeanData_S = nan(2,length(Groups));
BlankSTDData_A = nan(2,length(Groups));
BlankSTDData_P = nan(2,length(Groups));
BlankSTDData_S = nan(2,length(Groups));

for G = 1:length(Groups)
    ind = strcmp(Summary.Group, Groups(G));
    N = sum(ind);
    
    BlankPercentData_A(1,G) = mean(PercentBlankPeaks_A(ind),'omitnan');
    BlankPercentData_P(1,G) = mean(PercentBlankPeaks_P(ind),'omitnan');
    BlankPercentData_S(1,G) = mean(PercentBlankPeaks_S(ind),'omitnan');
    BlankPercentData_AP(1,G) = mean(PercentBlankPeaks_AP(ind),'omitnan');
    BlankPercentData_A(2,G) = std(PercentBlankPeaks_A(ind),'omitnan')/sqrt(N-1);
    BlankPercentData_P(2,G) = std(PercentBlankPeaks_P(ind),'omitnan')/sqrt(N-1);
    BlankPercentData_S(2,G) = std(PercentBlankPeaks_S(ind),'omitnan')/sqrt(N-1);
    BlankPercentData_AP(2,G) = std(PercentBlankPeaks_AP(ind),'omitnan')/sqrt(N-1);
    BlankMeanData_A(1,G) = mean(MeanBlankPeaks_A(ind),'omitnan');
    BlankMeanData_A(2,G) = std(MeanBlankPeaks_A(ind),'omitnan')/sqrt(N-1);
    BlankMeanData_P(1,G) = mean(MeanBlankPeaks_P(ind),'omitnan');
    BlankMeanData_P(2,G) = std(MeanBlankPeaks_P(ind),'omitnan')/sqrt(N-1);
    BlankMeanData_S(1,G) = mean(MeanBlankPeaks_S(ind),'omitnan');
    BlankMeanData_S(2,G) = std(MeanBlankPeaks_S(ind),'omitnan')/sqrt(N-1);
    BlankSTDData_A(1,G) = mean(STDBlankPeaks_A(ind),'omitnan');
    BlankSTDData_A(2,G) = std(STDBlankPeaks_A(ind),'omitnan')/sqrt(N-1);
    BlankSTDData_P(1,G) = mean(STDBlankPeaks_P(ind),'omitnan');
    BlankSTDData_P(2,G) = std(STDBlankPeaks_P(ind),'omitnan')/sqrt(N-1);
    BlankSTDData_S(1,G) = mean(STDBlankPeaks_S(ind),'omitnan');
    BlankSTDData_S(2,G) = std(STDBlankPeaks_S(ind),'omitnan')/sqrt(N-1);
end

figure; hold on

subplot(4,3,1); hold on
bar(1:length(Groups), BlankMeanData_A(1,:))
errorbar(1:length(Groups), BlankMeanData_A(1,:), BlankMeanData_A(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Activated')
ylabel('Mean')

subplot(4,3,4); hold on
bar(1:length(Groups), BlankSTDData_A(1,:))
errorbar(1:length(Groups), BlankSTDData_A(1,:), BlankSTDData_A(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
ylabel('STD')

subplot(4,3,2); hold on
bar(1:length(Groups), BlankMeanData_P(1,:))
errorbar(1:length(Groups), BlankMeanData_P(1,:), BlankMeanData_P(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Prolonged')
ylabel('Mean')

subplot(4,3,5); hold on
bar(1:length(Groups), BlankSTDData_P(1,:))
errorbar(1:length(Groups), BlankSTDData_P(1,:), BlankSTDData_P(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
ylabel('STD')

subplot(4,3,3); hold on
bar(1:length(Groups), BlankMeanData_S(1,:))
errorbar(1:length(Groups), BlankMeanData_S(1,:), BlankMeanData_S(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Suppressed')
ylabel('Mean')

subplot(4,3,6); hold on
bar(1:length(Groups), BlankSTDData_S(1,:))
errorbar(1:length(Groups), BlankSTDData_S(1,:), BlankSTDData_S(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
%title('Activated')
ylabel('STD')

subplot(4,3,7); hold on
bar(1:length(Groups), BlankPercentData_A(1,:))
errorbar(1:length(Groups), BlankPercentData_A(1,:), BlankPercentData_A(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
%title('Activated')
ylabel('Percent of trials')

subplot(4,3,8); hold on
bar(1:length(Groups), BlankPercentData_P(1,:))
errorbar(1:length(Groups), BlankPercentData_P(1,:), BlankPercentData_P(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
ylabel('Percent of trials')

subplot(4,3,9); hold on
bar(1:length(Groups), BlankPercentData_S(1,:))
errorbar(1:length(Groups), BlankPercentData_S(1,:), BlankPercentData_S(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
ylabel('Percent of trials')

subplot(4,3,10); hold on
bar(1:length(Groups), BlankPercentData_AP(1,:))
errorbar(1:length(Groups), BlankPercentData_AP(1,:), BlankPercentData_AP(2,:))
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Activated & Prolonged')
ylabel('Percent of trials')

sgtitle(StimInfo.StimType)

%% Plot figure Q2

Summary = ExtractedData.Summary.All;

%Groups = unique(Summary.Group);
Groups = ["PYR", "VIP", "NDNF"]; %Reorder

%Summarize data
ResponsiveData = cell(1,length(Groups));
ResponsiveMean = nan(1,length(Groups));
ReliableData = cell(1,length(Groups));
ReliableMean = nan(1,length(Groups));
ReliableOfResponsiveData = cell(1,length(Groups));
ReliableOfResponsiveMean = nan(1,length(Groups));

BlankResponsiveData = cell(1,length(Groups));
BlankResponsiveMean = nan(1,length(Groups));
BlankReliableData = cell(1,length(Groups));
BlankReliableMean = nan(1,length(Groups));
BlankReliableOfResponsiveData = cell(1,length(Groups));
BlankReliableOfResponsiveMean = nan(1,length(Groups));

for G = 1:length(Groups)
    ind = strcmp(Summary.Group, Groups(G));
    N = sum(ind);
    ResponsiveData{G} = ~isnan(Summary.R1(ind)); %If data has an R1 value, that means it had a peak
    ResponsiveMean(G) = sum(ResponsiveData{G})./N;
    ReliableData{G} = Summary.RF(ind);
    ReliableMean(G) = sum(ReliableData{G})./N;
    ReliableOfResponsiveData{G} =  ReliableData{G}(ResponsiveData{G} == 1);
    ReliableOfResponsiveMean(G) = sum(ReliableOfResponsiveData{G})./length(ReliableOfResponsiveData{G});
    
    BlankResponsiveData{G} = PercentResponsive(ind);
    BlankReliableData{G} = PercentReliableOfTotal(ind);
    
    if size(BlankTraces,1) <= nTrials
        BlankReliableOfResponsiveData{G} = BlankReliableData{G}(BlankResponsiveData{G} == 1);
        BlankResponsiveMean(G) = sum(BlankResponsiveData{G})./N;
        BlankReliableMean(G) = sum(BlankReliableData{G},'omitnan')./N;
        BlankReliableOfResponsiveMean(G) = sum(BlankReliableOfResponsiveData{G},'omitnan')./length(BlankReliableOfResponsiveData{G});
    else
        BlankReliableOfResponsiveData{G} = PercentReliableOfResponsive(ind);
        BlankResponsiveMean(G) = mean(BlankResponsiveData{G},'omitnan');
        BlankReliableMean(G) = mean(BlankReliableData{G},'omitnan');
        BlankReliableOfResponsiveMean(G) = mean(BlankReliableOfResponsiveData{G},'omitnan');
    end
end

figure; hold on

subplot(2,3,1)
bar(ResponsiveMean)
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Sound responsive')
ylabel('%')
ylim([0 1])

subplot(2,3,2)
bar(ReliableMean)
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Sound R&R')
ylim([0 1])

subplot(2,3,3)
bar(ReliableOfResponsiveMean)
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Sound Rel out of Resp')
ylim([0 1])

subplot(2,3,4)
bar(BlankResponsiveMean)
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Blank responsive')
ylabel('%')
ylim([0 1])

subplot(2,3,5)
bar(BlankReliableMean)
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Blank R&R')
ylim([0 1])

subplot(2,3,6)
bar(BlankReliableOfResponsiveMean)
set(gca,'XTick', 1:length(Groups))
set(gca,'XTickLabel', Groups)
title('Blank Rel out of Resp')
ylim([0 1])

sgtitle(StimInfo.StimType)

%% Store in ExtractedData, remove  and resave summary

return; %Make user intentionally save data

%%
%ExtractedData.Summary.All.EstimateBlankPercentResponsive = PercentResponsive;
%ExtractedData.Summary.All.EstimateBlankPercentReliable = PercentReliableOfTotal;
%ExtractedData.Summary.All.EstimateBlankPercentReliableOfResponsive = PercentReliableOfResponsive;
ExtractedData.Summary.All.PercentBlankPeaks_A = PercentBlankPeaks_A;
ExtractedData.Summary.All.PercentBlankPeaks_P = PercentBlankPeaks_P;
ExtractedData.Summary.All.PercentBlankPeaks_S = PercentBlankPeaks_S;
ExtractedData.Summary.All.PercentBlankPeaks_AP = PercentBlankPeaks_AP;
ExtractedData.Summary.All.MeanBlankPeaks_A = MeanBlankPeaks_A;
ExtractedData.Summary.All.MeanBlankPeaks_P = MeanBlankPeaks_P;
ExtractedData.Summary.All.MeanBlankPeaks_S = MeanBlankPeaks_S;
ExtractedData.Summary.All.STDBlankPeaks_A = STDBlankPeaks_A;
ExtractedData.Summary.All.STDBlankPeaks_P = STDBlankPeaks_P;
ExtractedData.Summary.All.STDBlankPeaks_S = STDBlankPeaks_S;

Summary = ExtractedData.Summary.All;
filename = ExtractedData.Filename;

%save([filename(1:end-4) '_withBlankResponsiveness.mat'], 'ExtractedData', '-v7.3');
save([filename(1:end-4) '_withBlankResponsiveness_All_Summary.mat'], 'Summary');
writetable(ExtractedData.Summary.All,[filename(1:end-4) '_withBlankResponsiveness_All_Summary.csv'],'Delimiter',',');