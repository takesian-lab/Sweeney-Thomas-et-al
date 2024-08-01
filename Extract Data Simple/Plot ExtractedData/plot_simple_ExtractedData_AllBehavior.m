function plot_simple_ExtractedData_AllBehavior(ExtractedData, UserOps)
% Plot behavior data (loco, whisker, pupil, etc)

%% Options

sortRasters = 1; % 0 for rasters to be in trial order, 1 to be sorted individually by behavior trace
sortBy = 'PeakAvg'; %'PeakAvg', 'StimLength' What value to sort traces by
zScoreBehaviorTraces = 1; % 0 to use raw behavior traces, 1 to z-score to baseline

%Return if not behavior data
if ~ismember(ExtractedData.StimInfo.StimProtocol, [7 13 9 50]) %Freq disc behavior, Maryse behavior, water, loco
    return;
end

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%I RECOMMEND NOT CHANGING THESE AND USING USEROPS INSTEAD
if nargin < 2
    UserOps = struct;
    UserOps.Loco                = 'All'; %All, Running, or NotRunning
    UserOps.ResponseType        = 'activated'; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
    UserOps.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
    UserOps.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
    UserOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1)
    UserOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
    UserOps.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    UserOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    UserOps.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
    UserOps.FibPhotChannel      = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
    UserOps.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress
end

disp('Plotting behavior data')

SubsetData = simple_subset_ExtractedData(ExtractedData, UserOps);

Ops = SubsetData.SubsetOps;
BlockInfo = SubsetData.BlockInfo;
TrialData = SubsetData.TrialData;
StimInfo = SubsetData.StimInfo;
Summary = SubsetData.Summary;
GroupList = Summary.Group;
Groups = unique(GroupList);

if isempty(TrialData)
    disp('No behavior data found to plot')
    return
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

%Useful for plotting
peakWin = 0.5; %Window for sorting in seconds
peakWinInFrames = peakWin*StimInfo.fs;
nBaselineFrames = StimInfo.nBaselineFrames;
xInFrames = nBaselineFrames + StimInfo.after_inFrames;
xInSeconds = xInFrames/StimInfo.fs;
graphx = linspace(0, xInSeconds, xInFrames); %changes trial length (x) in 180 frames to 6 seconds
xtick = 0:nBaselineFrames:xInFrames;
xticklabels = xtick./StimInfo.fs;

%% Concatenate behavioral variables from TrialData

variablesToConcat = {'Stim_V1', 'Stim_V2', 'Stim_Length', 'IsBlank', 'Loco', 'IsRunning', 'Licks', 'Pupil', 'Whisker',...
    'Outcomes', 'ReactionTime', 'HoldingPeriod', 'WaitPeriod', 'WaitCondition', 'Result', 'WaterDelivery'};

%Remove variables not in table (e.g. pupil might not be included)
variablesToConcat(~ismember(variablesToConcat,fields(TrialData))) = [];

%Remove columns from BlockInfo that we won't need
columnsToRemove = ["Framerate", "Baseline", "After_Stim", "conv_factorX", "conv_factorY", "refImg", "meanImgE"];
for c = 1:length(columnsToRemove)
    if ismember(columnsToRemove(c), BlockInfo.Properties.VariableNames)
        BlockInfo = removevars(BlockInfo, columnsToRemove(c));
    end
end

%Concatenate and store trial data
ConcatData = struct;

for v = 1:length(variablesToConcat)
    concat_variable = [];
    if v == 1; concat_BlockInfo = table; end
    for b = 1:size(TrialData,2)
        
        temp = TrialData(b).(variablesToConcat{v});
        
        %SPECIAL CASES
        if isempty(temp)
            error('Need to code for this scenario')
            %Data might be empty of there was no pupil or whisker trace for example
            %Need to make temp a matrix filled with nans that is nTrials x nFrames
        end
        
        %Make licks a full matrix
        if strcmp(variablesToConcat{v}, 'Licks')
            temp = full(temp);
        elseif any(strcmp(variablesToConcat{v}, {'Loco', 'Pupil', 'Whisker'}))           
            %Zscore all traces except for licks
            if zScoreBehaviorTraces
                baseline = temp(:,1:nBaselineFrames);
                std_baseline = std(baseline,[],2,'omitnan');
                std_baseline(std_baseline == 0) = 1; %Correct for all zero baselines to avoid dividing by 0
                %If baseline is all zeros, the trace will be (temp - 0)/1 = temp: so we are just keeping the raw trace
                zscore_temp = (temp - mean(baseline,2,'omitnan'))./std_baseline;
                temp = zscore_temp;
            end
        elseif size(temp,1) == 1
            %Make sure 1D variables are vertical
            temp = temp';
        end
        
        %CONCATENATE
        concat_variable = [concat_variable; temp];
        
        %ALSO SAVE BLOCK INFO
        if v == 1; concat_BlockInfo = [concat_BlockInfo; repmat(BlockInfo(b,:),size(temp,1),1)]; end
    end
    if v == 1; ConcatData.BlockInfo = concat_BlockInfo; end
    ConcatData.(variablesToConcat{v}) = concat_variable;
end

%% Separate trials into outcomes

switch StimInfo.StimProtocol
    
    case 7 %Freq disc behavior
        outcome = ConcatData.Stim_V2;
        result = strings(size(outcome));
        result(outcome == 1) = "Hit";
        result(outcome == 0) = "Miss";
        result(outcome == 3) = "Withhold";
        result(outcome == 4) = "False Alarm";
        result(outcome == 5) = "Catch";
        
    case 9 %Random water
        outcome = ConcatData.Stim_V2;
        result = strings(size(outcome));
        result(outcome == 1) = "Hit";
        result(outcome == 0) = "Miss";
        result(ConcatData.IsBlank == 1) = "Withhold";
        result(outcome == 4) = "False Alarm";
        
    case 13 %Maryse behavior
        result = ConcatData.Result;
        catch_index = eq(ConcatData.WaterDelivery, 0);
        result(catch_index) = "Catch";

    otherwise %All non-behavioral data
        result = strings(size(ConcatData.Stim_V1));
        result(:) = "Trial";
        
end

unique_results = unique(result);

%% Plot behavior
% There will be one column for each result (hit, miss, etc) and two rows for each behavior measure (licks, pupil, etc)

behaviorToPlot = {'Pupil', 'Whisker', 'Loco', 'Licks'};

%Remove behaviors that are not in ConcatData
behaviorToPlot(~ismember(behaviorToPlot,fields(ConcatData))) = [];

for G = 1:length(Groups) %One figure per group

    GroupRows = strcmp(ConcatData.BlockInfo.Group, Groups{G});
    
    h = figure; t = tiledlayout(length(behaviorToPlot)*2, length(unique_results));

    minMat = nan(length(behaviorToPlot),length(unique_results));
    maxMat = nan(length(behaviorToPlot),length(unique_results));
    count = 1;

    for r = 1:length(unique_results)
        for b = 1:length(behaviorToPlot)

            isResult = strcmp(result(GroupRows),unique_results{r});
            dataToPlot = ConcatData.(behaviorToPlot{b})(GroupRows,:);
            dataToPlot = dataToPlot(isResult,:);
            stimLength = ConcatData.Stim_Length(GroupRows);
            stimLength = stimLength(isResult);

            %Remove traces with nans
            nanrows = sum(isnan(dataToPlot),2) > 0;
            dataToPlot(nanrows,:) = [];
            stimLength(nanrows,:) = [];

            %Average
            ax(count) = nexttile(t, r + (b-1)*2*length(unique_results)); hold on
            y = mean(dataToPlot,1,'omitnan');
            minMat(b,r) = min(y);
            maxMat(b,r) = max(y);
            count = count + 1;

            plot(graphx, y)
            ylabel(behaviorToPlot{b})
            xlim([graphx(1) graphx(end)])
            if b == 1
                title(unique_results{r})
            end

            %Imagesc
            im(count) = nexttile(t, r + (b-1)*2*length(unique_results) + length(unique_results)); hold on
            if UserOps.smooth_rasters
                for d = 1:size(dataToPlot,1)
                    dataToPlot(d,:) = smooth(dataToPlot(d,:));
                end
            end
            if sortRasters
                switch sortBy
                    case 'StimLength'
                    [~, sortInd] = sort(stimLength, 'descend');
                    
                    case 'PeakAvg'
                    peakavg = mean(dataToPlot(:,nBaselineFrames+peakWinInFrames:nBaselineFrames+peakWinInFrames),2,'omitnan');
                    [~, sortInd] = sort(peakavg, 'ascend');
                end
                dataToPlot = dataToPlot(sortInd,:);
            end
            imagesc(dataToPlot)
            set(gca, 'XTick', xtick)
            set(gca, 'XTickLabel', xticklabels)
            vline(nBaselineFrames)
            ylabel('Trials')
            try
                xlim([1 size(dataToPlot,2)])
                ylim([1 size(dataToPlot,1)]);
            catch
                %nothin
            end
            if b == length(behaviorToPlot)
                xlabel('Seconds')
            end
        end
    end

    %Set y max and min to be the same for all plots
    minnestMat = min(minMat,[],2);
    maxestMat = max(maxMat,[],2);
    count = 1;
    for r = 1:length(unique_results)
        for b = 1:length(behaviorToPlot)
            try
                set(ax(count),'Ylim', [minnestMat(b,1) maxestMat(b,1)]);
            catch
                %do nothing
            end
            set(h, 'CurrentAxes', ax(count));
            vline(nBaselineFrames/StimInfo.fs);
            set(h, 'CurrentAxes', im(count+1));
            try
            caxis([minnestMat(b,1) maxestMat(b,1)]);
            catch
                %nothing
            end
            colorbar;
            count = count + 1;
        end
    end
    
    sgtitle(regexprep(Groups{G},'_',' ','emptymatch'))
end

%% Plot stim length and magnitude (bar graphs)
% USEFUL FOR LOOKING AT LOCO SPEED FOR LOCO BLOCKS
% 
% for r = 1:length(unique_results)
%     for b = 1:length(behaviorToPlot)
% 
%         stimValueMat = nan(1,length(Groups));
%         stimValueSEM = nan(1,length(Groups));
%         stimLengthMat = nan(1,length(Groups));
%         stimLengthSEM = nan(1,length(Groups));
%         
%         for G = 1:length(Groups) %One figure per group
% 
%             GroupRows = strcmp(ConcatData.BlockInfo.Group, Groups{G});
% 
%             isResult = strcmp(result(GroupRows),unique_results{r});
%             dataToPlot = ConcatData.(behaviorToPlot{b})(GroupRows,:);
%             dataToPlot = dataToPlot(isResult,:);
%             stimLength = ConcatData.Stim_Length(GroupRows);
%             stimLength = stimLength(isResult);
% 
%             %Remove traces with nans
%             nanrows = sum(isnan(dataToPlot),2) > 0;
%             dataToPlot(nanrows,:) = [];
%             stimLength(nanrows,:) = [];
% 
%             %Average
%             mean_DataToPlot = mean(dataToPlot,2,'omitnan');
%             
%             stimValueMat(G) = mean(mean_DataToPlot,'omitnan');
%             stimValueSEM(G) = std(mean_DataToPlot,'omitnan')/sqrt(length(mean_DataToPlot)-1);
%             stimLengthMat(G) = mean(stimLength,'omitnan');
%             stimLengthSEM(G) = std(stimLength,'omitnan')/sqrt(length(stimLength)-1);
%         end
%         
%         figure;
%         subplot(2,1,1); hold on
%         bar(1:length(Groups), stimValueMat)
%         errorbar(1:length(Groups), stimValueMat, stimValueSEM)
%         set(gca, 'XTick', 1:length(Groups))
%         set(gca, 'XTickLabel', Groups)
%         ylabel(['Average ' behaviorToPlot{b}])
% 
%         subplot(2,1,2); hold on
%         bar(1:length(Groups), stimLengthMat)
%         errorbar(1:length(Groups), stimLengthMat, stimLengthSEM)
%         set(gca, 'XTick', 1:length(Groups))
%         set(gca, 'XTickLabel', Groups)
%         ylabel(['Average stim length'])
%         
%         sgtitle(strcat(unique_results{r}, {' '}, behaviorToPlot{b}))
%     end
% end
          