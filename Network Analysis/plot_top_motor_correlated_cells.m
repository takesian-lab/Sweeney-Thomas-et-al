%plot_top_motor_correlated_cells
% Figure type 1: plot top N positively and negatively correlated cells per group
% Figure type 2: plot loco trace, +/- xcorr, and histographs for all blocks in order of % time spent running

%% Set options and paths

cctype = 1; %1 to use xcorr @ zero lag, 2 to use xcorr @max lag for determining +/- peaks and significance
nCellsToPlot = 20; %FOR CELL PLOTS
plotSignificantOnly = 0; %FOR BLOCK PLOTS

%SET PATHS
ExtractedData_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\Loco correlation analysis';
ExtractedData_file = 'ExtractedData_Spontaneous_20230520-183014_LocoAnalysis';
MotorCorrelations_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\NetworkAnalysis Loco Threshold with XYZR min vals';
MotorCorrelations_file = 'MotorCorrelations_Table_Spontaneous.mat';
CorrelationsSummary_file = 'Correlations_Summary_Spontaneous.mat';

%% Load files
cd(MotorCorrelations_path)
load(MotorCorrelations_file)
load(CorrelationsSummary_file)

cd(ExtractedData_path)
load(ExtractedData_file)

%% FIGURE 1: Find top + and - motor correlated cells per group

GroupList = MotorCorrelations_Matrix_All.MotorCellType;
Groups = unique(GroupList);

for G = 1:length(Groups)
    
    GroupData = MotorCorrelations_Matrix_All(strcmp(GroupList,Groups{G}),:);
    
    if cctype == 1
        correlations = GroupData.MotorCorrelations;
    elseif cctype == 2
        corrIsNegative = GroupData.MotorMaxMinSign == -1;
        correlations = GroupData.MotorMaxCorrelations;
        correlations(corrIsNegative) = GroupData.MotorMinCorrelations;
    end
    [pos_sorted_correlations, pos_sort_ind] = sort(correlations,'descend');
    [neg_sorted_correlations, neg_sort_ind] = sort(correlations,'ascend');
    
    topPosCells = GroupData(pos_sort_ind(1:nCellsToPlot),:);
    topNegCells = GroupData(neg_sort_ind(1:nCellsToPlot),:);
    
    %Plot figures    
    for f = 1:2 %pos and neg
        if f == 1
            currentCells = topPosCells;
            traceColor = 'm';
        elseif f == 2
            currentCells = topNegCells;
            traceColor = 'c';
        end

        for c = 1:height(currentCells)
            mouse = currentCells.MotorMouseName(c);
            block = currentCells.MotorBlockName(c);
            stim = currentCells.MotorStimName(c);
            cellNumber = currentCells.Cell_s2P_list(c);
            
            %Get cellOrder from CellList
            block_CellList = ExtractedData.CellList(strcmp(ExtractedData.CellList.Block, block),:);
            cellOrder = block_CellList.CellOrder(block_CellList.CellNumber == cellNumber);

            %Get full traces from TrialData
            b = find(strcmp(({ExtractedData.TrialData.Block})', block));
            cell_timestamp = ExtractedData.TrialData(b).Timestamp;
            loco = ExtractedData.TrialData(b).Full_Loco;
            loco(loco < 0.8) = 0;
            cell_F = ExtractedData.TrialData(b).F7(cellOrder,:);

            if isfield(ExtractedData.TrialData, 'Full_xoff')
                XY = sqrt(ExtractedData.TrialData(b).Full_xoff.^2 + ExtractedData.TrialData(b).Full_yoff.^2);
            else
                XY = [];
            end
            
            %Compute locomotor residuals
            %Perform linear regression between loco and XY trace and compute residuals
            XY_zscored = zscore(XY,'',2);
            loco_zscored = zscore(loco,'',2);
            p = polyfit(XY_zscored, loco_zscored, 1);
            LocoFit = polyval(p, XY_zscored);
            LocoResiduals = loco_zscored - LocoFit;
            LocoResiduals_zscored = zscore(LocoResiduals,'',2);
            
            %Get xcorr function from Correlations_Summary
            new_name = makeFieldNameFromBlock(block);
            Correlations = Correlations_Summary.(mouse).(stim).(new_name);
            cellRow = find(Correlations.CellNumbers.S2P{1,1} == cellNumber);
            xmat = Correlations.MotorCorrMat{1,1}(cellRow,:);
            lag = Correlations.MotorCorrLag{1,1}(1,:);
            
            %Plot
            figure;
            subplot(4,1,1); hold on
            %dff = (cell_F - mean(cell_F,'omitnan'))./(mean(cell_F,'omitnan'));
            cell_F_Z = zscore(cell_F,'',2);
            plot(cell_timestamp, smooth(cell_F_Z), traceColor, 'Linewidth', 2);
            xlim([cell_timestamp(1) cell_timestamp(end)])
            title(['Motor correlation = ' num2str(round(currentCells.MotorCorrelations(c),3)) ' ZTest = ' num2str(currentCells.MotorZTest(c))])
            ylabel('Z')
            
            subplot(4,1,2); hold on
            plot(cell_timestamp, smooth(loco), 'k', 'Linewidth', 2);
            plot(cell_timestamp, smooth(LocoResiduals_zscored), 'r', 'Linewidth', 2);
            xlim([cell_timestamp(1) cell_timestamp(end)])
            legend({'Raw Loco', 'Zscored Residuals'})
            title('Loco')
            ylabel('Speed (cm/s)')
            
            subplot(4,1,3); hold on
            if ~isempty(XY)
                plot(cell_timestamp, XY, 'k');
            end
            xlim([cell_timestamp(1) cell_timestamp(end)])
            title('XY shift')
            xlabel('Time (s)')
            ylabel('Pixels')
            
            subplot(4,1,4); hold on
            plot(lag, xmat)
            xlabel('Lag (s)')
            ylabel('xcorr')
            %ylim([-0.5 0.5])
            hline(0)

            sgtitle([regexprep(block,'_',' ','emptymatch') ' Cell Number ' num2str(cellNumber)])
        end
    end    
end

%% FIGURE 2: Plot block previews for motor corr.

plotTitle = {'Motor Correlations', 'XY Correlations', 'Residuals', 'Z Correlations'};
dataToPlot = {'MotorCorrMat', 'XYCorrMat', 'ResidualsCorrMat', 'ZCorrMat'}; %Add NoiseCorr later
dataXlabel = {'MotorCorrLag', 'XYCorrLag', 'ResidualsCorrLag', 'ZCorrLag'};
dataMatrix = {'MotorCorrelations_Matrix_All', 'MotorCorrelations_Matrix_All', 'MotorCorrelations_Matrix_All', 'MotorCorrelations_Matrix_All'};
dataTable = {'MotorCorr', 'XYCorr', 'ResidualsCorr', 'ZCorr'};
dataCorr = {'cc_zero', 'cc_zero', 'cc_zero', 'cc_zero'};
dataZTest = {'cc_zero_z', 'cc_zero_z', 'cc_zero_z', 'cc_zero_z'};
dataCorr2 = {'cc_max', 'cc_max', 'cc_max', 'cc_max'};
dataZTest2 = {'cc_max_z', 'cc_max_z', 'cc_max_z', 'cc_max_z'};
UseMotorTable = [1, 1, 1, 1];

for f = 1:length(dataToPlot) %One figure per data type

    % Sort Data According to Loco, Cell Type and Stim Type
    if UseMotorTable(f)
        CellColumn = 'MotorCellType';
        StimColumn = 'MotorStimName';
        MouseColumn = 'MotorMouseName';
        BlockColumn = 'MotorBlockName';
    else
        CellColumn = 'CellType';
        StimColumn = 'StimName';
        MouseColumn = 'MouseName';
        BlockColumn = 'BlockName';
    end

    % Sort by Cell Types
    GroupList = eval(dataMatrix{f}).(CellColumn);
    Groups = unique(GroupList);

    % Sort by Stim Type
    StimList = eval(dataMatrix{f}).(StimColumn);
    Stim = unique(StimList);
    if length(Stim) > 1
        warning('This code is currently not written for more than one stim type. All stim will be plotted together')
    end

   for G = 1:length(Groups) %One trace per group
              
       grouprows = strcmp(GroupList, Groups{G});

       blockcount = 1;
       GroupData = struct;
       GroupData.blocknames = {};
       GroupData.time_spent_running = [];
       GroupData.timestamps = {};
       GroupData.loco_traces = {};
       GroupData.xcorr_traces = {};
       GroupData.correlations = {};
       GroupData.all_correlations = {};

       %Find mice belonging to this group
       MouseList = eval(dataMatrix{f}).(MouseColumn);
       Mice = unique(MouseList(grouprows));

       for M = 1:length(Mice)
           mouserows = strcmp(MouseList, Mice{M});

           %Find blocks for this mouse
           BlockList = eval(dataMatrix{f}).(BlockColumn)(mouserows,:);
           Blocks = unique(BlockList);

           %Go into each Correlations_Summary.Mouse.Stim.Block to get the traces
           for B = 1:length(Blocks)

               %GET LOCO TRACES FROM TRIAL DATA
               b = find(strcmp(({ExtractedData.TrialData.Block})', Blocks(B)));
               cell_timestamp = ExtractedData.TrialData(b).Timestamp;
               loco = ExtractedData.TrialData(b).Full_Loco;
               
               %find times within trace with running for each analyzed (responsive) neuron
               loco(loco < 0.8) = 0; %Correct for noise floor
               block_time_spent_running = sum(loco > 0)/length(loco);

               %GET XCORR TRACES FROM CORRELATIONS_SUMMARY
               new_name = makeFieldNameFromBlock(Blocks(B));
               if ~isfield(Correlations_Summary.(Mice{M}).(Stim).(new_name), dataToPlot{f})
                   continue;
               end
               block_traces = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataToPlot{f}){1,1};
               lag = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataXlabel{f}){1,1};
               blockTable = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataTable{f}){1,1};
               if cctype == 1
                   block_correlations = blockTable.(dataCorr{f});
                   ZTest = blockTable.(dataZTest{f});
               elseif cctype == 2
                   block_correlations = blockTable.(dataCorr2{f});
                   ZTest = blockTable.(dataZTest2{f});
               end

               all_block_correlations = block_correlations;
               
               if plotSignificantOnly == 1
                   block_correlations = block_correlations(ZTest == 1,:);
                   block_traces = block_traces(ZTest == 1,:);
               end

               %Store data
               GroupData(blockcount).blocknames = new_name;
               GroupData(blockcount).time_spent_running = block_time_spent_running;
               GroupData(blockcount).all_correlations = all_block_correlations; %For histogram
               GroupData(blockcount).correlations = block_correlations;
               GroupData(blockcount).xcorr_traces = block_traces;
               GroupData(blockcount).timestamps = cell_timestamp;
               GroupData(blockcount).loco_traces = loco;
               blockcount = blockcount+1;
           end
       end
       
       %Sort GroupData by loco speed
       [~, sortind] = sort([GroupData.time_spent_running], 'descend');
       GroupData = GroupData(:,sortind);
       
       %Plot figure
       figure; tiledlayout(size(GroupData,2), 4)
       
       for B = 1:size(GroupData,2)
           
           %loco trace
           nexttile; hold on
           plot(GroupData(B).timestamps, GroupData(B).loco_traces, 'k');
           xlim([GroupData(B).timestamps(1) GroupData(B).timestamps(end)])
           ylim([0 30])
           maxspeed = round(max(GroupData(B).loco_traces),2);
           title(['Time running = ' num2str(round(GroupData(B).time_spent_running*100)) '% Max speed = ' num2str(maxspeed)])
       
           %Positive correlations
           nexttile; hold on
           traces = GroupData(B).xcorr_traces;
           correlations = GroupData(B).correlations;
           ind = correlations > 0;
           
           if any(ind)
               y = mean(traces(ind,:),1);
               y_sem = std(traces(ind,:),[],1,'omitnan')./sqrt(size(traces,1)-1);
               shadedErrorBar(1:length(y), y, y_sem)
               vline(length(y)/2)
           end
           set(gca,'xtick',[])
           xlim([1 length(y)])
           title([num2str(sum(ind)) ' positive'])
           
           %Negative correlations
           nexttile; hold on
           ind = correlations < 0;
           
           if any(ind)
               y = mean(traces(ind,:),1);
               y_sem = std(traces(ind,:),[],1,'omitnan')./sqrt(size(traces,1)-1);
               shadedErrorBar(1:length(y), y, y_sem)
               vline(length(y)/2)
           end
           set(gca,'xtick',[])
           xlim([1 length(y)])
           title([num2str(sum(ind)) ' negative'])
           
           %Histogram
           nexttile; hold on
           Ops.histogram_type = 'probability';
           Ops.bins = [-0.4:0.01:0.4];
           histogram(GroupData(B).all_correlations,Ops.bins, 'Normalization',Ops.histogram_type,'EdgeColor','k','FaceColor', 'k'); 
           histogram(correlations,Ops.bins, 'Normalization',Ops.histogram_type,'EdgeColor','g','FaceColor', 'g'); hold on;
           vline(0,'k')
           title(regexprep(GroupData(B).blocknames,'_',' ','emptymatch'))
       end
       sgtitle([Groups{G} ' - ' plotTitle{f}])
   end
end