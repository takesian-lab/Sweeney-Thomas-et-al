% Plot XCorr functions from noise/signal correlation data

% Argument(s):
% Correlations_Matrix_All - output from Network Analysis Pipeline
% Correlations_Summary - output from Network Analysis Pipeline

% Output
% Graphs of xcorr functions

%VERSION HISTORY
% v1 = archived 6/26/23 by MET
% v2 = for use with v2 simple_motor_correlation (added significant max/min determination)

%% User options

datapath = 'Z:\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\NetworkAnalysis Spontaneous Synch new analysis';
plotSignificantOnly = 0;
cctype = 1; %1 to use xcorr @ zero lag, 2 to use xcorr @max lag, 3 to use xcorr @min lag for determining +/- peaks and significance
plotMotor = 0; %0 if motor was not extracted during NetworkAnalysis run

%% Load data

cd(datapath)
disp('Loading files...')
load('Correlations_Summary_Spontaneous.mat')
load('Correlations_Table_Spontaneous.mat')
if plotMotor
    load('MotorCorrelations_Table.mat')
end

%% Plot xcorr functions for trace and motor corr. for all groups and loco types

plotTitle  = {'Trace Correlations',      'Motor Correlations',           'XY Correlations',              'Residuals',                    'Z Correlations'};
dataToPlot = {'XCorrMat',                'MotorCorrMat',                 'XYCorrMat',                    'ResidualsCorrMat',             'ZCorrMat'};
dataXlabel = {'XCorrLag',                'MotorCorrLag',                 'XYCorrLag',                    'ResidualsCorrLag',             'ZCorrLag'};
dataMatrix = {'Correlations_Matrix_All', 'MotorCorrelations_Matrix_All', 'MotorCorrelations_Matrix_All', 'MotorCorrelations_Matrix_All', 'MotorCorrelations_Matrix_All'};
dataTable  = {'XCorrTraces',             'MotorCorr',                    'XYCorr',                       'ResidualsCorr',                'ZCorr'};
[dataCorr{1:5}]   = deal('cc_zero');
[dataZTest{1:5}]  = deal('cc_zero_z');
[dataSign{1:5}]   = deal('cc_zero_sign');
[dataCorr2{1:5}]  = deal('cc_max');
[dataZTest2{1:5}] = deal('cc_max_z');
[dataSign2{1:5}]  = deal('cc_sign');
[dataCorr3{1:5}]  = deal('cc_min');
[dataZTest3{1:5}] = deal('cc_min_z');
[dataSign3{1:5}]  = deal('cc_sign');
UseMotorTable = [0, 1, 1, 1, 1];

for f = 1:length(dataToPlot) %One figure per data type
    figure; hold on

    % Sort Data According to Loco, Cell Type and Stim Type
    if UseMotorTable(f)
        
        %If no motor data to plot, skip figure
        if ~plotMotor
            continue;
        end
        
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

    % Sort by Locomotion
    if ~UseMotorTable(f)
        LocoList = eval(dataMatrix{f}).LocoType;
        Loco = unique(LocoList);
    else
        LocoList = repmat("All", size(eval(dataMatrix{f}),1), 1);
        Loco = unique(LocoList);
    end

    %For plot legend
    nMat = nan(3,length(Loco),length(Groups));

    for L = 1:length(Loco) %One subplot per loco type
        
        clear('ln');

       for G = 1:length(Groups) %One trace per group
           grouprows = strcmp(GroupList, Groups{G});

           traces = [];
           correlations = [];
           signs = [];

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
                   new_name = makeFieldNameFromBlock(Blocks(B));
                   if ~isfield(Correlations_Summary.(Mice{M}).(Stim).(new_name), dataToPlot{f})
                       continue;
                   end
                   block_traces = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataToPlot{f}){1,L};
                   lag = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataXlabel{f}){1,L};
                   blockTable = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataTable{f}){1,L};
                   if cctype == 1
                       block_correlations = blockTable.(dataCorr{f});
                       ZTest = blockTable.(dataZTest{f});
                       block_signs = double(block_correlations >= 0);
                       block_signs(block_signs == 0) = -1;
                   elseif cctype == 2
                       block_correlations = blockTable.(dataCorr2{f});
                       ZTest = blockTable.(dataZTest2{f});
                       block_signs = blockTable.(dataSign2{f});
                       ZTest(block_signs == -1) = 0; %don't count min xcorrs as significant
                   elseif cctype == 3
                       block_correlations = blockTable.(dataCorr3{f});
                       ZTest = blockTable.(dataZTest3{f});
                       block_signs = blockTable.(dataSign3{f});
                       ZTest(block_signs == 1) = 0; %don't count min xcorrs as significant
                   end
                   
                   if plotSignificantOnly == 1
                       block_correlations = block_correlations(ZTest == 1,:);
                       block_signs = block_signs(ZTest == 1,:);
                       block_traces = block_traces(ZTest == 1,:);
                   end

                   correlations = [correlations; block_correlations];
                   signs = [signs; block_signs];
                   traces = [traces; block_traces];
               end
           end

           %Plot figure
           ylabels = {'Both', 'Positive', 'Negative'};
           for ff = 1:3
               if ff == 1
                   ind = true(size(correlations));
               elseif ff == 2
                   ind = signs > 0;
               elseif ff == 3
                   ind = signs < 0;
               end

               %Record number of pairs
               nMat(ff,L,G) = sum(ind);

               if isempty(ind); continue; end
               
               subplot(3,length(Loco),L+((ff-1)*(length(Loco)))); hold on
               x = lag(1,:);
               y = mean(traces(ind,:),1, 'omitnan');
               y_sem = std(traces(ind,:),[],1,'omitnan')./sqrt(size(traces,1)-1);
               shadedErrorBar(x, y, y_sem);
               ln(ff,G) = plot(x,y,'LineWidth',2);
               xlabel('Lag (seconds)')
               ylabel(ylabels{ff})

               title(Loco{L})
           end
       end 
    end

    %Add legends
    for L = 1:length(Loco)
        for ff = 1:3
            groupLegend = Groups;
            for G = 1:length(Groups)
                groupLegend{G} = [groupLegend{G} ' ' num2str(nMat(ff,L,G)) ' pairs'];
            end
            subplot(3,length(Loco),L+((ff-1)*(length(Loco)))); hold on
            subtitle(groupLegend);
            if exist('ln', 'var'); legend(ln(ff,:), groupLegend, 'Location', 'Northwest'); end
        end
    end

    sgtitle(plotTitle{f})
end

%% Plot heat maps for motor corr.

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
L = 1; %TODO: add back in Loco if we add trace correlations to this section

for f = 1:length(dataToPlot) %One figure per data type
    figure; hold on

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

    %For plot legend
    nMat = nan(3,length(Loco),length(Groups));

   for G = 1:length(Groups) %One trace per group
       grouprows = strcmp(GroupList, Groups{G});

       traces = [];
       correlations = [];

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

               %alter block name to be compatible as field name
               new_name = makeFieldNameFromBlock(block_name);
               if ~isfield(Correlations_Summary.(Mice{M}).(Stim).(new_name), dataToPlot{f})
                   continue;
               end
               block_traces = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataToPlot{f}){1,L};
               lag = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataXlabel{f}){1,L};
               blockTable = Correlations_Summary.(Mice{M}).(Stim).(new_name).(dataTable{f}){1,L};
               if cctype == 1
                   block_correlations = blockTable.(dataCorr{f});
                   ZTest = blockTable.(dataZTest{f});
               elseif cctype == 2
                   block_correlations = blockTable.(dataCorr2{f});
                   ZTest = blockTable.(dataZTest2{f});
               end

               if plotSignificantOnly == 1
                   block_correlations = block_correlations(ZTest == 1,:);
                   block_traces = block_traces(ZTest == 1,:);
               end

               correlations = [correlations; block_correlations];
               traces = [traces; block_traces];
           end
       end

       %Plot figure
       ylabels = {'Both', 'Positive', 'Negative'};
       for ff = 1:3
           if ff == 1
               ind = true(size(correlations));
           elseif ff == 2
               ind = correlations > 0;
           elseif ff == 3
               ind = correlations < 0;
           end

           %Record number of pairs
           nMat(ff,L,G) = sum(ind);

           if ~any(ind); continue; end
           
           subplot(3,length(Groups),G+((ff-1)*(length(Groups)))); hold on

           %Sort traces by max value
           y = traces(ind,:);
           [maxvals, ~] = max(y,[],2);
           [~, sortind] = sort(maxvals, 'ascend');
           imagesc(y(sortind,:));
           xlim([0.5 size(y,2)])
           ylim([0.5 size(y,1)])
           xlabel('Lag (frames)')
           ylabel(ylabels{ff})
           caxis([-0.5 0.5])
           
           title([Groups{G} ' ncells = ' num2str(sum(ind))])
       end
    end

    %Add legends
    sgtitle(plotTitle{f})
end