function plot_simple_ExtractedData_Controls(ExtractedData, UserOps)
% Plot results of ExtractedData by control variables like age, sex, field, depth, and gcamp type
% Figure1 = Percent sound-responsive cells by group and control variable
% Figure2 = Scatter plot of peak and stim response values by group and control variable

%% Turn plots on or off

plotfig1 = 1; %Percent activated/suppressed/prolonged cells
plotfig2 = 0; %Scattter plots of peak and stim properties

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
    UserOps.plotCorrelations	= 0; %0 = don't plot dependent variable correlations, 1 = plot
    UserOps.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
    UserOps.FibPhotChannel      = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
    UserOps.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress
end

disp('Plotting control data')

%% Controls
controlColumns = {'GCaMP', 'Sex', 'Age', 'ImagingDepth', 'Field1', 'Field2'};
variableType = {'Nominal', 'Nominal', 'Continuous', 'Continuous', 'Nominal', 'Nominal'};
    
%% Get data from ExtractedData *FIRST TIME*
%This time we want to include all cells, even non-responsive ones, to make pie charts

Ops = UserOps; %Copy and overwrite some things
Ops.ResponseType        = ''; %No filtering
Ops.RF_Type             = ''; %No filtering
Ops.IsResponsive        = 2; %Both responsive and non-responsive

SubsetData = simple_subset_ExtractedData(ExtractedData, Ops);

Ops = SubsetData.SubsetOps;
StimInfo = SubsetData.StimInfo;
Summary = SubsetData.Summary;
GroupList = Summary.Group;
Groups = unique(GroupList);

%Add "All" column
if length(Groups) > 1
    Groups{end+1} = 'All';
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

%% FIGURE 1: Percent sound-responsive cells by group and control variable

if plotfig1
    
ResponseTypesToPlot = {'activated', 'prolonged', 'suppressed', 'none'};

%Plot one figure per control variable
for c = 1:length(controlColumns)
    
    %Skip control columns not present in data
    if ~ismember(controlColumns{c}, Summary.Properties.VariableNames)
        continue;
    end
        
    %Skip GCaMP control if already sorting by GCaMP
    if strcmp(controlColumns{c}, 'GCaMP') && Ops.sortbyGCAMP == 1
        continue;
    end
    
    ControlColumn = Summary.(controlColumns{c});
    ResponseType = Summary.ResponseType;
    MouseColumn = Summary.MouseID;
    TempGroupList = GroupList;
    
    %Remove nan or missing data
    remove_ind = ismissing(ControlColumn);
    ControlColumn(remove_ind) = [];
    ResponseType(remove_ind) = [];
    MouseColumn(remove_ind) = [];
    TempGroupList(remove_ind) = [];
    
    if strcmp(variableType{c}, 'Continuous') && isnumeric(ControlColumn) %Using the isnumeric condition allows these variables to be treated as nominal if they aren't numeric

        %Bin continuous variables into 4 bins
        minc = min(ControlColumn);
        maxc = max(ControlColumn);
        rangec = maxc - minc;
        
        %Skip if all the cells in the dataset have the same control (nothing to compare!)
        if rangec == 0
            continue;
        end
            
        binEdges = round(minc:rangec/4:maxc);
        
        edgeLabels = strings(1,4);
        for i = 1:length(edgeLabels)
            edgeLabels{i} = [num2str(binEdges(i)) '-' num2str(binEdges(i+1))];
        end

        mat1 = ControlColumn >= binEdges(1:4);
        mat2 = ControlColumn <= binEdges(2:end);
        mat3 = (mat1 + mat2) == 2;
        
        DiscreteControlColumn = strings(size(ControlColumn));
        for i = 1:length(ControlColumn)
            DiscreteControlColumn(i) = edgeLabels{mat3(i,:)==1};
        end
        
        ControlColumn = DiscreteControlColumn;
        Controls = edgeLabels; %Prevents reordering with unique fxn
    else
        Controls = unique(ControlColumn);
    end
 
    %Skip if all the cells in the dataset have the same control (nothing to compare!)
    if length(Controls) == 1
        continue;
    else
        %Add "All" columns
        Controls{end+1} = 'All';
    end
    
    figure('units','normalized','outerposition',[0 0 1 1]); hold on

    for g = 1:length(Groups)
        
        if g == length(Groups)
            currentData = ResponseType;
            currentControl = ControlColumn;
            currentMice = MouseColumn;
        else
            groupInd = strcmp(Groups{g}, TempGroupList);
            currentData = ResponseType(groupInd);
            currentControl = ControlColumn(groupInd);
            currentMice = MouseColumn(groupInd);
        end
        
        XTickLabel = strings(1,length(Controls));
        nMiceRepresented = nan(1,length(Controls));
        nDataToPlot = nan(length(Controls), length(ResponseTypesToPlot));
        percentDataToPlot = nan(length(Controls), length(ResponseTypesToPlot));
        for i = 1:length(Controls)
            if i < length(Controls)
                ind = strcmp(currentControl, Controls{i});
            else
                %All group includes all data
                ind = true(size(currentControl));
            end
            
            %Store number of mice included in each control
            nMiceRepresented(i) = length(unique(currentMice(ind))); 
            XTickLabel(i) = [Controls{i} ' (' num2str(nMiceRepresented(i)) ')'];
            
            tabData = tabulate(currentData(ind));
            
            for ii = 1:length(ResponseTypesToPlot)
                ind2 = find(strcmp(tabData(:,1), ResponseTypesToPlot{ii}));
                if ~isempty(ind2)
                    nDataToPlot(i,ii) = tabData{ind2,2};
                    percentDataToPlot(i,ii) = tabData{ind2,3};
                end
            end
        end

        %First row = percent out of 100 including None
        subplot(4,length(Groups), g); hold on
        bar(percentDataToPlot, 'stacked');
        ylabel('Percent')
        legend(ResponseTypesToPlot)
        title(GroupsLabel{g})
        set(gca, 'XTick', 1:1:length(Controls))
        set(gca, 'XTickLabel', XTickLabel, 'fontsize', 8)
        ylim([0 100])

        %Second row = count including None
        subplot(4,length(Groups), g + length(Groups)); hold on
        bar(nDataToPlot, 'stacked');
        ylabel('Count')
        title(GroupsLabel{g})
        set(gca, 'XTick', 1:1:length(Controls))
        set(gca, 'XTickLabel', XTickLabel, 'fontsize', 8)
        
        %Third row = percent out of 100 without None
        ResponseTypesToPlotWithoutNone = ResponseTypesToPlot;
        noneInd = strcmp('none', ResponseTypesToPlot);
        nDataToPlot(:,noneInd) = [];
        percentDataToPlot = 100*(nDataToPlot./sum(nDataToPlot,2,'omitnan')); %Recompute percent with none removed
        ResponseTypesToPlotWithoutNone(:,noneInd) = [];
        
        subplot(4,length(Groups), g + length(Groups)*2); hold on
        bar(percentDataToPlot, 'stacked');
        ylabel('Percent')
        title(GroupsLabel{g})
        set(gca, 'XTick', 1:1:length(Controls))
        set(gca, 'XTickLabel', XTickLabel, 'fontsize', 8)
        ylim([0 100])

        %Fourth row = count without None
        subplot(4,length(Groups), g + length(Groups)*3); hold on
        bar(nDataToPlot, 'stacked');
        ylabel('Count')
        title(GroupsLabel{g})
        set(gca, 'XTick', 1:1:length(Controls))
        set(gca, 'XTickLabel', XTickLabel, 'fontsize', 8)
    end
    sgtitle([StimInfo.StimType ' ' controlColumns{c}])
    drawnow
end

end

%% Get data from ExtractedData *SECOND TIME*
%This time we want to include just sound responsive cells to plot differences in peak and stim response variables

SubsetData = simple_subset_ExtractedData(ExtractedData, UserOps);

Ops = SubsetData.SubsetOps;
StimInfo = SubsetData.StimInfo;
Summary = SubsetData.Summary;
GroupList = Summary.Group;
Groups = unique(GroupList);
peaktype = SubsetData.peaktype;

%Add "All" column
if length(Groups) > 1
    Groups{end+1} = 'All';
end

if isempty(Summary)
    disp('No responsive cells found to plot for controls')
    return
end

%% FIGURE 2: Scatter plot of peak and stim response values by group and control variable
     
if plotfig2
    
%Peak variables:
variablesToPlot = {'_AUC', 'R2'}; %{'', '_3fr', '_win', '_avg', '_Onset', '_Latency', '_Width', '_AUC', 'R2'};

%Concatenate Peak or Trough to the beginning of each variables, except R2
for v = 1:length(variablesToPlot)
    if strcmp(variablesToPlot{v}, 'R2')
        continue;
    end
    variablesToPlot{v} = strcat(peaktype, variablesToPlot{v});
end

%Add additional variables specific to stim protocol
switch StimInfo.StimProtocol
    case 2 %RF
        variablesToPlot = [variablesToPlot, {'Sparseness', 'BWBF_I', 'BF_I', 'ISI'}]; %{'BF', 'BF_I', 'BWBF_I', 'CF', 'CF_I', 'BW20', 'dPrime', 'Sparseness','BestInt', 'BWInt', 'ISI'};
    case 3 %FM
        variablesToPlot = [variablesToPlot, {'Sparseness', 'BestSpeed_fit', 'BestAbsSpeed_fit', 'MI'}]; %{'BestSpeed', 'BestSpeed_fit', 'Sparseness', 'BestAbsSpeed', 'BestAbsSpeed_fit', 'SSI', 'MI', 'SMI', 'Slope'};
    case 5 %SAM
        variablesToPlot = [variablesToPlot, {'Sparseness', 'BestRate_fit', 'RSI', 'BestDepth_fit', 'DSI'}]; %{'BestRate', 'BestRate_fit', 'BestRate_fit_halfwidth', 'RSI', 'BestDepth', 'BestDepth_fit', 'BestDepth_fit_halfwidth', 'DSI', 'dprime', 'Sparseness'};
    case 6 %SAMfreq
        variablesToPlot = [variablesToPlot, {'Sparseness', 'BF_fit_halfwidth', 'BestDepth_fit', 'DSI'}]; %{'BF', 'BF_fit', 'BF_fit_halfwidth', 'BestDepth', 'BestDepth_fit', 'BestDepth_fit_halfwidth', 'DSI', 'dprime', 'Sparseness'};
end

%Remove variable columns that aren't in table
remove_ind = [];
for v = 1:length(variablesToPlot)
    if ~ismember(variablesToPlot{v}, Summary.Properties.VariableNames)
        remove_ind = [remove_ind, v];
    end
end
variablesToPlot(remove_ind) = [];

if isempty(variablesToPlot)
    disp('No peak or stim variables found to plot')
    return;
end

%Plot one figure per control variable
for c = 1:length(controlColumns)
    
    %Skip control columns not present in data
    if ~ismember(controlColumns{c}, Summary.Properties.VariableNames)
        continue;
    end
        
    %Skip GCaMP control if already sorting by GCaMP
    if strcmp(controlColumns{c}, 'GCaMP') && Ops.sortbyGCAMP == 1
        continue;
    end
    
    h = figure('units','normalized','outerposition',[0 0 1 1]); hold on
    tiledlayout(length(variablesToPlot),length(Groups))
        
    for v = 1:length(variablesToPlot)
        ControlColumn = Summary.(controlColumns{c});
        ResponseData = Summary.(variablesToPlot{v});
        MouseColumn = Summary.MouseID;
        TempGroupList = GroupList;

        %Remove nan or missing data
        remove_ind = ismissing(ControlColumn);
        ControlColumn(remove_ind) = [];
        ResponseData(remove_ind) = [];
        MouseColumn(remove_ind) = [];
        TempGroupList(remove_ind) = [];

        %Get min and max of ResponseData for ylim
        miny = min(ResponseData);
        maxy = max(ResponseData);
        
        if strcmp(variableType{c}, 'Continuous')  && isnumeric(ControlColumn)

            %Bin continuous variables into 4 bins
            minc = min(ControlColumn);
            maxc = max(ControlColumn);
            rangec = maxc - minc;
            
            %Skip if all the cells in the dataset have the same control (nothing to compare!)
            if rangec == 0
                if isvalid(h)
                    close(h)
                end
                continue;
            end
            
            binEdges = round(minc:rangec/4:maxc);

            edgeLabels = strings(1,4);
            for i = 1:length(edgeLabels)
                edgeLabels{i} = [num2str(binEdges(i)) '-' num2str(binEdges(i+1))];
            end

            mat1 = ControlColumn >= binEdges(1:4);
            mat2 = ControlColumn <= binEdges(2:end);
            mat3 = (mat1 + mat2) == 2;

            DiscreteControlColumn = strings(size(ControlColumn));
            for i = 1:length(ControlColumn)
                DiscreteControlColumn(i) = edgeLabels{mat3(i,:)==1};
            end

            ControlColumn = DiscreteControlColumn;
            Controls = edgeLabels; %Prevents reordering with unique fxn
        else
            Controls = unique(ControlColumn);
        end

        %Skip if all the cells in the dataset have the same control (nothing to compare!)
        if length(Controls) == 1
            if isvalid(h)
                close(h)
            end
            continue;
        end
        
        for g = 1:length(Groups)

            if g == length(Groups)
                currentData = ResponseData;
                currentControl = ControlColumn;
                currentMice = MouseColumn;
            else
                groupInd = strcmp(Groups{g}, TempGroupList);
                currentData = ResponseData(groupInd);
                currentControl = ControlColumn(groupInd);
                currentMice = MouseColumn(groupInd);
            end

            XTickLabel = strings(1,length(Controls));
            nMiceRepresented = nan(1,length(Controls));
            DataToPlot = cell(1,length(Controls));
            for i = 1:length(Controls)
                ind = strcmp(currentControl, Controls{i});
                
                %Store number of mice included in each control
                nMiceRepresented(i) = length(unique(currentMice(ind))); 
                XTickLabel(i) = [Controls{i} ' (' num2str(nMiceRepresented(i)) ')'];

                DataToPlot{i} = currentData(ind);
            end

            %First row = percent out of 100 including None
            nexttile; hold on
            if length(Controls) < 3
                colours = flipud(cbrewer('qual','Set1',3)); %Formula will print warning if ncol is less than 3
            else
                colours = flipud(cbrewer('qual','Set1',length(Controls)));
            end
            for d = 1:size(DataToPlot,2)
                meanY = mean(DataToPlot{d},'omitnan');
                seY = std(DataToPlot{d},'omitnan')/sqrt(length(DataToPlot{d})-1);
                scatter(zeros(length(DataToPlot{d}))+d, DataToPlot{d}, 25, colours(d,:), 'filled', 'jitter', 'on', 'jitterAmount', 0.05);
                line([d-0.2 d+0.2], [meanY meanY], 'Color', 'k', 'Linewidth', 1.5)
                errorbar(d, meanY, seY, 'k', 'LineWidth', 1.5)
                %line([d d], [meanY-stdY, meanY+stdY], 'Color', 'k', 'Linewidth', 1.5)
            end
            ylabel(regexprep(variablesToPlot{v},'_',' ','emptymatch'))
            xlim([0 length(Controls)+0.5])
            try
                ylim([miny maxy])
            catch
                %don't set ylim
            end
            title(GroupsLabel{g})
            set(gca, 'XTick', 1:1:length(Controls))
            set(gca, 'XTickLabel', XTickLabel, 'fontsize', 8)
        end
        
        sgtitle([StimInfo.StimType ' ' controlColumns{c}])
        drawnow
    end
end

end
