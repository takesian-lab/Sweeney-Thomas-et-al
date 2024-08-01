function plot_simple_ExtractedData(ExtractedData, UserOps)
%Plot summary figures for ExtractedData created by simple_extract_data

%V1 archived 2/7/2023
%V2 = current version, uses simple_subset_ExtractedData

%% User options

%Activity types to analyze
%'none' can be removed/added from this list to include in plots or not
Activity = {'activated', 'prolonged', 'suppressed', 'none'}; %'none'

%I RECOMMEND NOT CHANGING THESE AND USING USEROPS INSTEAD
if nargin < 2
    UserOps = struct;
    UserOps.Loco                 = 'All'; %All, Running, or NotRunning
    UserOps.ResponseType         = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
    UserOps.RF_Type              = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
    UserOps.IsResponsive         = 1; %1 for responsive (and reliable), 0 for not, 2 for both
    UserOps.FOV                  = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1)
    UserOps.sortbyGCAMP          = 0; %0 for groups, 1 for gcamp, 2 to combine
    UserOps.sortbyCondition      = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    UserOps.sortbyRedCell        = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    UserOps.plotCorrelations	 = 0; %0 = don't plot dependent variable correlations, 1 = plot
    UserOps.smooth_rasters       = 1; %0 = don't smooth, 1 = smooth
    UserOps.plotShadedErrorBars  = 0; %0 = don't include error bars, 1 = do
    UserOps.FibPhotChannel       = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
    UserOps.SuppressOutput       = 0; %0 to print Ops to command line, 1 to suppress
end

disp('Plotting activity rasters and pie charts')

cb = 2; % 1:use color brewer; 2: use TAK lab colors to make predefined color schemes for graphs
if cb >0
    % choose colors here:
    colorsForGraphs = {'Purples', 'Reds','Greens'};
    % if cb = 1 some color options are:
    %{'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
    % 'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Spectral'};
    % but you can use cbrewer.mat for more options
    % if cb = 2, you are using TAK lab color schemes and
    % your options are {'Blues','Greens','Purples', 'Reds'};
end

%% Get set up for plotting

%Loco types to analyze
Loco = fieldnames(ExtractedData.CellData);

%Account for data with no Group
if any(ismissing(ExtractedData.CellList.Group))
    ExtractedData.CellList.Group(ismissing(ExtractedData.CellList.Group)) = "none";
end

%Do this to fill out user ops with anything that might be missing
Temp = simple_subset_ExtractedData(ExtractedData, UserOps);
UserOps = Temp.SubsetOps;

%Groups to analyze
if UserOps.sortbyGCAMP
    Groups = unique(ExtractedData.CellList.GCaMP);
else
    Groups = unique(ExtractedData.CellList.Group);
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end
            
%X labels for figures
xlabel_increment = 0.5; %seconds
fs = ExtractedData.StimInfo.fs;
nFrames = ExtractedData.StimInfo.nFrames;
increment_in_frames = xlabel_increment*fs;
xtick = 0:increment_in_frames:nFrames;
xticklabel = 0:xlabel_increment:((1/fs)*nFrames);
soundtime = ExtractedData.StimInfo.baseline*fs;

%% Sort data to plot

[groupRasters, groupAverages, groupActivity, groupSTD] = deal({});

for L = 1:length(Loco)
    for g = 1:length(Groups)
        average_raster = nan(length(Activity),nFrames);
        raster_SEM = nan(length(Activity),nFrames);
        resorted_raster = [];
        activity_type = [];
            
        for a = 1:length(Activity)

            %*NEW* Use function to subset data by Ops
            Ops = UserOps; %Copy and overwrite some of the values in UserOps
            Ops.Loco = Loco{L};
            Ops.ResponseType = Activity{a};
            Ops.RF_Type = ''; %Include everything
            Ops.IsResponsive = 2; %Include everything
            Ops.SuppressOutput = 1; %Prevent ops from printing to the command line
            
            %USE THIS DATA WHICH HAS ALREADY BEEN SUBSET
            SubsetData = simple_subset_ExtractedData(ExtractedData, Ops);
            Ops = SubsetData.SubsetOps; %This will contain any other variables that were missing
            FinalTraces = SubsetData.FinalTraces;
            Summary = SubsetData.Summary;
    
            % -------------- get data for rasters -----
            %find rows in current group and activity type and sort by amplitude and/or latency
            currentRows = strcmp(Groups{g}, Summary.Group);
            current_raster = FinalTraces(currentRows,:);
            
            if isempty(current_raster)
                continue;
            end
            
            %Store smoothed average for plots
            if Ops.smooth_rasters
                average_raster(a,:) = smooth(mean(current_raster,1,'omitnan'));
            else
                average_raster(a,:) = mean(current_raster,1,'omitnan');
            end

            %Store SEM for error bars
            N = size(current_raster,1) - 1; if N == 0; N = 1; end
            raster_SEM(a,:) = smooth(std(current_raster,0,1,'omitnan')/sqrt(N));

            %Use combined peak data to sort by peaks for activated/prolonged and troughs for suppressed
            if ismember(Activity{a}, {'activated', 'prolonged'})
                current_peak_amplitude = Summary.Peak(currentRows);
                [~, sort_ind] = sort(current_peak_amplitude, 'descend');
            else
                current_peak_amplitude = Summary.Trough(currentRows);
                [~, sort_ind] = sort(current_peak_amplitude, 'ascend');
            end

            %Store resorted raster and activity for plotting
            resorted_raster = [resorted_raster; current_raster(sort_ind,:)];
            activity_type = [activity_type; Summary.ResponseType(currentRows)];
        end
        
        groupActivity{g,L} = activity_type;
        groupRasters{g,L}(:,:) = resorted_raster;
        groupAverages{g,L} = average_raster;
        groupSTD{g,L} = raster_SEM;
    end
end

%%   ----Pie chart with ALL cells (not just active ones) ----

for f = 1:2 %Plot 2 times, once with % and once with count labels

    figure; tiledlayout(length(Groups),length(Loco))

    count = 0;
    for g = 1:length(Groups)
        for L = 1:length(Loco)
            count = count +1;

            T = tabulate(groupActivity{g,L});
            
            if isempty(T)
                continue;
            else
                T = cell2table(T);
            end

            noneInd = find(strcmp(T.T1, 'none'));
            actInd = find(strcmp(T.T1, 'activated'));
            proInd = find(strcmp(T.T1, 'prolonged'));
            supInd = find(strcmp(T.T1, 'suppressed'));
            finalInd = [noneInd; actInd; proInd; supInd];

            X = T.T2(finalInd);
            finallabels = T.T1(finalInd);
            
            %add percent or count to the labels
            for ff = 1:length(finallabels)
                if f == 1 %Percent
                    finallabels{ff} = [finallabels{ff} ' (' num2str(round(T.T3(finalInd(ff)),1)) '%)'];
                elseif f == 2 %Count
                    finallabels{ff} = [finallabels{ff} ' (' num2str(T.T2(finalInd(ff))) ')'];
                end
            end

            nexttile(count)
            pie(X,finallabels)
            temp_title = strcat(GroupsLabel{g},{' '},Loco{L});
            title(temp_title,'FontSize',10)
        end
    end
end

%% generate colors for each group
if cb >0
    if size(colorsForGraphs,2) < height(Groups)
        for i = 1:height(Groups)
            ctemp = flipud(cbrewer('qual','Set1',3));
            if length(Activity) == 4
                ctemp(4,:) = [0 , 0, 0]; % gray for "none" condition
            end
            Graph_Colors{i} = ctemp;
        end
    else
        if cb ==1
            for i = 1:height(Groups)
                ctemp = flipud(cbrewer('seq',colorsForGraphs{i},3)); %  activity types (activated, prolonged, suppressed)
                if length(Activity) == 4
                    ctemp(4,:) = [0 , 0, 0]; % gray for "none" condition
                end
                Graph_Colors{i} = ctemp;
            end
        elseif cb ==2
            for i = 1:height(Groups)
                ctemp = TAK_lab_color(colorsForGraphs{i}); %  activity types (activated, prolonged, suppressed)
                if length(Activity) == 4
                    ctemp(4,:) = [0 , 0, 0]; % gray for "none" condition
                end
                Graph_Colors{i} = ctemp;
            end
            
        end
    end
end
%% PLOT FIGURE BY GROUP

for g = 1:length(Groups)
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]); hold on
    t = tiledlayout(length(Loco),3);
    count = 0;
    for L = 1:length(Loco)
        count = count +1;
        %Skip Loco conditions with no active cells
        if isempty(groupActivity{g,L})
            continue
        end
        %--------------------------------------------------------------
        ax1 = nexttile(count); % avergage curves tile
        lns = plot(groupAverages{g,L}', 'LineWidth', 2);
        if Ops.plotShadedErrorBars
            [lenY, lenX] = size(groupAverages{g,L});
            for a = 1:lenY
                shadedErrorBar(1:lenX, groupAverages{g,L}(a,:), groupSTD{g,L}(a,:));
            end
        end
        
        title(strcat({'Average activity for '}, Loco{L}));
        ax1.Colormap = colororder(Graph_Colors{g});
        set(ax1, 'XTick', xtick)
        set(ax1, 'XTickLabel', xticklabel)
        legend(lns, Activity)
        xlim([0 nFrames])
        ylim([-1 2.5]);
        vline(soundtime, 'k')
        xlabel('Seconds')
        ylabel('z-score')
        
        %Plot stacked chart ----------------------
        count = count + 1;
        ax2 = nexttile(count);
        T = tabulate(string(groupActivity{g,L}));
        
        %Check that all the activity types are included in bar graph
        newT = cell(length(Activity),3);
        for a = 1:length(Activity)
            newT{a,1} = Activity{a};
            activityInd = find(strcmp(Activity{a}, T(:,1)));
            if ~isempty(activityInd)
                newT{a,2} = T{activityInd,2};
                newT{a,3} = T{activityInd,3};
            else
                newT{a,2} = 0;
                newT{a,3} = 0;
            end
        end
        
        bardata = ([newT{:,2}]);
        bar([1;nan], [bardata; nan(size(bardata))], 'stacked')
        title(strcat({'N Cells for '}, Loco{L}));
        ax2.Colormap = colororder(Graph_Colors{g});
        set(ax2, 'YDir','reverse')
        set(ax2, 'XTick', []) %Turn x axis off
        ylim([0 sum(bardata)]) %Make bar graph full height of chart
        xlim([0.75 1.25])
        ylabel('Cells')
        
        % Rasters -----------------------------------------
        count = count +1;
        ax3 = nexttile(count); % raster
        if Ops.smooth_rasters
            yyy = groupRasters{g,L}(:,1:end-1);
            for yy = 1:size(yyy,1)
                yyy(yy,:) = smooth(yyy(yy,:),10);
            end
            imagesc(yyy);
        else
            imagesc(groupRasters{g,L}(:,1:end-1));
        end
        ylabel('Cells')
        h = colorbar;
        h.Title.String = 'Z';
        caxis([-1 2]);

        title(Loco{L});
        set(ax3, 'XTick', xtick)
        set(ax3, 'XTickLabel', xticklabel)
        xlim([0 nFrames])
        vline(soundtime, 'k')
        hline(cumsum(bardata)+0.5, 'k') %Plot lines at activity borders
        xlabel('Seconds')
        colormap(bluewhitered(256)); %TODO use: colormap(bluewhitered_TAKlab)
    end
    title(t,GroupsLabel{g});
    
end

%% PLOT FIGURE BY ACTIVITY TYPE

fig2 = figure('units','normalized','outerposition',[0 0 1 1]); hold on
tiledlayout(length(Loco)*2,length(Activity));
for a = 1:length(Activity)
    count = a-length(Activity);
    for L = 1:length(Loco)
        
        %----------------average curves-----------------------------------
        count = count+length(Activity);
        nexttile(count)
        for g = 1:height(Groups)
            ln(g) = plot(groupAverages{g,L}(a,:), 'LineWidth', 2); hold on
            ln(g).Color = Graph_Colors{1,g}(a,:);
            if Ops.plotShadedErrorBars
                shadedErrorBar(1:lenX, groupAverages{g,L}(a,:), groupSTD{g,L}(a,:));
            end
        end
        
        temp_title = strcat(Loco{L},{' '},Activity(a));
        
        %         colororder(Graph_Colors{1,g}(1,:))
        legend(ln, GroupsLabel)
        set(gca, 'XTick', xtick)
        set(gca, 'XTickLabel', xticklabel)
        xlim([0 nFrames])
        vline(soundtime, 'k')
        xlabel('Seconds')
        ylabel('z-score')
        title([temp_title ' ' 'z-score'])
        
        
        %------------------- raster --------------------------------
        count = count+length(Activity);
        nexttile(count)
        grp_Raster = [];
        grp_N = [];
        for g = 1:height(Groups)
            grp_Activity  = strcmp(groupActivity{g,L}, Activity{a})';
            grp_Raster = [grp_Raster; groupRasters{g,L}(grp_Activity,1:end-1)];
            grp_N = [grp_N; size(groupRasters{g,L}(grp_Activity,1:end-1),1)];
        end
        if Ops.smooth_rasters
            yyy = grp_Raster;
            for yy = 1:size(yyy,1)
                yyy(yy,:) = smooth(yyy(yy,:),10);
            end
            imagesc(yyy);
        else
            imagesc(grp_Raster);
        end
        ylabel('Cells')
        h = colorbar;
        h.Title.String = 'Z';
        caxis([-1, 2]);
%         colormap(bluewhitered_TAKlab)
        set(gca, 'XTick', xtick)
        set(gca, 'XTickLabel', xticklabel)
        xlim([0 nFrames])
        vline(soundtime, 'k')
        hline(cumsum(grp_N)+0.5, 'k') %Plot lines at activity borders
        xlabel('Seconds')
        colormap(bluewhitered(256)); %TODO use: colormap(bluewhitered_TAKlab)
    end
end

end