function plot_simple_ExtractedData_MaryseBehavior(ExtractedData, UserOps)

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%Return if not Maryse behavior data
if ExtractedData.StimInfo.StimProtocol ~= 13
    return;
end

return

disp('Plotting Maryse behavior data')

% SCRIPT NOT FINISHED

StimData = simple_subset_ExtractedData(ExtractedData);

%% Options

sortbyGCAMP = 0; %0 for groups, 1 for gcamp, 2 to combine
sortbyCondition = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
sortbyRedCell = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
FibPhotChannel = 2; %1 = blue, 2 = green, 3 = red

%% Get data from ExtractedData

CellList = ExtractedData.CellList;
CellData = ExtractedData.CellData.All;
CellDataByStim = ExtractedData.CellDataByStim.All;
FinalTraces = ExtractedData.FinalTraces.All;
StimAnalysis = ExtractedData.StimAnalysis.All;
%PeakData = [ExtractedData.CellDataByStim.All.PeakData];
Activity = {'activated', 'prolonged', 'suppressed', 'none'};% 'none'};

%X labels for figures
xlabel_increment = 0.5; %seconds
fs = ExtractedData.StimInfo.fs;
nFrames = ExtractedData.StimInfo.nFrames;
increment_in_frames = xlabel_increment*fs;
xtick = 0:increment_in_frames:nFrames;
xticklabel = 0:xlabel_increment:((1/fs)*nFrames);
soundtime = ExtractedData.StimInfo.baseline*fs;

%% Sort by Group vs. GCaMP vs. Condition

if sortbyGCAMP == 0
    GroupList = CellList.Group;
    if sortbyCondition > 0
        GroupList = split(GroupList,'_');
        GroupList = GroupList(:,sortbyCondition);
    end 
elseif sortbyGCAMP == 1
    GroupList = CellList.GCaMP;
elseif sortbyGCAMP == 2
    GroupList = CellList.Group;
    if sortbyCondition > 0
        GroupList = split(GroupList,'_');
        GroupList = GroupList(:,sortbyCondition);
    end 
    GroupList = strcat(GroupList, '_', CellList.GCaMP);
end

Groups = unique(GroupList);

%% Sort by FibPhot Channel or Red Cell

%Determine if FibPhot data
hasFibPhot = any(strcmp('FibPhot', ExtractedData.BlockInfo.AnalysisPath));
FibPhotChannels = ["blue", "green", "red"];
FibPhotChannelToAnalyze = FibPhotChannels{FibPhotChannel};

if hasFibPhot
    channels = ExtractedData.CellList.Channel;
    keep = strcmp(channels, FibPhotChannelToAnalyze);
else
    redCell = CellList.RedCell;

    if sortbyRedCell == 0
        keep = logical(ones(size(redCell)));
    elseif sortbyRedCell == 1
        keep = logical(redCell);
    elseif sortbyRedCell == 2
        keep = logical(~redCell);
    end
end
    
GroupList = GroupList(keep,:);
CellData = CellData(keep,:);
CellDataByStim = CellDataByStim(keep);
StimAnalysis = StimAnalysis(keep,:);
FinalTraces = FinalTraces(keep,:);

%% Plot rasters and figure out percentage of each activity type for each group
% Use final trace results for this (??)

percent_mat = zeros(length(Groups),length(Activity));
n_mat = zeros(length(Groups),length(Activity));

for g = 1:length(Groups)
    group_ind = strcmp(GroupList, Groups{g});
    
    activityList = CellData.ResponseType(group_ind);
    T = tabulate(activityList);
    
    for a = 1:length(Activity)
        if ~any(ismember(T(:,1), Activity{a}))
            continue
        end
        ind = strcmp(T(:,1), Activity{a});
        percent_mat(g,a) = T{ind,3};
        n_mat(g,a) = T{ind,2};
    end
end

%% RASTER AND PIE CHART
group_sort_idx = {};
%cellrows = 1:1286;

figure
for g = 1:length(Groups)
    group_ind = strcmp(GroupList, Groups{g});
    group_raster = [];
    group_rows = [];
    
    for a = 1:length(Activity)
        if strcmp(Activity{a}, 'none')
            continue;
        end
        act_ind = strcmp(CellData.ResponseType, Activity{a});
        overlap = (group_ind + act_ind) == 2;
        %cellrows_overlap = cellrows(overlap);
        
        currentTraces = FinalTraces(overlap,:);
        
        %sort traces and store in raster
        if strcmp(Activity{a}, 'suppressed')
            currentAUC = CellData.Trough_AUC(overlap);
        else
            currentAUC = CellData.Peak_AUC(overlap);
        end
        [~, sortind] = sort(currentAUC);
        currentTraces = currentTraces(sortind,:);
        %cellrows_sorted = cellrows_overlap(sortind);
        
        group_raster = [group_raster; currentTraces];
        %group_rows = [group_rows; cellrows_sorted'];
        
        Y = smooth(mean(currentTraces,1,'omitnan'));
        
        %PLOT MEAN TRACE
        subplot(4,length(Groups),g); hold on
        if size(Y,2) > 1
            shadedErrorBar(1:size(currentTraces,2), Y, smooth(std(currentTraces,[],1)./sqrt(size(currentTraces,1)-1)))
        end
        plot(Y);
        %ylim([-0.5 1])
        set(gca, 'XTick', xtick)
        set(gca, 'XTickLabel', xticklabel)
        xlim([0 nFrames])
        vline(soundtime, 'k')
        xlabel('Seconds')
        ylabel('z-score')
        title(Groups{g})
    end

    %PLOT RASTER
    subplot(4,length(Groups),g+length(Groups))
    group_raster_smoothed = group_raster;
    for i = 1:size(group_raster_smoothed,1)
        group_raster_smoothed(i,:) = smooth(group_raster(i,:));
    end
    imagesc(group_raster_smoothed)
    ylim([1 size(group_raster,1)])
    ylabel('Cells')
    set(gca, 'XTick', xtick)
    set(gca, 'XTickLabel', xticklabel)
    xlim([0 nFrames])
    vline(soundtime, 'k')
    %hline(cumsum(bardata)+0.5, 'k') %Plot lines at activity borders
    xlabel('Seconds')
    %colormap(gca, bluewhitered(256))
    c = colorbar;
    c.Title.String = 'Z';
    %caxis([0 1])
    
    %PLOT PIECHART
    subplot(4,length(Groups),g+(length(Groups)*2))
    pie(percent_mat(g,:),Activity)
    
    %PLOT COUNT
    subplot(4,length(Groups),g+(length(Groups)*3))
    bar(n_mat(g,:))
    text((1:length(Activity))-0.3, n_mat(g,:) + 7, num2str(n_mat(g,:)'))
    ylabel('N cells')
    set(gca, 'XTickLabel', Activity)
    
    group_sort_idx{g} = group_rows;
end

% TODO: Replace passive rasters with rasters sorted by active stim condition

%% Plots by stim condition

StimInfo =  ExtractedData.StimInfo;
V1 = StimInfo.V1;
V2 = [1 0 5];
V2_label = {'Go', 'No-Go', 'Catch'};
combined_stim = StimInfo.Combined;

stim_rasters = cell(length(V1),3)';

for s = 1:length(V1)
    for ss = 1:length(V2)
    
        stim_row = ((combined_stim(:,1) == V1(s)) + (combined_stim(:,2) == V2(ss))) == 2;
        temp_raster = [];
        
        for c = 1:size(CellDataByStim,2) %Loop through all channels
            if ~isempty(CellDataByStim(c).StimTracesAveraged)
                temp = CellDataByStim(c).StimTracesAveraged(stim_row,:);
                if ~all(isnan(temp))
                    temp_raster = [temp_raster; temp];
                end
            end
        end
        
        stim_rasters{ss,s} = temp_raster;
    end
end

figure('units','normalized','outerposition',[0 0 1 1]); hold on

for s = 1:length(V1)
    for ss = 1:length(V2)
        
        temp_raster = stim_rasters{ss,s};
        Y_n = sum(~all(isnan(temp_raster),2));
        if Y_n <=1
            Y_n = 1;
        end
        Y = mean(temp_raster,1,'omitnan');    
        Y_sem = std(temp_raster,[],1,'omitnan')./sqrt(Y_n-1);
        
        subplot(2*length(V2),length(V1),length(V1)*2*(ss-1) + s); hold on
        if ~isempty(Y)
            plot(Y)
            shadedErrorBar(1:length(Y), Y, Y_sem)
        end
        ylim([-1 10])
        xlim([0 nFrames])
        vline(soundtime, 'k')
        if ismember(s, 1:length(V1):length(V1)*4)
            ylabel(V2_label{ss})
        end
        set(gca, 'XTick', xtick)
        set(gca, 'XTickLabel', xticklabel)
        
        if ss == 1
            title([num2str(V1(s)) ' oct.'])
        end
        
        subplot(2*length(V2),length(V1),length(V1)*2*(ss-1) + s + length(V1)); hold on
        imagesc(temp_raster)
        xlim([0 nFrames])
        if ~isempty(temp_raster)
            ylim([0.5 size(temp_raster,1)+0.5])
        end
        vline(soundtime, 'k')
        if ismember(s, 1:length(V1):length(V1)*4)
            ylabel('Sessions')
        end
        set(gca, 'XTick', xtick)
        set(gca, 'XTickLabel', xticklabel)
        if (s + length(V1)) > length(V1)*4
            xlabel('Seconds')
        end
    end
end

%% TODO: Plot neurometric rasters
% 
% RF = ExtractedData.RF;
% absPeakAvg = squeeze(RF.F.absPeakAvg)';
% RFcolumns = RF.ColumnHeaders;
% MI_ind = find(strcmp(RFcolumns, 'Modulation Index'));
% MI = RF.F.data(:,MI_ind);
% 
% for g = 1:length(Groups)
%     currentCells = strcmpi(groupList,Groups(g));
%     
%     peak_F = [];
%     peak_S = [];
%     MI_F = [];
%     MI_S = [];
%     activity_type = [];
%         
%     
%     activeRows = ~strcmpi(activityList, 'none');
%     currentRows = and(currentCells, activeRows);
%     currentRows = rows(currentRows);
%     cellNumbers = find(currentRows); %cell numbers to use for comparison in locomotor 
% 
%     current_peak_F = absPeakAvg(currentRows,:);
%     current_MI = MI(currentRows);
%     
%     if g == 1
%         [b1, b] = sort(current_MI);
%         b1 = flipud(b1);
%         b = flipud(b);
%         
%         NDNF = current_peak_F;
% 
%     elseif g == 2
%         [c1, c] = sort(current_MI);
%         c1 = flipud(c1);
%         c = flipud(c);
%         
%         Pyr = current_peak_F;
%     end
% end
% 
% % GO TO EXCEL AND FIX NDNF and Pyr variables because of diff. freqs
% sortedNDNF = NDNF(b,:);
% sortedNDNF = sortedNDNF./(max(sortedNDNF,[],2));
%    
% sortedPyr = Pyr(c,:);
% sortedPyr = sortedPyr./(max(sortedPyr,[],2));
% 
% b1(b1 < -1) = nan;
% b1(b1 > 1) = nan;
% 
% c1(c1 < -1) = nan;
% c1(c1 > 1) = nan;
% 
% figure;
% subplot(1,6,1)
% imagesc(sortedNDNF)
% %g1 = colorbar;
% title('NDNF')
% ylabel('Cells')
% 
% subplot(1,6,2)
% imagesc(b1)
% h1 = colorbar;
% title('MI')
% ylim(h1,[-1 1])
% 
% subplot(1,6,3); hold on
% boxplot(b1)
% errorbar(nanmean(b1), std(b1)/(sqrt(length(b1)-1)));
% scatter(1, nanmean(b1))
% ylim([-.8 1])
% title(['Avg MI: ' num2str(nanmean(b1))])
% 
% subplot(1,6,4)
% imagesc(sortedPyr)
% title('Pyr')
% ylabel('Cells')
% 
% subplot(1,6,5)
% imagesc(c1)
% h2 = colorbar;
% title('MI')
% ylim(h2,[-1 1])
% 
% subplot(1,6,6); hold on
% boxplot(c1)
% errorbar(nanmean(c1), std(c1)/(sqrt(length(c1)-1)));
% scatter(1, nanmean(c1))
% %scatter(ones(length(c1)),c1)
% ylim([-.8 1])
% title(['Avg MI: ' num2str(nanmean(c1))])
% 
% colormap('hot')
% 
%         
% 
% normMax = data.neurometric_curve_max./max(data.neurometric_curve_max,[],2);
% normMean = data.neurometric_curve_mean./max(data.neurometric_curve_mean,[],2);
% 
% [~, max_I] = sort(data.modIndexMax);
% sorted_max = data.neurometric_curve_max(max_I,:);
% sorted_max = normMax(max_I,:);
% [~, mean_I] = sort(data.modIndexMean);
% sorted_mean = data.neurometric_curve_mean(mean_I,:);
% sorted_mean = normMean(mean_I,:);
% 
% 
% figure;
% 
% subplot(2,2,1)
% imagesc(sorted_max)
% ylabel('Sorted cells')
% set(gca, 'XTick', 1:length(alternating))
% set(gca, 'XTickLabel', alternating)
% c = colorbar;
% title('Max neurometric curve')
% ylabel(c, 'Max df/f response')
% 
% subplot(2,2,3); hold on
% scatter(ones(1,length(data.modIndexMax)), data.modIndexMax)
% ylim([-1 1])
% hline(0)
% hline(nanmean(data.modIndexMax), 'k')
% title(['Average mod index: ' num2str(nanmean(data.modIndexMax))])
% ylabel('Modulation index')
% 
% subplot(2,2,2)
% imagesc(sorted_mean)
% ylabel('Sorted cells')
% set(gca, 'XTick', 1:length(alternating))
% set(gca, 'XTickLabel', alternating)
% c = colorbar;
% title('Mean neurometric curve')
% ylabel(c,'Mean df/f response')
% 
% subplot(2,2,4); hold on
% scatter(ones(1,length(data.modIndexMean)), data.modIndexMean)
% ylim([-1 1])
% hline(0)
% hline(nanmean(data.modIndexMean), 'k')
% title(['Average mod index: ' num2str(nanmean(data.modIndexMean))])
% ylabel('Modulation index')
