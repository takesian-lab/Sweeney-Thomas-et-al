function plot_simple_ExtractedData_FreqDisc(ExtractedData, UserOps)
% Plot behavior data (loco, whisker, pupil, etc)

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%Return if not behavior data
if ~ismember(ExtractedData.StimInfo.StimProtocol, [7]) %Freq disc behavior
    return;
end

%I RECOMMEND NOT CHANGING THESE AND USING USEROPS INSTEAD
if nargin < 2
    UserOps = struct;
    UserOps.Loco                = ''; %All, Running, or NotRunning
    UserOps.ResponseType        = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
    UserOps.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
    UserOps.IsResponsive        = 2; %1 for responsive (and reliable), 0 for not, 2 for both
    UserOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1)
    UserOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
    UserOps.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    UserOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    UserOps.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
    UserOps.FibPhotChannel      = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
    UserOps.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress
    UserOps.ByMouse             = 1; %Sort by unique mice. Currently only utilized in FreqDisc (stimProtocol =7)
    UserOps.plotShadedErrorBars = 1; %%0 = don't plot, 1 = plot shaded error bars
end

disp('Plotting Frequency Discrimination behavior data')

SubsetData = simple_subset_ExtractedData(ExtractedData, UserOps); %
Ops = SubsetData.SubsetOps;
BlockInfo = SubsetData.BlockInfo;
TrialData = SubsetData.TrialData;
StimInfo = SubsetData.StimInfo;
Summary = SubsetData.Summary;
GroupList = Summary.Group;
Groups = unique(GroupList);
MouseList = BlockInfo.MouseID;
Mice = unique(MouseList);
BehaviorData = ExtractedData.BehaviorData; 
%TODO: edit simple_subset_ExtractedData function to include BehaviorData and FreqDiscAnalysis to also be filtered 

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
nBaselineFrames = StimInfo.nBaselineFrames;
xInFrames = nBaselineFrames + StimInfo.after_inFrames;
xInSeconds = xInFrames/StimInfo.fs;
graphx = linspace(0, xInSeconds, xInFrames); %changes trial length (x) in 180 frames to 6 seconds
xtick = 0:nBaselineFrames:xInFrames;
xticklabels = xtick./StimInfo.fs;

if  UserOps.FibPhotChannel < 4
    channel = UserOps.FibPhotChannel;
else
    channel = 1:3;
end


%% Concat data for plotting


for M = 1:length(Mice)
    %%
    mIX = find(strcmp(Mice(M),MouseList)); % blocks that correspond to this mouse
    
    cat_outcome     = [];
    cat_trials      = [];
    cat_dates       = [];
    
    for j = 1:length(mIX)
        tO = TrialData(mIX(j)).Stim_V2; % outcome
        tT = TrialData(mIX(j)).Trials;
        dD = SubsetData.BlockInfo.Date(mIX(j));
        tD = repelem(dD,length(tO))';
        
        if ~isempty(tT) % some mice may be missing functional data from a given day
        cat_outcome     = [cat_outcome;tO];
        cat_trials      = cat(2,cat_trials,tT);
        cat_dates       = [cat_dates;tD];
        end
        

    end
    % Behavior data is index the same way as BlockInfo 
    cat_dprime = BehaviorData.dprime(mIX);
    cat_dprime_catchonly = BehaviorData.dprime_catchonly(mIX);
    
    % if there are results with NaN as the outcome, remove them here:
    remvnan = find(isnan(cat_outcome));
    cat_outcome(remvnan) = [];
    cat_trials  (:,remvnan,:) = [];
    cat_dates (remvnan) = [];
    
    result = strings(size(cat_outcome));
    result(cat_outcome == 1) = "Hit";
    result(cat_outcome == 0) = "Miss";
    result(cat_outcome == 3) = "Withhold";
    result(cat_outcome == 4) = "False Alarm";
    result(cat_outcome == 5) = "Catch";
    uqO = unique(result); uqD = unique(cat_dates);
    %     day_colors = cbrewer('div','RdBu',numel(uqD));
    day_colors = cbrewer('div', 'Spectral', numel(uqD));
    
    
    for Ch = 1:length(channel)

        % set up data for rasters:
        for i = 1:length(uqO)
            trialsIX = find(strcmp(result,uqO(i)));
            tempD = cat_dates(trialsIX);
            tuqd = unique(tempD);
            
            % find 
            allDays = unique(cat_dates);
            allDayForStack = zeros(1,length(allDays));
            
            % grab colors for graphing:
            ccIX = find(ismember(uqD,tuqd));
            storeColors{i} = day_colors(ccIX,:);
            
            storeRaster{i} =  squeeze(cat_trials(channel(Ch),trialsIX,:));
            
            % find average across days and number of trials across days:
            for d = 1:length(tuqd)
                Dd = find(strcmp(tempD,tuqd(d)));
                
                %             daystart = daystart + daystart+numel(Dd);
                %             dayline = [dayline; daystart];
                DayTrace(d,:) = mean(squeeze(cat_trials(channel(Ch),trialsIX(Dd),:)),1,'omitnan');
                DaySEM(d,:) = std(squeeze(cat_trials(channel(Ch),trialsIX(Dd),:)),1,'omitnan')./sqrt(length(trialsIX(Dd)));
                DayNumel(d) = numel(Dd);
                
                dayix = find(contains(allDays,tuqd(d)));
                allDayForStack(dayix) = numel(Dd);
                
            end
            %         All_dayline{i} = dayline;
            All_DayTrace{i} = DayTrace;
            All_DaySEM{i} = DaySEM;
            DaysforStack{i} = DayNumel;
            DaysforStack_allDays{i} = allDayForStack;
            clear DayTrace DayNumel DaySEM allDayForStack
        end
        
        %% Plot behavior
        
        % There will be one column for each result (hit, miss, etc) and two rows for each behavior measure (curves and rasters)
        h = figure; t = tiledlayout(4, length(uqO));
        allTracesStacked = vertcat(All_DayTrace{:});
        y_limits = [min(allTracesStacked(:)), max(allTracesStacked(:))];
        for r = 1:length(uqO)
            for rr = 1:4
                if rr ==1 % average curves per day
                    nexttile(r)
                    traces = All_DayTrace{r};
                    sems = All_DaySEM{r};
                    colors = storeColors{r};

                    for p = 1:size(traces,1)
                        % if add shaded error bar 
                        if UserOps.plotShadedErrorBars == 1
                            [hl, ~] = boundedline(double(graphx),double(smooth(traces(p,:))'),double(sems(p,:)),'alpha','transparency',0.1,...
                                'cmap',[colors(p,:)]); %returns the handles the resulting line and patch object(s).
                            hold on
                            hl.LineWidth = 1; % set the line thicker

                        else % plot average line only
                            plot(graphx, smooth(traces(p,:)),'Color',colors(p,:)); hold on 
                        end 
                        
                        ylabel('zscore')
                        xlabel ('time (s)')
                        xlim([graphx(1) graphx(end)])
                        title(uqO(r))
%                         ylim(y_limits)
                    end
                    
                    
                    
                elseif rr == 2 % rasters across days
                    nexttile (r+length(uqO))
                    
                    % smooth rasters
                    if UserOps.smooth_rasters
                        yyy = storeRaster{r}(:,1:end-1);
                        for yy = 1:size(yyy,1)
                            yyy(yy,:) = smooth(yyy(yy,:),10);
                        end
                        imagesc(yyy);
                    else
                        imagesc(storeRaster{r});
                    end
                    caxis([-2, 2]);
                    
                    
                    set(gca, 'XTick', xtick)
                    set(gca, 'XTickLabel', xticklabels)
                    vline(nBaselineFrames)
                    ylabel('Trials')
                    xlabel ('time')
                    colormap(bluewhitered(256)); %TODO use:colormap(bluewhitered_TAKlab)
                    %         hline(All_dayline{i});
                elseif   rr ==3
                    nexttile(r+(length(uqO)*2))
                    
                    bardata = flip([DaysforStack{r}])';
                    colormp = flip(colors);
                    %         bar([1;nan], [bardata; nan(size(bardata))], 'stacked')
                    b =  bar([1], [bardata], 'stacked','FaceColor','flat');
                    for k = 1:length(bardata)
                        b(k).CData = colormp(k,:);
                    end
                    
                else %rr == 4 % dprime (per day) for each trial of this outcome type
                    nexttile(r+(length(uqO)*3))

                   if strcmp(uqO(r),string('Catch'))% catch trials
                        thisdp = cat_dprime_catchonly;
                    else 
                        thisdp = cat_dprime;

                    end 
                    
                    bardata = DaysforStack_allDays{r};
                    thisDPrime_plot = [];
                    
                    for b = 1:length(bardata)
                        thisDPrime_plot = [thisDPrime_plot,ones(1,bardata(b))* thisdp(b)];
                    end
                    
                    
                    plot(thisDPrime_plot,[1:length(thisDPrime_plot)],'LineWidth',2,'Color','k');
                    set(gca,'YDir', 'reverse')
                    if strcmp(uqO(r),string('Catch')) % catch trials
                        xlabel('dprime: only catch trials')
                    else 
                        xlabel('dprime')
                    end 
                    clear thisDPrime_plot thisdp bardata
                    xlim([0 4])
                    
                end
                
            end
            
        end
    end
    title(t,strcat(Mice(M),(' '),'Channel',num2str(channel(Ch))));
    
    
end% loop through mice
end







