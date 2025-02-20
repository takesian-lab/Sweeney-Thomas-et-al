function plot_simple_ExtractedData_SAM_APS(ExtractedData, UserOps)
% Plot output of SAM stim analysis for ExtractedData by activated/prolonged/suppressed

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%Return if not SAM data
if ExtractedData.StimInfo.StimProtocol ~= 5
    return;
end

disp('Plotting SAM data')

Activity = {'activated', 'prolonged', 'suppressed'};

if nargin < 2
    UserOps = struct;
    UserOps.Loco                = 'All'; %All, Running, or NotRunning
    UserOps.ResponseType        = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
    UserOps.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
    UserOps.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
    UserOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1)
    UserOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
    UserOps.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    UserOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    UserOps.plotCorrelations	= 0; %0 = don't plot dependent variable correlations, 1 = plot
    UserOps.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
    UserOps.FibPhotChannel      = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
    UserOps.SuppressOutput      = 1; %0 to print Ops to command line, 1 to suppress
end

UserOps.ResponseType = ''; %Regardless of what user supplies, take all ResponseTypes
StimData = simple_subset_ExtractedData(ExtractedData, UserOps);

Ops = StimData.SubsetOps;
StimInfo = StimData.StimInfo;
Summary = StimData.Summary;
StimAnalysis = StimData.StimAnalysis;
GroupList = Summary.Group;
Groups = unique(GroupList);
ActivityList = Summary.ResponseType;

if isempty(Summary)
    disp('No responsive cells found to plot for SAM')
    return
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

%Return if stim is not a matrix (i.e. when we look at SAM 0% only)
if isequal(size(StimInfo.Combined), [1,2])
    return
end

%% Plot Rate measures

variablesToPlot = {'BestRate','BestRate_fit_halfwidth', 'RSI', 'dprime', 'Sparseness', 'R2'};
responseToPlot = {'BestDepth', 'BestDepth', 'BestDepth', 'BestDepth', 'BestDepth', 'BestDepth'}; %For heatmap (All depth rows)
        
%Remove variables not in table (e.g. R2 might not be included)
variablesToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];
responseToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];       

stackedHeatMaps = {};
stackedCount = {};

for g = 1:length(Groups)

    figure('units','normalized','outerposition',[0 0 1 1]); hold on
    tiledlayout(length(Activity)*2,length(variablesToPlot))
    
    for a = 1:length(Activity)
        isGroup = strcmp(Groups{g},GroupList);
        isActivity = strcmp(Activity{a},ActivityList);
        isGroupandActivity = (isGroup + isActivity) == 2;
        tempStim = StimAnalysis(isGroupandActivity,:);
        
        response_mat = tempStim.Response;

        %Flip sign of responses for suppressed stim for heat map
        if strcmp(Activity{a},'suppressed')
            response_mat = -response_mat;
            subtitle_tag = ' - sign flipped';
        else
            subtitle_tag = '';
        end
        
        %Normalize response_mats
        max_resp = max(response_mat,[],2);
        norm_response_mat = response_mat./max_resp;

        for v = 1:length(variablesToPlot)
            vname = regexprep(variablesToPlot{v},'_',' ','emptymatch');

            %HISTOGRAM
            nexttile((a-1)*2*length(variablesToPlot) + v); hold on
            y = tempStim.(variablesToPlot{v});
            if isempty(y); continue; end
            histogram(y)%,'facecolor','g','facealpha',0.2,'edgecolor','none');
            vline(mean(y, 'omitnan'), 'r')
            vline(mean(y, 'omitnan') - std(y))
            vline(mean(y, 'omitnan') + std(y))
            title(Activity{a})
            subtitle(strcat({'Mean '}, vname, {' = '}, num2str(mean(y, 'omitnan'))))
            xlabel(vname)
            ylabel('N cells')
            legend(['N = ' num2str(length(y))])

            %HEAT MAPS
            nexttile((a-1)*2*length(variablesToPlot) + length(variablesToPlot) + v); hold on

            %For heatmat, get values represented in responseToPlot
            %i.e. instead of plotting all 64 responses in RF, just plot the intensity row or frequency column specified
            [respmatind, respmat] = deal(nan(length(y),length(StimInfo.V1)));
            if strcmp(responseToPlot{v}, 'BestDepth')
                V2 = tempStim.(responseToPlot{v});
                %Fill respmat with responses for that intensity row
                for r = 1:size(V2,1)
                    respmatind(r,:) = find(StimInfo.Combined(:,2) == V2(r));
                    respmat(r,:) = norm_response_mat(r,respmatind(r,:));
                    
                    %Special for SAM, fill 0% rows/columns (otherwise it shows up as NaN in the heat maps)
                    zeroInd = find(((StimInfo.Combined(:,1) == 0) + (StimInfo.Combined(:,2) == 0)) == 2);
                    if ~isempty(zeroInd)
                        if V2(r) == 0
                            %If BestDepth is 0, take average of responses at other rates
                            tempRF = nan(size(StimInfo.RF_V1));
                            RFind = sub2ind(size(tempRF), StimInfo.RF_Map(:,1), StimInfo.RF_Map(:,2));
                            tempRF(RFind) = norm_response_mat(r,:);
                            %Remove zero row and column
                            tempRF(:,StimInfo.V1 == 0) = [];
                            tempRF(StimInfo.V2 == 0,:) = [];
                            %Get average at all other rates
                            respmat(r, StimInfo.V1 ~= 0) = mean(tempRF,1,'omitnan');
                        else
                           %If BestDepth is not 0, fill in response to 0% at missing value
                           respmatind(r,StimInfo.V1==0) = zeroInd;
                           respmat(r,:) = norm_response_mat(r,respmatind(r,:));
                        end
                    end
                end
            else
                error('Response type has not been coded for yet')
            end

            [~, sortind] = sort(y, 'ascend');
            sortedmat = respmat(sortind,:);
            imagesc(sortedmat)
            xlim([0.5 size(sortedmat,2)+0.5])
            ylim([0.5 size(sortedmat,1)+0.5])
            set(gca, 'XTick', 1:size(sortedmat,2))
            set(gca, 'XTickLabel', StimInfo.V1) 
            h = colorbar;
            h.Title.String = 'AUC';
            caxis([0 1]); %Normalized responses should be between 0 to 1
            xlabel('Rate (Hz)')
            ylabel('N cells')
            title([responseToPlot{v} ' sorted by ' vname subtitle_tag])
            
            %Store stacked heat maps
            if a == 1
                stackedHeatMaps{g,v} = flipud(sortedmat);
                stackedCount{g,v} = size(sortedmat,1);
            else
                stackedHeatMaps{g,v} = [stackedHeatMaps{g,v}; flipud(sortedmat)];
                stackedCount{g,v} = [stackedCount{g,v}; size(sortedmat,1)];
            end
        end
    end

    tak_suptitle([StimInfo.StimType ' ' Groups{g}])
    drawnow
end

%% STACKED HEAT MAPS
for g = 1:length(Groups)
    figure('units','normalized','outerposition',[0 0 1 1]); hold on
    tiledlayout(1,length(variablesToPlot))
    
    for v = 1:length(variablesToPlot)
        vname = regexprep(variablesToPlot{v},'_',' ','emptymatch');

        nexttile; hold on
        imagesc(flipud(stackedHeatMaps{g,v}))
        %Draw white line showing border between activity types
        %Imagesc plots lines from bottom up, so need to flip count
        count = flipud(stackedCount{g,v});
        for a = 1:length(Activity)-1
            N = sum(count(1:a));
            hline(N,'w')
        end
        xlim([0.5 size(stackedHeatMaps{g,v},2)+0.5])
        ylim([0.5 size(stackedHeatMaps{g,v},1)+0.5])
        set(gca, 'XTick', 1:size(stackedHeatMaps{g,v},2))
        set(gca, 'XTickLabel', StimInfo.V1) 
        h = colorbar;
        h.Title.String = 'AUC';
        caxis([0 1]); %Normalized responses should be between 0 to 1
        xlabel('Rate (Hz)')
        ylabel('N cells')
        title([responseToPlot{v} ' sorted by ' vname])
    end
    tak_suptitle([StimInfo.StimType ' ' Groups{g}])
    drawnow
end

%% Plot everything else

variablesToPlot = {'BestDepth', 'BestDepth_fit_halfwidth', 'DSI', 'dprime', 'Sparseness', 'R2'};
responseToPlot = {'BestRate', 'BestRate', 'BestRate', 'BestRate', 'BestRate', 'BestRate'}; %For heatmap (All BF freq column)
        
%Remove variables not in table (e.g. R2 might not be included)
variablesToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];
responseToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];       

stackedHeatMaps = {};
stackedCount = {};

for g = 1:length(Groups)
    
    figure('units','normalized','outerposition',[0 0 1 1]); hold on
    tiledlayout(length(Activity)*2,length(variablesToPlot))

    for a = 1:length(Activity)
        isGroup = strcmp(Groups{g},GroupList);
        isActivity = strcmp(Activity{a},ActivityList);
        isGroupandActivity = (isGroup + isActivity) == 2;
        tempStim = StimAnalysis(isGroupandActivity,:);

        response_mat = tempStim.Response;

        %Flip sign of responses for suppressed stim for heat map
        if strcmp(Activity{a},'suppressed')
            response_mat = -response_mat;
            subtitle_tag = ' - sign flipped';
        else
            subtitle_tag = '';
        end
        
        %Normalize response_mats
        max_resp = max(response_mat,[],2);
        norm_response_mat = response_mat./max_resp;

        for v = 1:length(variablesToPlot)
            vname = regexprep(variablesToPlot{v},'_',' ','emptymatch');

            %HISTOGRAM
            nexttile((a-1)*2*length(variablesToPlot) + v); hold on
            y = tempStim.(variablesToPlot{v});
            if isempty(y); continue; end
            histogram(y)%,'facecolor','g','facealpha',0.2,'edgecolor','none');
            vline(mean(y, 'omitnan'), 'r')
            vline(mean(y, 'omitnan') - std(y))
            vline(mean(y, 'omitnan') + std(y))
            title(Activity{a})
            subtitle(strcat({'Mean '}, vname, {' = '}, num2str(mean(y, 'omitnan'))))
            xlabel(vname)
            ylabel('N cells')
            legend(['N = ' num2str(length(y))])

            %HEAT MAPS
            nexttile((a-1)*2*length(variablesToPlot) + length(variablesToPlot) + v); hold on

            %For heatmat, get values represented in responseToPlot
            %i.e. instead of plotting all 64 responses in RF, just plot the intensity row or frequency column specified
            [respmatind, respmat] = deal(nan(length(y),length(StimInfo.V2)));
            if strcmp(responseToPlot{v}, 'BestRate')
                V1 = tempStim.(responseToPlot{v});
                %Fill respmat with responses for that intensity row
                for r = 1:size(V1,1)
                    respmatind(r,:) = find(StimInfo.Combined(:,1) == V1(r));
                    respmat(r,:) = norm_response_mat(r,respmatind(r,:));
                    
                    %Special for SAM, fill 0% rows/columns (otherwise it shows up as NaN in the heat maps)
                    zeroInd = find(((StimInfo.Combined(:,1) == 0) + (StimInfo.Combined(:,2) == 0)) == 2);
                    if ~isempty(zeroInd)
                        if V1(r) == 0
                            %If BestDepth is 0, take average of responses at other rates
                            tempRF = nan(size(StimInfo.RF_V1));
                            RFind = sub2ind(size(tempRF), StimInfo.RF_Map(:,1), StimInfo.RF_Map(:,2));
                            tempRF(RFind) = norm_response_mat(r,:);
                            %Remove zero row and column
                            tempRF(:,StimInfo.V1 == 0) = [];
                            tempRF(StimInfo.V2 == 0,:) = [];
                            %Get average at all other rates
                            respmat(r, StimInfo.V2 ~= 0) = mean(tempRF,2,'omitnan');
                        else
                           %If BestDepth is not 0, fill in response to 0% at missing value
                           respmatind(r,StimInfo.V2==0) = zeroInd;
                           respmat(r,:) = norm_response_mat(r,respmatind(r,:));
                        end
                    end
                end
            else
                error('Response type has not been coded for yet')
            end

            %flip responses L/R so that depth is ordered from 0 to 100
            respmat = fliplr(respmat);
            
            [~, sortind] = sort(y, 'ascend');
            sortedmat = respmat(sortind,:);
            imagesc(sortedmat)
            xlim([0.5 size(sortedmat,2)+0.5])
            ylim([0.5 size(sortedmat,1)+0.5])
            set(gca, 'XTick', 1:size(sortedmat,2))
            set(gca, 'XTickLabel', flipud(StimInfo.V2)) 
            h = colorbar;
            h.Title.String = 'AUC';
            caxis([0 1]); %Normalized responses should be between 0 to 1
            xlabel('Depths (%)')
            ylabel('N cells')
            title([responseToPlot{v} ' sorted by ' vname subtitle_tag])
            
            %Store stacked heat maps
            if a == 1
                stackedHeatMaps{g,v} = flipud(sortedmat);
                stackedCount{g,v} = size(sortedmat,1);
            else
                stackedHeatMaps{g,v} = [stackedHeatMaps{g,v}; flipud(sortedmat)];
                stackedCount{g,v} = [stackedCount{g,v}; size(sortedmat,1)];
            end
        end
    end
    
    tak_suptitle([StimInfo.StimType ' ' Groups{g}])
    drawnow
end

%% STACKED HEAT MAPS
for g = 1:length(Groups)
    figure('units','normalized','outerposition',[0 0 1 1]); hold on
    tiledlayout(1,length(variablesToPlot))
    
    for v = 1:length(variablesToPlot)
        vname = regexprep(variablesToPlot{v},'_',' ','emptymatch');

        nexttile; hold on
        imagesc(flipud(stackedHeatMaps{g,v}))
        %Draw white line showing border between activity types
        %Imagesc plots lines from bottom up, so need to flip count
        count = flipud(stackedCount{g,v});
        for a = 1:length(Activity)-1
            N = sum(count(1:a));
            hline(N,'w')
        end
        xlim([0.5 size(stackedHeatMaps{g,v},2)+0.5])
        ylim([0.5 size(stackedHeatMaps{g,v},1)+0.5])
        set(gca, 'XTick', 1:size(sortedmat,2))
        set(gca, 'XTickLabel', flipud(StimInfo.V2)) 
        h = colorbar;
        h.Title.String = 'AUC';
        caxis([0 1]); %Normalized responses should be between 0 to 1
        xlabel('Depths (%)')
        ylabel('N cells')
        title([responseToPlot{v} ' sorted by ' vname])
    end
    tak_suptitle([StimInfo.StimType ' ' Groups{g}])
    drawnow
end
