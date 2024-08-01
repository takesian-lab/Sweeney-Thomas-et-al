function plot_simple_ExtractedData_FM(ExtractedData, UserOps)
% Plot output of FM stim analysis for ExtractedData

%Version control:
%V1: archived 2/1/2023
%V2: current version

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%Return if not FM data
if ExtractedData.StimInfo.StimProtocol ~= 3
    return;
end

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

disp('Plotting FM data')

StimData = simple_subset_ExtractedData(ExtractedData, UserOps);

Ops = StimData.SubsetOps;
StimInfo = StimData.StimInfo;
Summary = StimData.Summary;
StimAnalysis = StimData.StimAnalysis;
GroupList = Summary.Group;
Groups = unique(GroupList);

if isempty(Summary)
    disp('No responsive cells found to plot for FM')
    return
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

% Reshape stim list
speeds = StimInfo.Combined(:,1).*StimInfo.Combined(:,2);
if ~issorted(speeds(speeds > 0),'ascend') || ~issorted(speeds(speeds < 0), 'descend')
    error('Speeds are in unexpected format for reordering')
end
reorder_idx = [flipud(find(speeds < 0)); find(speeds > 0)];
speeds = speeds(reorder_idx);
abs_speeds = abs(speeds);
unique_speeds = unique(abs_speeds);

%% Plot BestSpeed, BestSpeed_fit, and fit_halfwidth

variablesToPlot = {'BestSpeed', 'Sparseness', 'BestAbsSpeed', 'SSI'}; %fit_width, abs_fit_width, BestSpeed_fit & BestAbsSpeed_fit (same as BestSpeed now and BestAbsSpeed now)
        
%Remove variables not in table (e.g. R2 might not be included)
variablesToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];
   
figure('units','normalized','outerposition',[0 0 1 1]); hold on
tiledlayout(length(Groups)*2,length(variablesToPlot))

for g = 1:length(Groups)
    isGroup = strcmp(Groups{g},GroupList);
    tempStim = StimAnalysis(isGroup,:);
    
    response_mat = tempStim.Response;
    abs_response_mat = tempStim.Abs_Response;
    
    %Normalize response_mats
    max_resp = max(response_mat,[],2);
    norm_response_mat = response_mat./max_resp;
    
    max_abs_resp = max(abs_response_mat,[],2);
    norm_abs_response_mat = abs_response_mat./max_abs_resp;
    
    for v = 1:length(variablesToPlot)
        vname = regexprep(variablesToPlot{v},'_',' ','emptymatch');
        
        %HISTOGRAM
        nexttile((g-1)*2*length(variablesToPlot) + v); hold on
        y = tempStim.(variablesToPlot{v});
        histogram(y)%,'facecolor','g','facealpha',0.2,'edgecolor','none');
        vline(mean(y, 'omitnan'), 'r')
        vline(mean(y, 'omitnan') - std(y))
        vline(mean(y, 'omitnan') + std(y))
        title(GroupsLabel{g})
        subtitle(strcat({'Mean '}, vname, {' = '}, num2str(mean(y, 'omitnan'))))
        xlabel(vname)
        ylabel('N cells')
        legend(['N = ' num2str(length(y))])
        
        %HEAT MAPS
        nexttile((g-1)*2*length(variablesToPlot) + length(variablesToPlot) + v); hold on
        if any(strcmp(variablesToPlot{v}, {'BestSpeed', 'BestSpeed_fit', 'Sparseness'}))
            respmat = norm_response_mat; %Plot all speeds from - to +
            xlab = 'oct/sec';
            xtick = speeds;
        else %SSI, BestAbsSpeed
            respmat = norm_abs_response_mat; %Plot absolute speeds
            xlab = 'abs oct/sec';
            xtick = unique_speeds;
        end
        [~, sortind] = sort(y, 'ascend');
        sortedmat = respmat(sortind,:);
        imagesc(sortedmat)
        xlim([0.5 size(sortedmat,2)+0.5])
        ylim([0.5 size(sortedmat,1)+0.5])
        set(gca, 'XTick', 1:size(sortedmat,2))
        set(gca, 'XTickLabel', num2str(xtick)) 
        h = colorbar;
        h.Title.String = 'AUC';
        caxis([0 1]); %Normalized responses should be between 0 to 1
        xlabel(xlab)
        ylabel('N cells')
        title(['Resp. sorted by ' vname])
    end
end

tak_suptitle(StimInfo.StimType)
drawnow

%% Plot MI, Slope, SSI, and Sparseness

variablesToPlot = {'MI', 'SMI', 'Slope', 'R2'};
        
%Remove variables not in table (e.g. R2 might not be included)
variablesToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];
   
figure('units','normalized','outerposition',[0 0 1 1]); hold on
tiledlayout(length(Groups)*2,length(variablesToPlot))

for g = 1:length(Groups)
    isGroup = strcmp(Groups{g},GroupList);
    tempStim = StimAnalysis(isGroup,:);

    response_mat = tempStim.Response;
    abs_response_mat = tempStim.Abs_Response;
    
    %Normalize response_mats
    max_resp = max(response_mat,[],2);
    norm_response_mat = response_mat./max_resp;
    
    max_abs_resp = max(abs_response_mat,[],2);
    norm_abs_response_mat = abs_response_mat./max_abs_resp;
    
    for v = 1:length(variablesToPlot)
        vname = regexprep(variablesToPlot{v},'_',' ','emptymatch');

        %HISTOGRAM
        nexttile((g-1)*2*length(variablesToPlot) + v); hold on
        y = tempStim.(variablesToPlot{v});
        histogram(y) %,'facecolor','g','facealpha',0.2,'edgecolor','none');
        vline(mean(y, 'omitnan'), 'r')
        vline(mean(y, 'omitnan') - std(y))
        vline(mean(y, 'omitnan') + std(y))
        title(GroupsLabel{g})
        subtitle(strcat({'Mean '}, vname, {' = '}, num2str(mean(y, 'omitnan'))))
        xlabel(vname)
        ylabel('N cells')
        
        %Also plot significant values
        figlegend = {['Total N = ' num2str(length(y))]};
        if strcmp(variablesToPlot{v}, 'MI') || strcmp(variablesToPlot{v}, 'Slope')
            
            p = tempStim.([variablesToPlot{v} '_p']);
            y_p = y(p < 0.05);
            if ~isempty(y_p)
                histogram(y_p,10,'facecolor','r')%,'facealpha',0.2,'edgecolor','none');
                figlegend{end + 1} = ['Sig shuffled N = ' num2str(length(y_p))];
            end
            
            p_bs = tempStim.([variablesToPlot{v} '_p_bs']);
            y_p_bs = y(p_bs < 0.05);
            if ~isempty(y_p_bs)
                histogram(y_p_bs,10,'facecolor','k')%,'facealpha',0.2,'edgecolor','none');
                figlegend{end + 1} = ['Sig bootstrap N = ' num2str(length(y_p_bs))];
            end
        end
        legend(figlegend)
        
        %HEAT MAPS
        nexttile((g-1)*2*length(variablesToPlot) + length(variablesToPlot) + v); hold on
        if any(strcmp(variablesToPlot{v}, {'MI'}))
            respmat = norm_response_mat; %Plot all speeds from - to +
            xlab = 'oct/sec';
            xtick = speeds;
        else %SMI, Slope
            respmat = norm_abs_response_mat; %Plot absolute speeds
            xlab = 'abs oct/sec';
            xtick = unique_speeds;
        end
        [~, sortind] = sort(y, 'ascend');
        sortedmat = respmat(sortind,:);
        imagesc(sortedmat)
        xlim([0.5 size(sortedmat,2)+0.5])
        ylim([0.5 size(sortedmat,1)+0.5])
        set(gca, 'XTick', 1:size(sortedmat,2))
        set(gca, 'XTickLabel', num2str(xtick)) 
        h = colorbar;
        h.Title.String = 'AUC';
        caxis([0 1]); %Normalized responses should be between 0 to 1
        xlabel(xlab)
        ylabel('N cells')
        title(['Resp. sorted by ' vname])
    end
end

tak_suptitle(StimInfo.StimType)
drawnow

%% Plot correlations

if Ops.plotCorrelations

    variablesToPlot = {'BestSpeed', 'Sparseness', 'BestAbsSpeed', 'SSI', 'MI', 'SMI', 'Slope', 'R2'}; %fit_width
        
    %Remove variables not in table (e.g. R2 might not be included)
    variablesToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];

    for g = 1:length(Groups)
        figure('units','normalized','outerposition',[0 0 1 1]); hold on
        tiledlayout(length(variablesToPlot),length(variablesToPlot))

        isGroup = strcmp(Groups{g},GroupList);
        tempStim = StimAnalysis(isGroup,:);

        for v = 1:length(variablesToPlot)
            for vv = 1:length(variablesToPlot)
                vname1 = regexprep(variablesToPlot{v},'_',' ','emptymatch');
                vname2 = regexprep(variablesToPlot{vv},'_',' ','emptymatch');

                x = tempStim.(variablesToPlot{v});
                y = tempStim.(variablesToPlot{vv});

                %Compute regression line (https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html)
                X = [ones(length(x),1) x];
                b = X\y;
                regline = X*b;
                Rsq = 1 - sum((y - regline).^2)/sum((y - mean(y, 'omitnan')).^2);

                nexttile; hold on
                scatter(x, y)
                plot(x, regline, 'r')
                xlabel(vname1)
                ylabel(vname2)
                title(['Rsq = ' num2str(round(Rsq,2))])
            end
        end
        tak_suptitle([StimInfo.StimType ' ' GroupsLabel{g}])
        drawnow
    end
end
