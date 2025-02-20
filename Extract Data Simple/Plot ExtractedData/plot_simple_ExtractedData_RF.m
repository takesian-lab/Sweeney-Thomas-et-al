function plot_simple_ExtractedData_RF(ExtractedData, UserOps)
% Plot output of RF stim analysis for ExtractedData

%Version control:
%V1: archived 2/1/2023
%V2: current version

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%Return if not RF data
if ExtractedData.StimInfo.StimProtocol ~= 2
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

disp('Plotting RF data')

StimData = simple_subset_ExtractedData(ExtractedData, UserOps);

Ops = StimData.SubsetOps;
StimInfo = StimData.StimInfo;
Summary = StimData.Summary;
StimAnalysis = StimData.StimAnalysis;
GroupList = Summary.Group;
Groups = unique(GroupList);

if isempty(Summary)
    disp('No responsive cells found to plot for RF')
    return
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

%% Plot BF and CF measures

variablesToPlot = {'BF', 'BF_I', 'BWBF_I', 'CF', 'CF_I', 'BW20', 'dPrime', 'Sparseness', 'R2'};
responseToPlot = {'BF_I', 'BF_I', 'BF_I', 'CF_I', 'CF_I', 'CF_I+20', 'BF_I', 'BF_I', 'BF_I'}; %For heatmap (All intensity rows)
        
%Remove variables not in table (e.g. R2 might not be included)
variablesToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];
responseToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];       
        
figure('units','normalized','outerposition',[0 0 1 1]); hold on
tiledlayout(length(Groups)*2,length(variablesToPlot))

for g = 1:length(Groups)
    isGroup = strcmp(Groups{g},GroupList);
    tempStim = StimAnalysis(isGroup,:);
    
    response_mat = tempStim.Response;
    
    %Normalize response_mats
    max_resp = max(response_mat,[],2);
    norm_response_mat = response_mat./max_resp;
    
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
        
        %For heatmat, get values represented in responseToPlot
        %i.e. instead of plotting all 64 responses in RF, just plot the intensity row or frequency column specified
        [respmatind, respmat] = deal(nan(length(y),length(StimInfo.V1)));
        if strcmp(responseToPlot{v}, 'BF_I') || strcmp(responseToPlot{v}, 'CF_I')
            ints = tempStim.(responseToPlot{v});
            %Fill respmat with responses for that intensity row
            for r = 1:size(ints,1)
                respmatind(r,:) = find(StimInfo.Combined(:,2) == ints(r));
                respmat(r,:) = norm_response_mat(r,respmatind(r,:));
            end
        elseif strcmp(responseToPlot{v}, 'CF_I+20')
            ints = tempStim.CF_I + 20;
            ints(ints > max(StimInfo.V2)) = nan;
            %Fill respmat with responses for that intensity row
            for r = 1:size(ints,1)
                if ~isnan(ints(r))
                    respmatind(r,:) = find(StimInfo.Combined(:,2) == ints(r));
                    respmat(r,:) = norm_response_mat(r,respmatind(r,:));
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
        xlabel('Frequency (kHz)')
        ylabel('N cells')
        title([responseToPlot{v} ' sorted by ' vname])
    end
end

tak_suptitle(StimInfo.StimType)
drawnow

%% Plot everything else

variablesToPlot = {'BestInt', 'BWInt', 'ISI', 'dPrime', 'Sparseness', 'R2'};
responseToPlot = {'BF', 'BF', 'BF', 'BF', 'BF', 'BF'}; %For heatmap (All BF freq column)
        
%Remove variables not in table (e.g. R2 might not be included)
variablesToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];
responseToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];       

figure('units','normalized','outerposition',[0 0 1 1]); hold on
tiledlayout(length(Groups)*2,length(variablesToPlot))

for g = 1:length(Groups)
    isGroup = strcmp(Groups{g},GroupList);
    tempStim = StimAnalysis(isGroup,:);
    
    response_mat = tempStim.Response;
    
    %Normalize response_mats
    max_resp = max(response_mat,[],2);
    norm_response_mat = response_mat./max_resp;
    
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
        
        %For heatmat, get values represented in responseToPlot
        %i.e. instead of plotting all 64 responses in RF, just plot the intensity row or frequency column specified
        [respmatind, respmat] = deal(nan(length(y),length(StimInfo.V2)));
        if strcmp(responseToPlot{v}, 'BF')
            BF = tempStim.BF;
            %Fill respmat with responses for that intensity row
            for r = 1:size(BF,1)
                respmatind(r,:) = find(StimInfo.Combined(:,1) == BF(r));
                respmat(r,:) = norm_response_mat(r,respmatind(r,:));
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
        set(gca, 'XTickLabel', StimInfo.V2) 
        h = colorbar;
        h.Title.String = 'AUC';
        caxis([0 1]); %Normalized responses should be between 0 to 1
        xlabel('Intensities (dB)')
        ylabel('N cells')
        title([responseToPlot{v} ' sorted by ' vname])
    end
end

tak_suptitle(StimInfo.StimType)
drawnow

%% Plot correlations

if Ops.plotCorrelations
    variablesToPlot = {'BF', 'BF_I', 'BWBF_I', 'CF', 'CF_I', 'BW20', 'BestInt', 'BWInt', 'ISI', 'dPrime', 'Sparseness'};
        
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
                Rsq = 1 - sum((y - regline).^2)/sum((y - mean(y)).^2);

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

