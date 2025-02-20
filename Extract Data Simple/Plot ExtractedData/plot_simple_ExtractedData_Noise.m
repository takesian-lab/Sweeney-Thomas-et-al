function plot_simple_ExtractedData_Noise(ExtractedData, UserOps)
% Plot output of Noise stim analysis for ExtractedData
% For noise with varying intensities only

%Version control:
%V1: current version

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%Return if not Noise data
if ExtractedData.StimInfo.StimProtocol ~= 10
    return;
end

%We are expecting data that only varies in the intensity dimension
if ~strcmp(ExtractedData.StimInfo.Parameters{1}, 'Intensity') || (numel(ExtractedData.StimInfo.V2) > 1)
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

disp('Plotting Noise data')

StimData = simple_subset_ExtractedData(ExtractedData, UserOps);

Ops = StimData.SubsetOps;
StimInfo = StimData.StimInfo;
Summary = StimData.Summary;
StimAnalysis = StimData.StimAnalysis;
GroupList = Summary.Group;
Groups = unique(GroupList);

if isempty(Summary)
    disp('No responsive cells found to plot for Noise')
    return
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

%% Plot Best Intensity rasters

variablesToPlot = {'BestInt', 'Raw_BestInt', 'ISI', 'Raw_ISI', 'FRA_area'};
        
%Remove variables not in table (e.g. R2 might not be included)
variablesToPlot(~ismember(variablesToPlot,StimAnalysis.Properties.VariableNames)) = [];
   
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
        
        [~, sortind] = sort(y, 'ascend');
        sortedmat = norm_response_mat(sortind,:);
        imagesc(sortedmat)
        xlim([0.5 size(sortedmat,2)+0.5])
        ylim([0.5 size(sortedmat,1)+0.5])
        set(gca, 'XTick', 1:size(sortedmat,2))
        set(gca, 'XTickLabel', StimInfo.V1) 
        h = colorbar;
        h.Title.String = 'AUC';
        caxis([0 1]); %Normalized responses should be between 0 to 1
        xlabel('Intensity (dB)')
        ylabel('N cells')
        title(['Sorted by ' vname])
    end
end

tak_suptitle(StimInfo.StimType)
drawnow

%% Plot correlations

if Ops.plotCorrelations
    variablesToPlot = {'BestInt', 'Raw_BestInt', 'ISI', 'Raw_ISI', 'FRA_area'};
        
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

%% Summary figure
%TODO: Plot by Group

intsToPlot = [40, 60, 80];
strIntsToPlot = {'40', '60', '80'};
BestInt = StimAnalysis.Raw_BestInt;

figure;

%Pie charts
subplot(3,3,1);
tab = nan(size(intsToPlot));
for i = 1:length(intsToPlot)
    tab(i) = sum(BestInt == intsToPlot(i));
end
pie(tab, strIntsToPlot)
title('Best Int')

%ISI scatter
subplot(3,3,2:3); hold on
for i = 1:length(intsToPlot)
    x = zeros(1,tab(i))+i;
    y = StimAnalysis.ISI(BestInt == intsToPlot(i));
    scatter(x,y,50,'jitter','on');
    meany = mean(y,'omitnan');
    semy = std(y,'omitnan')./(sqrt(length(y)-1));
    line([i-0.25 i+0.25], [meany meany], 'color','k','linewidth',2)
    line([i i], [meany-semy meany+semy], 'color','k','linewidth',2)
end
title('ISI vs. Best Int')
ylabel('ISI')
set(gca,'XTick',1:length(intsToPlot))
set(gca,'XTickLabel', strIntsToPlot)

%Plot dose response curves for cells with best intensity at each dB
for i = 1:length(intsToPlot)
    subplot(3,3,3+i); hold on
    tempData = StimAnalysis(BestInt == intsToPlot(i),:);
    tempRespMat = tempData.Response;
    
    respMat = [];
    for ii = 1:length(intsToPlot)
        respMat = [respMat, tempRespMat(:,StimInfo.V1' == intsToPlot(ii))];
    end
    
    %Normalize
    min_respMat = min(respMat,[],2);
    sign_min = sign(min_respMat);
    respMat(sign_min < 0,:) = respMat(sign_min < 0,:) + abs(min_respMat(sign_min < 0));
    
    max_respMat = max(respMat,[],2);
    respMat = respMat./max_respMat;
    
    mean_mat = mean(respMat,1,'omitnan');
    sem_mat = std(respMat,[],1,'omitnan')./(sqrt(height(respMat)-1));
    
    for ii = 1:length(intsToPlot)
        %scatter(zeros(1,height(respMat))+ii, respMat(:,ii),50,'jitter','on')
    end
    errorbar(1:length(intsToPlot),mean_mat,sem_mat,'Color','k','Linewidth',2)
    
    ylim([0 1.1])
    xlim([0.5 length(intsToPlot)+0.5])
    set(gca,'XTick',1:length(intsToPlot))
    set(gca,'XTickLabel',strIntsToPlot)
    ylabel('Norm. AUC')
    title(['Cells with best int ' num2str(intsToPlot(i)) 'dB'])
    subtitle(['N = ' num2str(height(respMat))])
end


%Plot dose response curves again without normalizing
for i = 1:length(intsToPlot)
    subplot(3,3,6+i); hold on
    tempData = StimAnalysis(BestInt == intsToPlot(i),:);
    tempRespMat = tempData.Response;
    
    respMat = [];
    for ii = 1:length(intsToPlot)
        respMat = [respMat, tempRespMat(:,StimInfo.V1' == intsToPlot(ii))];
    end
    
    mean_mat = mean(respMat,1,'omitnan');
    sem_mat = std(respMat,[],1,'omitnan')./(sqrt(height(respMat)-1));
    
    for ii = 1:length(intsToPlot)
        %scatter(zeros(1,height(respMat))+ii, respMat(:,ii),50,'jitter','on')
    end
    errorbar(1:length(intsToPlot),mean_mat,sem_mat,'Color','k','Linewidth',2)
    
    ylim([0 60])
    xlim([0.5 length(intsToPlot)+0.5])
    set(gca,'XTick',1:length(intsToPlot))
    set(gca,'XTickLabel',strIntsToPlot)
    ylabel('Raw AUC')
    xlabel('dB')
end

sgtitle('Intensity Tuning')
