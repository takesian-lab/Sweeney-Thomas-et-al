function plot_simple_ExtractedData_PeakData(ExtractedData, UserOps)
% Plot output of PeakData and sorted final traces

%Version control:
%V1: archived 2/1/2023
%V2: current version

%% Get data from ExtractedData

%DO NOT CHANGE THESE: USE USEROPS
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
    UserOps.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress
end

disp('Plotting peak data')

%Change what data is plotted in this function (e.g. group, gcamp, red cell, etc)
SubsetData = simple_subset_ExtractedData(ExtractedData, UserOps);

Ops = SubsetData.SubsetOps;
StimInfo = SubsetData.StimInfo;
Summary = SubsetData.Summary;
FinalTraces = SubsetData.FinalTraces;
GroupList = Summary.Group;
Groups = unique(GroupList);
stimname = 'Average traces';
peaktype = SubsetData.peaktype;

%If no ResponseType specified, default to sorting by peaks instead of troughs
if strcmp(peaktype,'None')
    peaktype = 'Peak';
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

%X labels for figures
xlabel_increment = 0.5; %seconds
fs = StimInfo.fs;
nFrames = StimInfo.nFrames;
increment_in_frames = xlabel_increment*fs;
xtick = 0:increment_in_frames:nFrames;
xticklabel = 0:xlabel_increment:((1/fs)*nFrames);
soundtime = StimInfo.baseline*fs;

if isempty(Summary)
    disp('No responsive cells found to plot for peak data')
    return
end
%% Plot Graph

variablesToPlot = {'', '_3fr', '_win', '_avg', '_Onset', '_Latency', '_Width', '_AUC', 'R2'};
sortdirection = {'ascend', 'ascend', 'ascend', 'ascend', 'descend', 'descend', 'ascend', 'ascend', 'ascend'};

%Concatenate Peak or Trough to the beginning of each variables, except R2
for v = 1:length(variablesToPlot)
    if strcmp(variablesToPlot{v}, 'R2')
        continue;
    end
    variablesToPlot{v} = strcat(peaktype, variablesToPlot{v});
end
      
%Remove variables not in table (e.g. R2 might not be included)
variablesToPlot(~ismember(variablesToPlot,Summary.Properties.VariableNames)) = [];
        
figure('units','normalized','outerposition',[0 0 1 1]); hold on

tiledlayout(length(Groups)*2,length(variablesToPlot))

for g = 1:length(Groups)
    isGroup = strcmp(Groups{g},GroupList);
    tempStim = Summary(isGroup,:);
    
    respmat = FinalTraces(isGroup,:);
    
    if Ops.smooth_rasters
        for i = 1:size(respmat,1)
            respmat(i,:) = smooth(respmat(i,:),5);
        end
    end
    
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
        
        [~, sortind] = sort(y, sortdirection{v});
        sortedmat = respmat(sortind,:);
        imagesc(sortedmat)
        xlim([0.5 size(sortedmat,2)])
        ylim([0.5 size(sortedmat,1)])
        vline(soundtime, 'k')
        set(gca, 'XTick', xtick)
        set(gca, 'XTickLabel', xticklabel) 
        h = colorbar;
        h.Title.String = 'Z';
        caxis([0 2]); %Might have to adjust this
        xlabel('Time (s)')
        ylabel('N cells')
        title(['Sorted by ' vname])
    end
end

tak_suptitle(stimname)
drawnow

%% Plot correlations

if Ops.plotCorrelations
    variablesToPlot = {'', '_3fr', '_win', '_avg', '_Onset', '_Latency', '_Width', '_AUC', 'R2'};

    %Concatenate Peak or Trough to the beginning of each variables, except R2
    for v = 1:length(variablesToPlot)
        if strcmp(variablesToPlot{v}, 'R2')
            continue;
        end
        variablesToPlot{v} = strcat(peaktype, variablesToPlot{v});
    end        
    
    %Remove variables not in table (e.g. R2 might not be included)
    variablesToPlot(~ismember(variablesToPlot,Summary.Properties.VariableNames)) = [];
   
    for g = 1:length(Groups)
        figure('units','normalized','outerposition',[0 0 1 1]); hold on
        tiledlayout(length(variablesToPlot),length(variablesToPlot))

        isGroup = strcmp(Groups{g},GroupList);
        tempStim = Summary(isGroup,:);

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
        tak_suptitle([stimname ' ' GroupsLabel{g}])
        drawnow
    end
end

