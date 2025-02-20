function plot_simple_ExtractedData_Reliability(ExtractedData, UserOps)
%Plot reliability data for ExtractedData created by simple_extract_data

%V1 archived 2/8/2023
%V2 = current version
%% Options

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

cb = 2; % 1:use color brewer; 2: use TAK lab colors to make predefined color schemes for graphs
if cb >0
    % choose colors here:
    colorsForGraphs = {'Reds', 'Blues'};
    % if cb = 1 some color options are:
    %{'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
    % 'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Spectral'};
    % but you can use cbrewer.mat for more options
    % if cb = 2, you are using TAK lab color schemes and
    % your options are {'Blues','Greens','Purples', 'Reds'};
end

%% Check for reliability data

if ~(ExtractedData.Ops.use_reliability)
    %No reliability data to plot
    return
end

if ~isfield(ExtractedData.Ops, 'reliability_type')
    ExtractedData.Ops.reliability_type = "Pearsons"; %Accommodate old ExtractedData files
end

%% Plot reliability data

Rfields = {'R1', 'R2'};

ReliabilityLabels = {' (All)', ' (IsRF)'};

Activity = {'all', 'activated', 'prolonged', 'suppressed'};

Loco = fieldnames(ExtractedData.CellData);

if UserOps.sortbyGCAMP
    Groups = unique(ExtractedData.Summary.(Loco{1}).GCaMP);
else
    Groups = unique(ExtractedData.Summary.(Loco{1}).Group);
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

for L = 1:length(Loco) %One figure per loco type
    
    figure; tiledlayout(length(Activity), length(Rfields))

    for a = 1:length(Activity)
        
        Ops = UserOps; %Copy and overwrite some things
        Ops.Loco = Loco{L};
        Ops.RF_Type = ''; %Include everything
        Ops.IsResponsive = 1; %Include everything
        Ops.SuppressOutput = 1;
            
        if strcmp(Activity{a}, 'all')
            Ops.ResponseType = ''; %Include everything
        else
            Ops.ResponseType = Activity{a}; 
        end
        
        SubsetData = simple_subset_ExtractedData(ExtractedData, Ops);
        Summary = SubsetData.Summary;

        %% Plot figures
        for r = 1:length(Rfields)
            
            nexttile; hold on

            for g = 1:length(Groups)
                groupInd = strcmp(Summary.Group, Groups{g});
                
                Y = Summary.(Rfields{r})(groupInd,:);
                X = zeros(size(Y)) + g;
                meanY = mean(Y, 'omitnan');

                scatter(X, Y, 50, 'jitter', 'on')
                line([g-0.25 g+0.25], [meanY meanY], 'Color', 'k', 'Linewidth', 2)
            end

            %Plot options
            set(gca, 'XTick', 1:length(Groups))
            set(gca, 'XTickLabel', GroupsLabel)
            xlim([0.5 length(Groups) + 0.5])
            ylim([0 1])
            if r == 1; ylabel([Activity{a} ' R']); end
            if a == 1; title([regexprep(Rfields{r},'_',' ','emptymatch') ReliabilityLabels{r}]); end
        end
        tak_suptitle(Loco{L})
    end
end  
    