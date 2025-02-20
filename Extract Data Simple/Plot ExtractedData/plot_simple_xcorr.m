function [fig1, fig2, fig3, fig4] = plot_simple_xcorr(CellDataByStim, StimInfo, Ops)

%% Default options

if nargin < 3
    Ops = struct;
    Ops.smTime = 0;
end

%% Get variables from StimInfo
V1 = StimInfo.V1;
V2 = StimInfo.V2;
CombinedStim = StimInfo.Combined;
nPlots = length(V1)*length(V2);

StimTraces = CellDataByStim.StimTraces;
blanks = CellDataByStim.BlankTraces;
    
%Optional smooth traces with a moving window
if Ops.smTime > 0
    smoothval = ceil(Ops.smTime/(1/Ops.fs)); %window size
    for a = 1:length(StimTraces)
        StimTraces{a} = single(smoothdata(StimTraces{a}, 2,'movmean',smoothval));
    end
    blanks = single(smoothdata(blanks, 2,'movmean',smoothval));
end

%% plot traces

fig1 = figure;
h1 = tight_subplot(length(V2),length(V1), 0.04, [0.06 0.1], 0.06);

for a = 1:nPlots
    axes(h1(a)); hold on
    
    r_trials = StimTraces{a};
    
    plot(r_trials')
    plot(mean(r_trials,1),'k','LineWidth', 1.5)
    vline(StimInfo.nBaselineFrames, 'r')
    
    if a <= length(V1) %Top row
        title([num2str(CombinedStim(a,1)) ' ' StimInfo.Units{1}]);
    end
    if a > nPlots - length(V1) %Bottom row
        xlabel('Frames');
    end
    
    if ismember(a,[1:length(V1):nPlots]) %First column
        ylabel([num2str(CombinedStim(a,2)) ' ' StimInfo.Units{2}]);
    end
end

%% plot blanks

if ~isempty(blanks)
    fig2 = figure; hold on
    plot(blanks')
    plot(mean(blanks,1),'k','LineWidth', 2)
    vline(StimInfo.nBaselineFrames, 'r')
    title('Blank trials')
    xlabel('Frames');
    ylabel('Z-Score');
else
    fig2 = [];
end

%% plot histogram of xcorr blanks compared to xcorr metric

variables = {'cc_max', 'cc_zero'};

for v = 1:2
    if v == 1
        fig3 = figure;
    elseif v == 2
        fig4 = figure;
    end
        
    h3 = tight_subplot(length(V2),length(V1), 0.04, [0.06 0.1], 0.06);

    for a = 1:nPlots
        axes(h3(a)); hold on

        mode             = CellDataByStim.XcorrData.Mode(a);
        R                = CellDataByStim.XcorrData.(strcat(variables{v}))(a);
        z_test           = CellDataByStim.XcorrData.(strcat(variables{v},'_z'))(a);
        z_R              = CellDataByStim.XcorrData.(strcat(variables{v},'_r'))(a);
        reliability_dist = CellDataByStim.XcorrData.(strcat(variables{v},'_dist')){a};

        %X and Y labels
        if a <= length(V1) %Top row
            subtitle(strcat('.'), 'FontSize',5); %Placeholder for formatting
            title([num2str(CombinedStim(a,1)) ' ' StimInfo.Units{1}]);
        end
        if a > nPlots - length(V1) %Bottom row
            xlabel('Frames');
        end

        if ismember(a,[1:length(V1):nPlots]) %First column
            ylabel([num2str(CombinedStim(a,2)) ' ' StimInfo.Units{2}]);
        end

        %Skip stim that reliability was not computed for
        if isnan(z_test)
            continue;
        end

        %Set histogram and text colors
        if strcmp(mode,'Shifted')
            colour = 'b';
            if z_test
                textcolour = ' {\color{blue}Z=} ';
            else
                textcolour = ' {\color{black}Z=} ';
            end
        elseif strcmp(mode,'Blanks')
            colour = 'g';
            if z_test
                textcolour = ' {\color{green}Z=} ';
            else
                textcolour = ' {\color{black}Z=} ';
            end
        end

        %Plot histogram
        min_g = min(min(reliability_dist), R)-0.1;
        max_g = max(max(reliability_dist), R)+0.1;

        h = histogram(reliability_dist,min_g:.01:max_g,'facecolor',colour,'facealpha',0.2,'edgecolor','none');
        [maxcount, whichbin] = max(h.Values); hold on;

        %Plot line for R value
        plot([R R],[0 maxcount], '--r'); hold on;

        %Test result
        subtitle(strcat(regexprep(variables{v},'_',' ','emptymatch'), '= ', sprintf('%.2f',R), textcolour, sprintf('%.2f', z_R)), 'FontSize',8);
    end
end
