function [fig1] = simple_plot_PeakData(StimInfo, PeakData, varargin)
%Plot RF using various measures

%% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% Use custom RF values instead of PeakData.HasPeak
addOptional(p, 'RF', []); 

% Change opacity of IsRF mask (1 = no mask)
addOptional(p, 'PercentOpaque', 0.3); 

% Use blue to red colormap instead of MATLAB colors
addOptional(p, 'bluewhitered', 1); 

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% Setup

RF_Map = StimInfo.RF_Map;
nanmat = nan(max(RF_Map));
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

variables = {'Combined 3fr', 'Final RF',...
    'HasPeak', 'Peak', 'Peak_3fr', 'Peak_win', 'Peak_avg', 'Peak_Onset', 'Peak_Latency', 'Peak_Width', 'Peak_AUC'...
    'HasTrough', 'Trough', 'Trough_3fr', 'Trough_win', 'Trough_avg', 'Trough_Onset', 'Trough_Latency', 'Trough_Width', 'Trough_AUC'};

colorbartitles = {'Z', '0 or 1',...
    '0 or 1', 'Z', 'Z', 'Z', 'Z', 'sec', 'sec', 'sec', 'Z'...
    '0 or 1', 'Z', 'Z', 'Z', 'Z', 'sec', 'sec', 'sec', 'Z'};

q = [1, 6,...
    11 2 3 4 5 8 9 10 7, ...
    16 12 13 14 15 18 19 20 17]; %plot order
    
%% Settings to mask non-responsive stim

IsRF = nanmat;
if isempty(ops.RF)
    IsRF(RF_ind) = PeakData.IsResponsive;
else
    IsRF(RF_ind) = ops.RF;
end

RFmask = ones(size(nanmat))*ops.PercentOpaque;
RFmask(IsRF == 1) = 1;

%% Plot

fig1 = figure('units','normalized','outerposition',[0 0 1 1]); hold on

for v = 1:length(variables)
    
    vmat = nanmat;
    if strcmp(variables{v}, 'Final RF')
        vmat = IsRF;
    elseif strcmp(variables{v}, 'Combined 3fr')
        %Generate RF from peaks and troughs
        temp_RF = nan(length(PeakData.Peak),1);
        for i = 1:length(PeakData.Peak)
            response_type = PeakData.ResponseType(i);
            if strcmp(response_type, 'suppressed')
                temp_RF(i) = PeakData.Trough_3fr(i);
            else
                temp_RF(i) = PeakData.Peak_3fr(i);
            end
        end
        vmat(RF_ind) = temp_RF;
    else
        vmat(RF_ind) = PeakData.(variables{v});
    end

    subplot(4,5,q(v))
    h = imagesc(vmat);
    set(h, 'AlphaData', RFmask)
    title(regexprep(variables{v},'_',' ','emptymatch')) %Replace underscore with a space in plot title
    set(gca,'YTick',1:length(StimInfo.V2))
    set(gca,'YTickLabel',num2str(StimInfo.V2))
    set(gca,'XTick',1:length(StimInfo.V1))
    set(gca,'XTickLabel',num2str(StimInfo.V1))
    
    if ismember(v, [1 2 3 12])
        ylabel([StimInfo.Parameters{2} ' (' StimInfo.Units{2} ')'])        
    end
    
    if ismember(v, [12, 17:20])
        xlabel([StimInfo.Parameters{1} ' (' StimInfo.Units{1} ')'])
    end
    
    if ops.bluewhitered
        colormap(gca, bluewhitered(256)) %TODO: replace with bluewhitered_TAKlab
    end
    c = colorbar;
    c.Title.String = colorbartitles{v};
    
    if strcmp(colorbartitles{v}, '0 or 1')
        caxis([0 1])
    end
end