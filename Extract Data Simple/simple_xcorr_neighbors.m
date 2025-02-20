function [New_RF, NeighborData, fig1, fig2] = simple_xcorr_neighbors(CellDataByStim, StimInfo, plot_figures, Ops)
% This function maximizes detection of a response within a receptive field 
% by measuring response reliability (using xcorr) compared to neighboring stimuli 

% Argument(s): 
% CellDataByStim - from simple_extract_data
% StimInfo
% plot_figures
% Ops.reliability_shuffles - user-defined (in setup) # of shuffles, suggested 1000 
% Ops.neighbor_p_value - alpha level
% Ops.smTime - moving window (s) for smoothing raw traces. Set to [] for no smoothing

% Returns NeighborData with:
% - V1, V2: record stim info
% - Visited: whether reliability was computed for that stim condition
% - Best neighbor: If condition was visited more than once, only results
% for its best neighbor (lowest p value) will be saved
% - Ztest - result of the test
% - P - p-value of the test
% - R -

% Argument(s): 


% Returns:

%% Default options

if nargin < 4
    Ops = struct;
    Ops.reliability_shuffles = 1000;
    Ops.neighbor_p_value = 0.05;
    Ops.adjust_neighbor_p_value = 1;
    Ops.smTime = 0;
end

r_shuffles = Ops.reliability_shuffles;
maxlag_in_s = 0.1;

%% Setup

plot_troubleshooting_figures = 0;

%Variables will return as NaN if neighbors are not tested
RF = CellDataByStim.RF;
New_RF = RF;
NeighborData = table;
NeighborData.V1 = StimInfo.Combined(:,1);
NeighborData.V2 = StimInfo.Combined(:,2);
NeighborData.Visited = false(size(CellDataByStim.RF));
NeighborData.BestNeighbor = nan(size(CellDataByStim.RF));
NeighborData.Ztest = nan(size(CellDataByStim.RF));
NeighborData.P = nan(size(CellDataByStim.RF));
NeighborData.R = nan(size(CellDataByStim.RF));
            
%Get info from StimInfo
freqs = StimInfo.V1;
ints = StimInfo.V2;
size_field = [length(ints), length(freqs)];
fig1 = []; fig2 = [];

%Do not continue if there are no responsive stim conditions or no neighbors
if ~any(CellDataByStim.RF) || isequal(size_field, [1 1])
    return
end

%Optional smooth traces with a moving window
StimTraces = CellDataByStim.StimTraces;
if Ops.smTime > 0
    smoothval = ceil(Ops.smTime/(1/StimInfo.fs)); %window size
    for a = 1:length(StimTraces)
        StimTraces{a} = single(smoothdata(StimTraces{a}, 2,'movmean',smoothval));
    end
end

fs = StimInfo.fs;
soundtime = StimInfo.baseline*fs;
maxlag = floor(fs*maxlag_in_s);

%% Pearson Correlation with adjacent trials

for v = 1:length(StimInfo.Combined)

    %Skip stim condition if not responsive/reliable
    if ~RF(v)
        continue
    end
    
    %Find adjacent points
    RF_ind = StimInfo.RF_Map(v,:);
    adj_ind = get_adjacent_points(RF_ind(1), RF_ind(2), size_field, 8, 'list', plot_troubleshooting_figures);

    %Skip neighbors that have no trials (e.g. because of loco), were already deemed Reliable, are not Responsive, or were already found to be a Reliable neighbor
    skip_neighbor = (cellfun(@isempty,StimTraces(adj_ind)) + RF(adj_ind) + ~CellDataByStim.PeakData.IsResponsive(adj_ind) + New_RF(adj_ind)) > 0;
    adj_ind(skip_neighbor) = [];
    
    if isempty(adj_ind)
        continue
    end
    
    %If valid neighbors are found, prepare to compute reliability
    trials = StimTraces{v};  
    
    %Remove baseline
    trials = trials(:,soundtime+1:end);
    
    N1 = size(trials,1); % number of trials for stim
    N2 = size(trials,2); % number of frames/trial
    %N3 = factorial(N1)/(2*factorial(N1-2)); % number of pairs (N!/(2*(N-2))!)
    
    %Replicate r_trials number of times needed for r_shuffles
    repeated_trials = repmat(trials,[r_shuffles,1]);
    
    %Shift the replicated trials by random frame
    rand_shifting = randsample(1:N2,N1*r_shuffles,'true')';  % random spots within traces to start shift (have to use replacement to do this outside of loop)
    shifted_trials = cell2mat(arrayfun(@(x) circshift(repeated_trials(x,:),[1 rand_shifting(x)]),(1:numel(rand_shifting))','un',0));

    % Measure reliability of stimulus responses compared to neighbors
    for a = 1:size(adj_ind) % loop through neighbors
        aa = adj_ind(a);
        neighbor_trials = StimTraces{aa};
        neighbor_trials = neighbor_trials(:,soundtime+1:end);
        N4 = size(neighbor_trials,1); % number of trials for each neighbor
        
        %Concatenate trials to perform xcorr
        all_trials = [trials; neighbor_trials];
        
        %Figure out how many unique pairs of trials there are
        pairs = [];
        count = 1;
        for n1 = 1:size(all_trials,1)
            for n2 = 1:size(all_trials,1)
                pairs(count,1:2) = [n1, n2];
                count = count + 1;
            end
        end
        
        %We only want pairs that include a trial from each set
        set_rows = (pairs(:,1) <= N1) + (pairs(:,2) >= N1) == 2;
        set_pairs = pairs(set_rows, :);
        nPairs = size(set_pairs,1);

        %Compute xcorr
        [cc_mat, lag] = xcorr(all_trials', maxlag, 'coeff');
        
        %Transpose cc_mats and only keep set pairs
        cc_mat = cc_mat(:,set_rows)';

        cc_max = mean(max(cc_mat,[],2),1);
        cc_zero = mean(cc_mat(:,lag == 0),1);

        % Measure reliability of shifted traces as control and determine Z-score for neighbor

        %Replicate neighor_trials number of times needed for r_shuffles
        repeated_neighbor_trials = repmat(neighbor_trials,[r_shuffles,1]);
        
        %Shift the replicated trials of neighbor by random frame
        rand_shifting = randsample(1:N2,N4*r_shuffles,'true')';  % random spots within traces to start shift (have to use replacement to do this outside of loop)
        shifted_trials_neighbor = cell2mat(arrayfun(@(x) circshift(repeated_neighbor_trials(x,:),[1 rand_shifting(x)]),(1:numel(rand_shifting))','un',0));

        %Size of shifted trials
        N5 = size(shifted_trials,1);
        N6 = size(shifted_trials_neighbor,1);

        %Pull out comparisons between N1 shifted trials from neighbors to compute the mean Pearsons & do this r_shuffles number of times
        max_shifted = nan(r_shuffles,1);
        zero_shifted = nan(r_shuffles,1);

        index = 1;
        index_neighbors = 1;
        for u = 1:r_shuffles
            shuffle_index = index:index+N1-1;
            shuffle_index_neighbors = index_neighbors:index_neighbors+N4-1;

            temp_trials = shifted_trials(shuffle_index,:);
            temp_neighbor_trials = shifted_trials_neighbor(shuffle_index_neighbors,:);
            temp_all_trials = [temp_trials; temp_neighbor_trials];
            
            %Compute xcorr
            [temp_cc_mat, temp_lag] = xcorr(temp_all_trials', maxlag, 'coeff');

            %Transpose cc_mats and only keep set pairs
            temp_cc_mat = temp_cc_mat(:,set_rows)';

            max_shifted(u) = mean(max(temp_cc_mat,[],2),1);
            zero_shifted(u) = mean(temp_cc_mat(:,temp_lag == 0),1);

            index = index+N1;
            index_neighbors = index_neighbors+N4;
        end

        % perform z test to determine if reliability metric R is from same
        % distribution as reliability measures of shifted trials,
        % comparing one response to a stimulus to a that of a neighbor
        if Ops.adjust_neighbor_p_value
            p_value = Ops.neighbor_p_value./length(adj_ind);
        else
            p_value = Ops.neighbor_p_value;
        end
        
        R = cc_zero; %cc_max
        reliability_shifted = zero_shifted; %max_shifted
        
        %Here, do a two-tailed test because sometimes a stimulus is NEGATIVELY correlated with its neighbor (-R) and that helps us pull out suppressed responses
        [z_test_shift, z_P_shift, z_CI_shift, z_R_shift] = ztest(R, mean(reliability_shifted), std(reliability_shifted),"Alpha",p_value);

        % if any neighbor shows correlated activity, flag it 'reliable'
        if z_test_shift == 1
            New_RF(aa) = 1; %Update value in New_RF to save
        end
        
        % save correlation data
        if ~NeighborData.Visited(aa)
            % if neighbor hasn't been visited yet, store correlation values
            NeighborData.Visited(aa) = 1;
            NeighborData.BestNeighbor(aa) = v;
            NeighborData.Ztest(aa) = z_test_shift;
            NeighborData.P(aa) = z_P_shift;
            NeighborData.R(aa) = z_R_shift;
        else
            % if neighbor has been visited, compare to current correlation values and store whichever is better
            [~,min_ind] = min([z_P_shift, NeighborData.P(aa)]); % find which p value is lower
            if min_ind == 1 % new neighbor is better
                NeighborData.BestNeighbor(aa) = v;
                NeighborData.Ztest(aa) = z_test_shift;
                NeighborData.P(aa) = z_P_shift;
                NeighborData.R(aa) = z_R_shift;
            end
        end

        % plot histogram of shifted reliability compared to relibiality metric
        if plot_troubleshooting_figures
            figure; hold on
            min_g = min(min(reliability_shifted), R)-0.1;
            max_g = max(max(reliability_shifted), R)+0.1;
            h = histogram(reliability_shifted,min_g:.01:max_g,'facecolor','g','facealpha',0.2,'edgecolor','none');
            [maxcount, whichbin] = max(h.Values);
            plot([R R],[0 maxcount], '--r');
            if z_P_shift < p_value
                colour = 'green';
            else
                colour = 'black';
            end
            subtitle(strcat('R= ', sprintf('%.2f',R), ...
                ' {\color{', colour, '}Z=} ', sprintf('%.2f', z_R_shift)), 'FontSize',10);
        end
    end  
end

if plot_figures

    %Plot RF maps with old and new IsRF
    RF_Map = StimInfo.RF_Map;
    nanmat = nan(max(RF_Map));
    RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));
    [old_map, new_map] = deal(nanmat);
    old_map(RF_ind) = CellDataByStim.RF;
    new_map(RF_ind) = New_RF;
    new_map(new_map + old_map == 1) = 0.5;
    
    fig1 = figure;
    subplot(2,1,1)
    imagesc(old_map)
    xlabel([StimInfo.Parameters{1} ' (' StimInfo.Units{1} ')'])
    ylabel([StimInfo.Parameters{2} ' (' StimInfo.Units{2} ')'])
    set(gca,'XTick',1:length(StimInfo.V1))
    set(gca,'YTick',1:length(StimInfo.V2))
    set(gca,'XTickLabel',num2str(StimInfo.V1))
    set(gca,'YTickLabel',num2str(StimInfo.V2))
    title('Old RF')
    
    subplot(2,1,2)
    imagesc(new_map)
    xlabel([StimInfo.Parameters{1} ' (' StimInfo.Units{1} ')'])
    ylabel([StimInfo.Parameters{2} ' (' StimInfo.Units{2} ')'])
    set(gca,'XTick',1:length(StimInfo.V1))
    set(gca,'YTick',1:length(StimInfo.V2))
    set(gca,'XTickLabel',num2str(StimInfo.V1))
    set(gca,'YTickLabel',num2str(StimInfo.V2))
    title('New RF with neighbors')
    
    %If you want to look at peak data, you can run the following functions
%   simple_plot_PeakData(StimInfo, CellDataByStim.PeakData, 'RF', RF); 
%   simple_plot_PeakData(StimInfo, CellDataByStim.PeakData, 'RF', New_RF);

    fig2 = figure('units','normalized','outerposition',[0 0 1 1]); hold on
    stimToPlot = find(CellDataByStim.RF);
    tiledlayout(length(stimToPlot),9) %9 = 8 max neighbors + current stimToPlot
    for vv = 1:length(stimToPlot)
        v = stimToPlot(vv);
        %Find adjacent points
        RF_ind = StimInfo.RF_Map(v,:);
        adj_ind = get_adjacent_points(RF_ind(1), RF_ind(2), size_field, 8, 'list', 0);
        
        nexttile; hold on
        traces = StimTraces{v};
        meantrace = mean(traces,1);
        plot(traces','Linewidth',0.5)
        plot(meantrace,'k','Linewidth',2)
        ylabel(['Stim index ' num2str(v)])
        if vv == 1; title('Responsive Condition'); end
        subtitle([num2str(round(NeighborData.V1(v),1)), 'kHz ', num2str(round(NeighborData.V2(v),1)), 'dB'])
            
        for a = 1:8
            if a > length(adj_ind)
                nexttile;
                continue
            end
            aa = adj_ind(a);
            nexttile; hold on
            traces = StimTraces{aa};
            meantrace = mean(traces,1);
            plot(traces','Linewidth',0.5)
            ylabel(['Stim index ' num2str(aa)])
            if vv == 1; title(['Neighbor ' num2str(a)]); end
            subtitle([num2str(round(NeighborData.V1(aa),1)), 'kHz ', num2str(round(NeighborData.V2(aa),1)), 'dB'])
            if CellDataByStim.RF(aa) == 1
                plot(meantrace,'m','Linewidth',2)
            elseif NeighborData.Ztest(aa) == 1 && NeighborData.BestNeighbor(aa) == v
                plot(meantrace,'g','Linewidth',2)
            else
                plot(meantrace,'k','Linewidth',2)
            end
        end
    end            
end
