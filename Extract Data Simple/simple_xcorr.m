function [XcorrData] = simple_xcorr(traces, blanks, StimInfo, Ind, Ops)
% This function outputs a response reliability measure based on xcorr
% The goal is that xcorr can quantify reliability even given stim responses that are not well time-locked
% Based on trace_correlation originally written by Anne

% Argument(s): 
% traces - individual trial responses to one unique stimulus
% blanks - blank trials from cell being analyzed
% StimInfo - from simple_extract_data
% Ind - stimulus index corresponding to StimInfo.Combined
% Ops.use_reliability_shift - 'Blanks' (unless no blanks), or 'Shifted' to use circle-shifted traces
% Ops.reliability_p_value - alpha level
% Ops.reliability_shuffles - # of shuffles, suggested 1000 
% Ops.smTime - moving window (s) for smoothing raw traces. Set to [] for no smoothing

% Returns XcorrData with:
%   Mode - Blanks or Shifted
%   cc_max = average max of the xcorr functions for all unique trial combinations
%   cc_min = average min of the xcorr functions...
%   lag_max = average lag time of max (s)
%   lag_min = average lag time of min (s)
%   cc_zero = average value at zero time lag of the xcorr funtions
%   " " _z = Z-test result (0 or 1) for comparing each of the above values against control distribution
%   " " _p = p-value of the test
%   " " _r = z-score from the test

%% Setup

%If Ops not specified, set default Ops
if nargin < 5
    Ops = struct;
    Ops.use_reliability_shift = 'Blanks';
    Ops.reliability_p_value = 0.05;
    Ops.reliability_shuffles = 1000;
    Ops.smTime = 0;
    Ops.ComputeReliabilityControls
end

maxlag_in_s = 0.1; %Time window for computing correlations (in seconds)
plot_correlations = 0; %Figure with all pairwise correlations (for troubleshooting)

%% Initialize xcorr table

XcorrData = table;
XcorrData.Mode = "";
XcorrData.cc_max = single(nan);
XcorrData.cc_min = single(nan);
XcorrData.lag_max = single(nan);
XcorrData.lag_min = single(nan);
XcorrData.cc_zero = single(nan);
XcorrData.cc_max_z = single(nan);
XcorrData.cc_max_p = single(nan);
XcorrData.cc_max_r = single(nan);
XcorrData.cc_min_z = single(nan);
XcorrData.cc_min_p = single(nan);
XcorrData.cc_min_r = single(nan);
XcorrData.lag_max_z = single(nan);
XcorrData.lag_max_p = single(nan);
XcorrData.lag_max_r = single(nan);
XcorrData.lag_min_z = single(nan);
XcorrData.lag_min_p = single(nan);
XcorrData.lag_min_r = single(nan);
XcorrData.cc_zero_z = single(nan);
XcorrData.cc_zero_p = single(nan);
XcorrData.cc_zero_r = single(nan);
XcorrData.cc_max_dist = cell(1,1);
XcorrData.cc_zero_dist = cell(1,1);

%If trials are all NaNs, all zeroes, empty, or there is only one trial -- this will skip function
if ~any(traces)
    return
end

if size(traces,1) <= 1
    return
end

%% Get data from StimInfo

stimlabel = StimInfo.Combined_Label{Ind};
xlabel_increment = 0.5; %seconds
fs = StimInfo.fs;
nFrames = StimInfo.nFrames;
nBaselineFrames = StimInfo.baseline*fs;

%Remove baseline and determine max lag in frames
nFrames = nFrames - nBaselineFrames;
traces = traces(:,nBaselineFrames+1:end);
blanks = blanks(:,nBaselineFrames+1:end);
maxlag = floor(fs*maxlag_in_s);

%If maxlag is longer than trial duration, set equal to trial duration
if maxlag > nFrames 
    maxlag = nFrames;
end

%X labels for figures
increment_in_frames = xlabel_increment*fs;
xtick = 0:increment_in_frames:nFrames;
xticklabel = 0:xlabel_increment:((1/fs)*nFrames);

%% Reliability mode
% 0 = use blanks, 1 = use shifted data

N1 = size(traces,1); % number of trials
N2 = size(traces,2); % number of frames/trial
N4 = size(blanks,1); % number of blank trials

%If no blank trials or less blank trials than real trials, default to shifted data instead
if strcmp(Ops.use_reliability_shift, 'Blanks')   
    if isempty(blanks) || N4 <= N1
        XcorrData.Mode = "Shifted";
    else
        XcorrData.Mode = "Blanks";
    end
else
    XcorrData.Mode = "Shifted";
end

%Optional smooth traces with a moving window
if Ops.smTime > 0
    smoothval = ceil(Ops.smTime/(1/StimInfo.fs)); %window size
    traces = single(smoothdata(traces, 2,'movmean',smoothval));
    blanks = single(smoothdata(blanks, 2,'movmean',smoothval));
end

%% Figure out how many unique pairs there are

pairs = [];
count = 1;
for n1 = 1:size(traces,1)
    for n2 = 1:size(traces,1)
        pairs(count,1:2) = [n1, n2];
        count = count + 1;
    end
end
unique_pairs = unique(sort(pairs')', 'rows');
unique_pairs(unique_pairs(:,1) == unique_pairs(:,2),:) = [];
nPairs = size(unique_pairs,1); %equivalent to nchoosek(size(traces,1),2)

unique_rows = nan(nPairs,1);
for i = 1:nPairs
    unique_rows(i) = find(sum((pairs - unique_pairs(i,:)) == 0, 2) == 2);
end

%% Make template tables that will be used to store intermediate xcorr results and control distributions

template = table;
template.cc_max = nan(nPairs,1);     %xcorr max values
template.cc_min = nan(nPairs,1);     %xcorr min values
template.lag_max = nan(nPairs,1);    %latency corresponding to xcorr max
template.lag_min = nan(nPairs,1);    %latency corresponding to xcorr min
template.cc_zero = nan(nPairs,1);    %xcorr value at 0 time lag

shuffle_template = table;
shuffle_template.cc_full = cell(Ops.reliability_shuffles,1);     %full xcorr
shuffle_template.cc = cell(Ops.reliability_shuffles,1);          %xcorr
shuffle_template.cc_max = nan(Ops.reliability_shuffles,1);       %xcorr max values
shuffle_template.cc_min = nan(Ops.reliability_shuffles,1);       %xcorr min values
shuffle_template.lag_max = nan(Ops.reliability_shuffles,1);      %latency corresponding to xcorr max
shuffle_template.lag_min = nan(Ops.reliability_shuffles,1);      %latency corresponding to xcorr min
shuffle_template.cc_zero = nan(Ops.reliability_shuffles,1);      %xcorr value at 0 time lag

%% Loop through each trial pair and save cross-correlations

%Compute cross-correlogram for all trials at once
[cc_mat_full, lag_full] = xcorr(traces', size(traces,2), 'coeff');
[cc_mat, lag] = xcorr(traces', maxlag, 'coeff');

%Transpose cc_mats and only keep unique pairs
cc_mat_full = cc_mat_full(:,unique_rows)';
cc_mat = cc_mat(:,unique_rows)';

%From this point onwards, we don't need to compute lag anymore because it will always be the same
lag_full_in_s = lag_full*(1/fs); %Convert to seconds for plotting
lag_in_s = lag*(1/fs); %Convert to seconds for plotting

%Compute xcorr measures for each trial
T = template; %T stands for trials
[T.cc_max, max_ind] = max(cc_mat,[],2,'omitnan');
[T.cc_min, min_ind] = min(cc_mat,[],2,'omitnan');
T.lag_max = lag_in_s(max_ind)';
T.lag_min = lag_in_s(min_ind)';
T.cc_zero = cc_mat(:,lag == 0);

%Compute average final trace for plotting
cc_full = mean(cc_mat_full,1,'omitnan');
cc_full_std = std(cc_mat_full,[],1,'omitnan');
cc = mean(cc_mat,1,'omitnan');

%Store mean xcorr measures 
XcorrData.cc_max = mean(T.cc_max, 'omitnan');
XcorrData.cc_min = mean(T.cc_min, 'omitnan');
XcorrData.lag_max = mean(abs(T.lag_max), 'omitnan'); %Take mean abs value here
XcorrData.lag_min = mean(abs(T.lag_min), 'omitnan'); %Take mean abs value here
XcorrData.cc_zero = mean(T.cc_zero, 'omitnan');

%% Optional plot to verify correlations

if plot_correlations
    
    nPlots = 10;
    nBins = ceil(nPairs/nPlots);
    
    for b = 1:nBins
        currentBin = (b-1)*(nPlots) + (1:nPlots);
        currentBin(currentBin > nPairs) = [];
        
        figure; tiledlayout(nPlots+1,2)

        %first two plots will be the average results
        nexttile; hold on
        plot(mean(traces,1), 'k', 'Linewidth', 2)
        shadedErrorBar(1:size(traces,2), mean(traces,1), std(traces,[],1))
        xlim([1 size(traces,2)])
        set(gca, 'Xtick', xtick)
        set(gca, 'Xticklabel', xticklabel)
        title('Average trace')

        nexttile; hold on
        plot(lag_full_in_s(1,:), cc_full, 'k', 'Linewidth', 1)
        plot(lag_in_s(1,:), cc, 'k', 'Linewidth', 3)
        shadedErrorBar(lag_full_in_s, cc_full, cc_full_std)
        xlim([lag_full_in_s(1) lag_full_in_s(end)])
        ylim([-0.5 0.5])
        vline(lag_in_s(1))
        vline(lag_in_s(end))
        title(['Average max = ' num2str(XcorrData.cc_max) ', lag ' num2str(XcorrData.lag_max) ', cc_0 = ' num2str(XcorrData.cc_zero)]);

        %all the other plots will be each pair
        for i = 1:length(currentBin)
            ii = currentBin(i);
            n1 = unique_pairs(ii,1);
            n2 = unique_pairs(ii,2);

            nexttile; hold on
            plot(traces(n1,:), 'm');
            plot(traces(n2,:), 'b');
            xlim([1 size(traces,2)])
            set(gca, 'Xtick', xtick)
            set(gca, 'Xticklabel', xticklabel)
            title(['Traces ' num2str(n1) ' & ' num2str(n2)])
            if i == nPlots; xlabel('Time (s)'); end

            nexttile; hold on
            plot(lag_full_in_s, cc_mat_full(ii,:), 'k', 'Linewidth', 1);
            plot(lag_in_s, cc_mat(ii,:), 'k', 'Linewidth', 3);
            scatter(T.lag_max(ii), T.cc_max(ii),100,'r')
            xlim([lag_full_in_s(1) lag_full_in_s(end)])
            ylim([-0.5 0.5])
            vline(lag_in_s(1))
            vline(lag_in_s(end))
            title(['max = ' num2str(T.cc_max(ii)) ', lag ' num2str(T.lag_max(ii)) ', cc_0 = ' num2str(T.cc_zero(ii))]);
            if i == nPlots; xlabel('Lag (s)'); end
        end

        tak_suptitle(stimlabel)
    end
end
   
%% Return function now if we don't want to compute controls
% We will have the R value and Pearsons values stored

if ~Ops.ComputeReliabilityControls
    return
end

%% Compare correlations to blank trials

if strcmp(XcorrData.Mode, 'Blanks')
    
    %Figure out how many possible blank pairs we have
    blank_pairs = [];
    count = 1;
    for n1 = 1:size(blanks,1)
        for n2 = 1:size(blanks,1)
            blank_pairs(count,1:2) = [n1, n2];
            count = count + 1;
        end
    end
    unique_blank_pairs = unique(sort(blank_pairs')', 'rows');
    unique_blank_pairs(unique_blank_pairs(:,1) == unique_blank_pairs(:,2),:) = [];
    nBlankPairs = size(unique_blank_pairs,1); %equivalent to nchoosek(size(blanks,1),2)

    unique_blank_rows = nan(nBlankPairs,1);
    for i = 1:nBlankPairs
        unique_blank_rows(i) = find(sum((blank_pairs - unique_blank_pairs(i,:)) == 0, 2) == 2);
    end

    %Compute cross-correlogram for all blank trials at once
    blank_cc_mat_full = xcorr(blanks', size(blanks,2), 'coeff');
    blank_cc_mat = xcorr(blanks', maxlag, 'coeff');

    %Transpose cc_mats and only keep unique pairs
    blank_cc_mat_full = blank_cc_mat_full(:,unique_blank_rows)';
    blank_cc_mat = blank_cc_mat(:,unique_blank_rows)';

    %Compute nShuffle x nPair xcorr functions and make distribution of averages
    B = shuffle_template; %B stands for blanks
    
    for u = 1:Ops.reliability_shuffles
        %Choose nPairs randomly among the possible pairs
        %randperm samples without replacement so we cannot accidentally choose the same pair more than once
        randpairs = randperm(nBlankPairs,nPairs);
        ctrl_cc_mat_full = blank_cc_mat_full(randpairs,:);
        ctrl_cc_mat = blank_cc_mat(randpairs,:);
        
        %Compute xcorr measures for each trial
        C = template; %C stands for control
        [C.cc_max, max_ind] = max(ctrl_cc_mat,[],2,'omitnan');
        [C.cc_min, min_ind] = min(ctrl_cc_mat,[],2,'omitnan');
        C.lag_max = lag_in_s(max_ind)';
        C.lag_min = lag_in_s(min_ind)';
        C.cc_zero = ctrl_cc_mat(:,lag == 0);

        %Store full_cc and cc for plotting
        B.cc_full{u} = mean(ctrl_cc_mat_full,1,'omitnan');
        B.cc{u} = mean(ctrl_cc_mat,1,'omitnan');
        
        %Store xcorr measures for average final trace
        B.cc_max(u) = mean(C.cc_max, 'omitnan');
        B.cc_min(u) = mean(C.cc_min, 'omitnan');
        B.lag_max(u) = mean(abs(C.lag_max), 'omitnan');
        B.lag_min(u) = mean(abs(C.lag_min), 'omitnan');
        B.cc_zero(u) = mean(C.cc_zero, 'omitnan');
    end
end

%% Compare vs. circle shifted trials

if strcmp(XcorrData.Mode, 'Shifted')
    
    %Replicate trials number of times needed for nShuffles
    repeated_trials = repmat(traces,[Ops.reliability_shuffles,1]);

    %Shift the replicated trials by random frame 
    rand_shifting = randsample(1:N2,N1*Ops.reliability_shuffles,'true')';  % random spots within traces to start shift (have to use replacement to do this outside of loop)
    shifted_trials = cell2mat(arrayfun(@(x) circshift(repeated_trials(x,:),[1 rand_shifting(x)]),(1:numel(rand_shifting))','un',0));

    %Compute nShuffle x nPair xcorr functions and make distribution of averages
    B = shuffle_template; %B stands for blanks
    
    shuffle_ind = 1;
    for u = 1:Ops.reliability_shuffles
        shuffle_index = shuffle_ind:shuffle_ind+N1-1;
        shuffled_trials = shifted_trials(shuffle_index,:);
        
        %Compute cross-correlogram for all blank trials at once
        ctrl_cc_mat_full = xcorr(shuffled_trials', size(shuffled_trials,2), 'coeff');
        ctrl_cc_mat = xcorr(shuffled_trials', maxlag, 'coeff');

        %Transpose cc_mats and only keep unique pairs
        ctrl_cc_mat_full = ctrl_cc_mat_full(:,unique_rows)';
        ctrl_cc_mat = ctrl_cc_mat(:,unique_rows)';
    
        %Compute xcorr measures for each trial
        C = template; %C stands for control
        [C.cc_max, max_ind] = max(ctrl_cc_mat,[],2,'omitnan');
        [C.cc_min, min_ind] = min(ctrl_cc_mat,[],2,'omitnan');
        C.lag_max = lag_in_s(max_ind)';
        C.lag_min = lag_in_s(min_ind)';
        C.cc_zero = ctrl_cc_mat(:,lag == 0);

        %Store full_cc and cc for plotting
        B.cc_full{u} = mean(ctrl_cc_mat_full,1,'omitnan');
        B.cc{u} = mean(ctrl_cc_mat,1,'omitnan');
        
        %Store xcorr measures for average final trace
        B.cc_max(u) = mean(C.cc_max, 'omitnan');
        B.cc_min(u) = mean(C.cc_min, 'omitnan');
        B.lag_max(u) = mean(abs(C.lag_max), 'omitnan');
        B.lag_min(u) = mean(abs(C.lag_min), 'omitnan');
        B.cc_zero(u) = mean(C.cc_zero, 'omitnan');
        
        shuffle_ind = shuffle_ind + N1;
    end
end

%% Perform z test on all measures to determine if trace correlations is from same distribution as shuffled trials

%Compute average final trace for plotting
ctrl_cc_full = mean(cell2mat(B.cc_full),1,'omitnan');
ctrl_cc_full_std = std(cell2mat(B.cc_full),[],1,'omitnan');
ctrl_cc = mean(cell2mat(B.cc),1,'omitnan');
    
variables = {'cc_max', 'cc_min', 'cc_zero', 'lag_max', 'lag_min'};
tail = {'right', 'both', 'right', 'both', 'both'}; %tail for distribution. We have hypothesis cc_max and cc_zero will be higher, we don't have hypothesis for others
    
for v = 1:length(variables)

    val = XcorrData.(variables{v});
    ctrl_dist = B.(variables{v});
    
    [z,p,CI,R] = ztest(val, mean(ctrl_dist), std(ctrl_dist), 'Alpha', Ops.reliability_p_value, 'tail', tail{v});

    %Store results
    XcorrData.(strcat(variables{v}, '_z')) = single(z);
    XcorrData.(strcat(variables{v}, '_p')) = single(p);
    XcorrData.(strcat(variables{v}, '_r')) = single(R);

    %Store full dist for max and zero
    if strcmp(variables{v}, 'cc_max') || strcmp(variables{v}, 'cc_zero')
        XcorrData.(strcat(variables{v}, '_dist')) = {ctrl_dist};
    end
end
    
%% Summary plot
if plot_correlations
    figure; tiledlayout(3,3)

    %first plot will be the raw traces for the real trials
    nexttile; hold on
    plot(traces')
    plot(mean(traces,1), 'k', 'Linewidth', 2)
    xlim([1 size(traces,2)])
    set(gca, 'Xtick', xtick)
    set(gca, 'Xticklabel', xticklabel)
    title('Raw traces')
    
    %second plot will be the average results for the real traces vs. blanks
    nexttile; hold on
    shadedErrorBar(1:size(blanks,2), mean(blanks,1), std(blanks,1))
    shadedErrorBar(1:size(traces,2), mean(traces,1), std(traces,[],1), 'lineProps','c')
    plot(mean(traces,1), 'c', 'Linewidth', 1)
    xlim([1 size(traces,2)])
    set(gca, 'Xtick', xtick)
    set(gca, 'Xticklabel', xticklabel)
    title('Average trace vs. blanks')

    %third plot will be the average xcorrs
    nexttile; hold on
    shadedErrorBar(lag_full_in_s, cc_full, cc_full_std, 'lineProps','c')
    shadedErrorBar(lag_full_in_s, ctrl_cc_full, ctrl_cc_full_std)
    plot(lag_in_s, ctrl_cc, 'k', 'Linewidth', 3)
    plot(lag_full_in_s, cc_full, 'c', 'Linewidth', 1)
    plot(lag_in_s, cc, 'c', 'linewidth', 3)
    xlim([lag_full_in_s(1) lag_full_in_s(end)])
    ylim([-0.5 0.5])
    vline(lag_in_s(1))
    vline(lag_in_s(end))
    title('Average xcorr vs. blanks')

    %all the other plots will be histograms of the xcorr measures
    for v = 1:length(variables)
        nexttile; hold on
        
        val = XcorrData.(variables{v});
        ctrl_dist = B.(variables{v});
        
        %Set histogram and text colors
        z_test = XcorrData.(strcat(variables{v}, '_z'));
        p = XcorrData.(strcat(variables{v}, '_p'));
        r = XcorrData.(strcat(variables{v}, '_r'));
        
        if strcmp(XcorrData.Mode,'Shifted')
            colour = 'b';
            if z_test == 1
                textcolour = ' {\color{blue}Z=} ';
            else
                textcolour = ' {\color{black}Z=} ';
            end
        elseif strcmp(XcorrData.Mode,'Blanks')
            colour = 'g';
            if z_test == 1
                textcolour = ' {\color{green}Z=} ';
            else
                textcolour = ' {\color{black}Z=} ';
            end
        end
    
        min_g = min(min(ctrl_dist), val)-0.1;
        max_g = max(max(ctrl_dist), val)+0.1;
        %h = histogram(ctrl_dist,min_g:.01:max_g,'facecolor','b','facealpha',0.2,'edgecolor','none');
        h = histogram(ctrl_dist,'facecolor',colour,'facealpha',0.2,'edgecolor','none');
        [maxcount, whichbin] = max(h.Values); hold on;
        plot([val val],[0 maxcount], '--r');
        title(regexprep(variables{v},'_',' ','emptymatch'))
        
        %Test result
        subtitle(strcat(sprintf('%.2f',val), textcolour, sprintf('%.2f', r)), 'FontSize',8);
    end
    tak_suptitle(stimlabel)
end
  
end

