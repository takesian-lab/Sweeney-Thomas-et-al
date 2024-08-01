function [XCorrTraces, cc_mat, lag_in_s] = simple_xcorr_network(traces,XCorrOps)

% This function outputs cross correlations of full traces of all neurons in a network 
% Based on trace_correlation originally written by Anne & simple_xcorr
% written by Maryse Thomas

% Argument(s): 
% traces - full raw traces from each cell in block
% maxlag_in_s - max lag (in seconds) for cross correlation
% fs - recording framerate
% shuffles - number of shuffles for control distribution for stats to
% determine significance of correlation

% Returns XcorrTraces with:
%   cc_max = average max of the xcorr functions for all unique trial combinations
%   cc_min = average min of the xcorr functions...
%   lag_max = average lag time of max (s)
%   lag_min = average lag time of min (s)
%   cc_zero = average value at zero time lag of the xcorr funtions
%   " " _z = Z-test result (0 or 1) for comparing each of the above values against control distribution
%   " " _p = p-value of the test
%   " " _r = z-score from the test


%% If trials are all NaNs, all zeroes, or empty -- this will skip function
if ~any(traces)
    return
end

%% Get data from StimInfo
maxlag = floor(XCorrOps.fs*XCorrOps.maxlag_in_s);

N1 = size(traces,1); % number of traces
N2 = size(traces,2); % number of frames/trace

% Optional smooth traces with a moving window
if XCorrOps.smTime > 0
    smoothval = ceil(XCorrOps.smTime/(1/XCorrOps.fs)); %window size
    traces = single(smoothdata(traces, 2,'movmean',smoothval));
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
template.cc_max_z = nan(nPairs,1);    %xcorr z test compared to cntrl
template.cc_max_p = nan(nPairs,1);    %xcorr p value compared to cntrl
template.cc_max_r = nan(nPairs,1);    %xcorr r value compared to cntrl
template.cc_zero_z = nan(nPairs,1);    %xcorr z test @ 0 lag compared to cntrl
template.cc_zero_p = nan(nPairs,1);    %xcorr p value @ 0 lag compared to cntrl
template.cc_zero_r = nan(nPairs,1);    %xcorr r value @ 0 lagcompared to cntrl

%% Convert raw traces to zscore or df/f
if XCorrOps.zscore
    traces = zscore(traces,'',2);
else
    mean_F = mean(traces,2); % Convert raw  traces to df/f
    traces = (traces - mean_F)./mean_F; %(total-mean)/mean
end

%% Loop through each cell pair and save cross-correlations

%Compute cross-correlogram for all cell pairs at once
[cc_mat, lag] = xcorr(traces', maxlag, 'coeff');

%Transpose cc_mats and only keep unique pairs
cc_mat = cc_mat(:,unique_rows)';

%From this point onwards, we don't need to compute lag anymore because it will always be the same
%lag_full_in_s = lag_full*(1/fs); %Convert to seconds for plotting
lag_in_s = lag*(1/XCorrOps.fs); %Convert to seconds for plotting

%Compute xcorr measures for each trial
XCorrTraces = template; %T stands for trials
[XCorrTraces.cc_max, max_ind] = max(cc_mat,[],2,'omitnan');
[XCorrTraces.cc_min, min_ind] = min(cc_mat,[],2,'omitnan');
XCorrTraces.lag_max = lag_in_s(max_ind)';
XCorrTraces.lag_min = lag_in_s(min_ind)';
XCorrTraces.cc_zero = cc_mat(:,lag == 0);

%Compute average final trace for plotting
cc = mean(cc_mat,1,'omitnan');

%Store mean xcorr measures 
XcorrMeanData.cc_max = mean(XCorrTraces.cc_max, 'omitnan');
XcorrMeanData.cc_min = mean(XCorrTraces.cc_min, 'omitnan');
XcorrMeanData.lag_max = mean(abs(XCorrTraces.lag_max), 'omitnan'); %Take mean abs value here
XcorrMeanData.lag_min = mean(abs(XCorrTraces.lag_min), 'omitnan'); %Take mean abs value here
XcorrMeanData.cc_zero = mean(XCorrTraces.cc_zero, 'omitnan');

%% Compare vs. circle shifted traces
    
%Replicate trials number of times needed for nShuffles
repeated_traces = repmat(traces,[XCorrOps.shuffles,1]);

%Shift the replicated traces by random frame
rand_shifting = randsample(1:N2,N1*XCorrOps.shuffles,'true')';  % random spots within traces to start shift (have to use replacement to do this outside of loop)
shifted_traces = cell2mat(arrayfun(@(x) circshift(repeated_traces(x,:),[1 rand_shifting(x)]),(1:numel(rand_shifting))','un',0));

%Compute nShuffle x nPair xcorr functions and make distribution of averages
C = struct;
shuffle_ind = 1;
for u = 1:XCorrOps.shuffles
    %shuffle_index = shuffle_ind:shuffle_ind+N1-1;
    shuffle_index = shuffle_ind:shuffle_ind-1+N1;
    shuffled_traces = shifted_traces(shuffle_index,:);

    %Compute cross-correlogram for all traces at once
    ctrl_cc_mat = xcorr(shuffled_traces', maxlag, 'coeff');

    %Transpose cc_mats and only keep unique pairs
    ctrl_cc_mat = ctrl_cc_mat(:,unique_rows)';

    %Compute xcorr measures for each trial
    B = template;
    [B.cc_max, max_ind] = max(ctrl_cc_mat,[],2,'omitnan');
    [B.cc_min, min_ind] = min(ctrl_cc_mat,[],2,'omitnan');
    B.lag_max = lag_in_s(max_ind)';
    B.lag_min = lag_in_s(min_ind)';
    B.cc_zero = ctrl_cc_mat(:,lag == 0);

    C.cc_max(u,:) = B.cc_max';
    C.cc_zero(u,:) = B.cc_zero';
    shuffle_ind = shuffle_ind + N1;
end

%% Perform z test on all measures to determine if trace correlations is from same distribution as shuffled trials

%Compute average final trace for plotting
variables = {'cc_max', 'cc_zero'};
% tail = {'both', 'both', 'both', 'both', 'both'}; %tails for distribution
    
for v = 1:length(variables)
    for pairs = 1:nPairs

    val = XCorrTraces.(variables{v})(pairs);
    ctrl_dist = C.(variables{v})(:,pairs);

    [z,p,CI,R] = ztest(val, mean(ctrl_dist), std(ctrl_dist), 'Alpha', XCorrOps.p_value); %, 'tail', tail{v});

%Store results
    XCorrTraces.(strcat(variables{v}, '_z'))(pairs,:) = z;
    XCorrTraces.(strcat(variables{v}, '_p'))(pairs,:) = p;
    XCorrTraces.(strcat(variables{v}, '_r')) (pairs,:)= R;
    end
end
 

%% Optional plot to verify correlations

if XCorrOps.plot_correlations
    
    plot_only_sig = 1;
    
    %first two plots will be the average traces and average correlations
    figure; tiledlayout(1,3) 
    nexttile([1 2]); hold on
    plot(mean(traces,1), 'k', 'Linewidth', 2)
    shadedErrorBar(1:size(traces,2), mean(traces,1), std(traces,[],1))
    xlim([1 size(traces,2)])
    title('Average trace')

    nexttile; hold on
    plot(lag_in_s(1,:), mean(cc_mat,1), 'k', 'Linewidth', 3)
    xlim([lag_in_s(1) lag_in_s(end)])
    %ylim([-0.5 0.5])
    title(['Average max = ' num2str(XcorrMeanData.cc_max) ', lag ' num2str(XcorrMeanData.lag_max) ', cc_0 = ' num2str(XcorrMeanData.cc_zero)]);

    %all the other plots will be each pair
    nPlots = 10;
    nBins = ceil(nPairs/nPlots);

    for b = 1:nBins
        currentBin = (b-1)*(nPlots) + (1:nPlots);
        currentBin(currentBin > nPairs) = [];
        figure; tiledlayout(nPlots,3)
        
        for i = 1:length(currentBin)
            ii = currentBin(i);

            if XCorrTraces.cc_zero_z(ii) == 0 && plot_only_sig == 1
                print 'Skip plotting non-significant pairs'
            else      
                n1 = unique_pairs(ii,1);
                n2 = unique_pairs(ii,2);

                nexttile; hold on % first tile is z-scored traces
                plot(traces(n1,:), 'm');
                plot(traces(n2,:), 'b');
                xlim([1 size(traces,2)])
                title(['Traces ' num2str(n1) ' & ' num2str(n2)])
                if i == nPlots; xlabel('Time (s)'); end

                nexttile; hold on % second tile is correlations across lags
                plot(lag_in_s, cc_mat(ii,:), 'k', 'Linewidth', 3);
                scatter(XCorrTraces.lag_max(ii), XCorrTraces.cc_max(ii),100,'r')
                ylim([-0.5 0.5])
                title(['max = ' num2str(XCorrTraces.cc_max(ii)) ', lag ' num2str(XCorrTraces.lag_max(ii)) ', cc_0 = ' num2str(XCorrTraces.cc_zero(ii))]);
                if i == nPlots; xlabel('Lag (s)'); end

                nexttile; hold on % third tile is distribution of control shuffles compared to zero-lag correlation of the pair
                min_g = min(min(C.cc_zero(:,ii)), XCorrTraces.cc_zero(ii))-0.1;
                max_g = max(max(C.cc_zero(:,ii)), XCorrTraces.cc_zero(ii))+0.1;
                histFig = histogram(C.cc_zero(:,ii),min_g:.01:max_g,'facecolor','b','facealpha',0.2,'edgecolor','none');
                [maxcount, whichbin] = max(histFig.Values); hold on;
                plot([XCorrTraces.cc_zero(ii) XCorrTraces.cc_zero(ii)],[0 maxcount], '--r');
                title(['Correlations at Zero Lag, r = ' num2str(XCorrTraces.cc_zero_r(ii)) ', p = ' num2str(XCorrTraces.cc_zero_p(ii))]);
            end
        end
    end
end

end

