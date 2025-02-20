function [ReliabilityData] = simple_reliability(trials, blanks, StimInfo, Ops)
% This function outputs a response reliability measure (how reliably does the
% response encode the stimulus?)

% Reliability measure is based on Sadeh and Cleopath (BioRxiv, 2021)
% It evaluates the Pearson's correlation for each neuron for each distinct
% stimulus across all pairs of repetitions

% We additionally included an analysis of shuffled data (number of shuffles determined in
% user data)using either 10 random blank trials or using trials that are circularly shifted 
% using the matlab function (circshift) starting at random points. We can determine the
% 'significance' of our reliability metric as the z_score and p of our reliability metric from 
% the shuffled distributions. 

% Data is kept in 'single' format wherever possible to reduce computation load

% Argument(s): 
% trials - individual responses to one unique stimulus
% blanks - blanks from cell being analyzed
% Ops.use_reliability_shift - 'Blanks' (unless no blanks), or 'Shifted' to use circle-shifted traces
% Ops.reliability_p_value - alpha level
% Ops.reliability_shuffles - # of shuffles, suggested 1000 
% Ops.smTime - moving window (s) for smoothing raw traces. Set to [] for no smoothing

% Returns ReliabilityData with:

% RELIABILITY METRIC
% Mode - blank or shifted
% R - reliability of responses measured by Pearson's correlation

% STATISTICS FOR COMPARISON AGAINST BLANK OR SHIFTED TRIALS
% IsReliable = 1 means R is from different distribution (rejects null)
% z_P = p-value of the test, R against control distribution
% z_CI = confidence interval
% z_R = z-score from control distriubtion
% R_dist = control distribution of reliability measures 

% TO PLOT: Use plot_simple_reliability script

%% Default options

if nargin < 4
    Ops = struct;
    Ops.use_reliability_shift = 'Blanks';
    Ops.reliability_p_value = 0.025;
    Ops.reliability_shuffles = 1000;
    Ops.smTime = 0;
    Ops.ComputeReliabilityControls = 1;
    Ops.PlotPearsonFigure = 0;
end

if ~isfield(Ops, 'PlotPearsonFigure')
    Ops.PlotPearsonFigure = 0;
end

if ~isfield(Ops, 'ComputeReliabilityControls')
    Ops.ComputeReliabilityControls = 1;
end

r_shuffles = Ops.reliability_shuffles;

%% Setup
N1 = size(trials,1); % number of trials
N2 = size(trials,2); % number of frames/trial
N4 = size(blanks,1); % number of blank trials

%% Reliability mode
% 0 = use blanks, 1 = use shifted data

%If no blank trials or less blank trials than real trials, default to shifted data instead
if strcmp(Ops.use_reliability_shift, 'Blanks')   
    if isempty(blanks) || N4 <= N1
        Ops.use_reliability_shift = 'Shifted';
    end
end

mode = Ops.use_reliability_shift;

%% Initialize reliability table

ReliabilityData = struct;
ReliabilityData.Mode = mode;
ReliabilityData.R = single(nan); 
ReliabilityData.IsReliable = single(nan); 
ReliabilityData.z_P = single(nan); 
ReliabilityData.z_CI = single(nan(1,2));
ReliabilityData.z_R = single(nan); 
ReliabilityData.R_dist = single(nan); 
ReliabilityData.Pearsons = [];
    
%If trials are all NaNs, all zeroes, empty, or there is only 1 trial -- this will skip function
if ~any(trials)
    return
end

if size(trials,1) <= 1
    return
end

%% Get data from StimInfo

fs = StimInfo.fs;
soundtime = StimInfo.baseline*fs;

%Optional smooth traces with a moving window
if Ops.smTime > 0
    smoothval = ceil(Ops.smTime/(1/fs)); %window size
    trials = single(smoothdata(trials, 2,'movmean',smoothval));
    blanks = single(smoothdata(blanks, 2,'movmean',smoothval));
end

%Remove baseline (this improves correlation values)
trials = trials(:,soundtime+1:end); %Trials and blanks input from simple_extract_data are already single
blanks = blanks(:,soundtime+1:end);

%% Measure reliability of stimulus responses

%Compute Pearson's correlation
Pearsons = corrcoef(trials','Rows','pairwise'); %Pairwise flag needed to omit nans
Pearson_AB = Pearsons(find(tril(Pearsons,-1))); %Reshape and keep values below the diagonal
R = mean(Pearson_AB);
ReliabilityData.R = R;
ReliabilityData.Pearsons = Pearsons;
% note: the normalization for reliability is different than Sadeh and Cleopath - 
% think it was a typo in the formula in their paper

if R > 0.3
    disp('')
end

if Ops.PlotPearsonFigure
    figure;
    index = 1;
    N3 = factorial(N1)/(2*factorial(N1-2)); % number of pairs (N!/(2*(N-2))!)

    for i = 1:N1-1
        A = trials(i,:);
        for j = 1:N1-i 
            B = trials(i+j,:); % pairwise comparison with following traces
            Pearsons = corrcoef(A,B);
            Pearson_AB(index) = Pearsons(2,1);
            subplot(N3,1,index); plot(A, 'r'); hold on; plot(B, 'b');
            title(strcat({'Pearsons Trial '}, num2str(i), {' vs. '}, num2str(i+j), {' = '}, num2str(Pearson_AB(index))));
            % To check, this Pearson's formula gives the same numbers
        %         A_norm = (A-mean(A))./std(A);
        %         B_norm = (B-mean(B))./std(B);
        %         C = A_norm.*B_norm; % sum of multiplied frames
        %         P_AB(index) = (1/(N2-1)).*sum(C); % pairwise pearson's
            index = index+1;
        end
    end
    
    R = mean(Pearson_AB);
end

%% Return function now if we don't want to compute controls
% We will have the R value and Pearsons values stored

if ~Ops.ComputeReliabilityControls
    return
end

%% Measure reliability of blanks as control and determine Z-score

if strcmp(mode, 'Blanks')
    
    %Compute Pearson's correlation for all possible combinations of blanks
    Pearsons = corrcoef(blanks','Rows','pairwise'); %Pairwise flag needed to omit nans

    %Generate list of random traces
    rand_blank_trial_list = randsample(1:N4,N1*r_shuffles,'true'); %have to use replacement to do this outside of loop
    shuffle_inds = sort(repmat((1:r_shuffles)',[N1,1])); %Indices corresponding to shuffle iteration

    %Pick N1 random traces to compute the mean Pearsons from & do this r_shuffles number of times
    reliability_blanks = nan(r_shuffles,1);
    for u = 1:r_shuffles
        rand_blank_trials = rand_blank_trial_list(shuffle_inds == u);
        %Make sure that there are no duplicate trials
        if length(unique(rand_blank_trials)) < N1
            rand_blank_trials = randsample(1:N4,N1,'false');
        end

        %Random index of Pearson's matrix
        pairs = nchoosek(rand_blank_trials,2); % all combinations
        ind = sub2ind(size(Pearsons),pairs(:,1),pairs(:,2)); %

        Pearson_AB_blanks = Pearsons(ind);
        reliability_blanks(u) = mean(Pearson_AB_blanks);
    end
    
    % perform one-tailed z test to determine if reliability metric R is from same distribution as reliability measures of blank trials
    % tail is 'right' because we have hypothesis R will be greater than blank distribution
    [z_test_blank, z_P_blank, z_CI_blank, z_R_blank] = ztest(R, mean(reliability_blanks), std(reliability_blanks), 'Alpha', Ops.reliability_p_value, 'tail', 'right');
    
    % eliminate those with R<0 from passing test
    if R<0
        z_test_blank = 0;
    end 
    % the above z-score is the same as how we calculate z scores (below)
    % z_score = (R-mean(reliability_blanks))./std(reliability_blanks)
    
    % Store results
    ReliabilityData.IsReliable = single(z_test_blank);
    ReliabilityData.z_P = single(z_P_blank);
    ReliabilityData.z_CI = single(z_CI_blank);
    ReliabilityData.z_R = single(z_R_blank);
    ReliabilityData.R_dist = single(reliability_blanks);
end

%% Measure reliability of shifted traces as control and determine Z-score

if strcmp(mode, 'Shifted')

    %FYI: In an earlier version of the code, I tried replicating the trials r_shuffles times and performing all of the shifting and corrcoef
    %computation ahead of this loop, however this was actually slower than doing it in a loop because the corrcoef function got very slow for that many trials
    
    reliability_shifted = zeros(r_shuffles,1,'single');
    for u = 1:r_shuffles
        %Shift the trials by random frame 
        rand_shifting = randsample(1:N2,N1,'true')';  % random spots within traces to start shift (have to use replacement to do this outside of loop)
        shifted_trials = cell2mat(arrayfun(@(x) circshift(trials(x,:),[1 rand_shifting(x)]),(1:numel(rand_shifting))','un',0));

        %Compute Pearson's correlation
        Pearson_shuffle = corrcoef(shifted_trials','Rows','pairwise'); %Pairwise flag needed to omit nans
        Pearson_AB_shifted = Pearson_shuffle(find(tril(Pearson_shuffle,-1))); 
        reliability_shifted(u) = mean(Pearson_AB_shifted);
    end

    % tail is 'right' because we have hypothesis R will be greater than blank distribution
    [z_test_shift, z_P_shift, z_CI_shift, z_R_shift] = ztest(R, mean(reliability_shifted), std(reliability_shifted), 'Alpha', Ops.reliability_p_value, 'tail', 'right');
    
    % eliminate those with R<0 from passing test
    if R<0
        z_test_shift = 0;
    end 

    % Store results
    ReliabilityData.IsReliable = single(z_test_shift);
    ReliabilityData.z_P = single(z_P_shift);
    ReliabilityData.z_CI = single(z_CI_shift);
    ReliabilityData.z_R = single(z_R_shift);
    ReliabilityData.R_dist = single(reliability_shifted);
end

%FYI: To plot result use plot_simple_reliability script
                    
end

