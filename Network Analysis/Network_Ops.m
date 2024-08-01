%% Anne's Ops for Network Analysis

%% Compile Ops
CompileOps = struct; 
CompileOps.stimTypes = {'FM','NoiseITI','RF', 'SAMNoise','SAMfreq','Spontaneous'};
CompileOps.data_path = 'Z:\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\Network Analysis_NoiseCorr';
CompileOps.avg_across_stim = 1; % run script to average across all stim (1) or skip (0), runtime is 1-2hrs

%% Network Ops
NetworkOps = struct;

% Setup user-specific file paths
NetworkOps.data_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice';
NetworkOps.save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\Network Analysis_NoiseCorr';
NetworkOps.matchFile_path = ''; %Path and filename ending in .xlsx

% General setup
NetworkOps.use_responsive_cells = 1; % 0 = use all cell, 1 = use only responsive cells to given stimulus (n/a for spont analysis)
NetworkOps.stimTypes = {'Spontaneous'}; %, 'NoiseITI', 'FM', 'SAMfreq','SAMnoise'}; %, 'Spontaneous', 'NoiseITI', 'RF', 'FM', 'SAM', 'SAMfreq'};  % input all stim types to analyze together {'Spontaneous','NoiseITI', 'RF','FM','SAM','SAMfreq'}
NetworkOps.bins = [-1:0.025:1]; % bins for histograms
NetworkOps.histogram_type = 'probability'; % histogram type, either 'probability' or 'count'
NetworkOps.histogram_line = 0; % 0 is standard histogram, 1 is with lines only
NetworkOps.save_figures = 1; % 0 = don't save figures, 1 = save figures
NetworkOps.locoThresh = 0.7057; % threshold of wheel speed to identify running bouts versus stationary periods
NetworkOps.shuffles = 1000; % number of shuffles to determine significance for correlation analysis

% Setup for XCorr (full trace) correlations
NetworkOps.Xcorr.ComputeXcorr = 1; %0 or 1 to compute trace correlations
NetworkOps.XCorr.smTime = 0; % Optional smooth traces with a moving window
NetworkOps.XCorr.plot_correlations = 0; % Figure with all pairwise correlations (for troubleshooting)
NetworkOps.XCorr.plot_only_sig = 0; % Plot only significant correlations (zero-lag)
NetworkOps.XCorr.p_value = 0.05; % p value to determine if correlations are significant
NetworkOps.XCorr.maxlag_in_s = 20; % maximum lag to evaluate cross correlations on traces (in seconds)
NetworkOps.XCorr.shuffles = 1000; % number of shuffles to determine significance for trace correlation analysis
NetworkOps.XCorr.zscore = 1; % 0 = df/f, 1 = zscore

% Setup for MotorCorr (correlations with locomotion)
NetworkOps.MotorCorr.ComputeMotorCorr = 1; % 0 or 1 to compute loco correlations
NetworkOps.MotorCorr.smTime = 0; % Optional smooth traces with a moving window
NetworkOps.MotorCorr.plot_correlations = 0;
NetworkOps.MotorCorr.p_value = 0.05;
NetworkOps.MotorCorr.maxlag_in_s = 20; 
NetworkOps.MotorCorr.shuffles = 1000;
NetworkOps.MotorCorr.zscore = 1; % 0 = df/f, 1 = zscore

% Setup for pairwise cell signal & noise correlations
NetworkOps.NoiseCorr.ComputeNoiseCorr = 1; % 0 or 1 to compute noise and signal correlations using pairwise Pearson correlation
NetworkOps.NoiseCorr.ComputeNoiseXCorr = 1; % 0 or 1 to compute noise and signal correlations using pairwise XCorr (max across lags)
NetworkOps.NoiseCorr.trial_lim = 3; % Minimum number of trials that have to exist to find correlations
NetworkOps.NoiseCorr.shuffles = 1000;
NetworkOps.NoiseCorr.p_value = 0.05;
NetworkOps.NoiseCorr.ComputeControls = 1; %1 to compute controls (with use_correlation_shift as mode), 0 to skip (faster)
NetworkOps.NoiseCorr.use_correlation_shift = 0;  % 0 = blank shuffled (unless no blanks), 1 = shifted shuffles
NetworkOps.NoiseCorr.plot_traces = 0; % don't plot traces
NetworkOps.NoiseCorr.plot_histograms = 0; % don't plot histograms
NetworkOps.NoiseCorr.only_red_cells = 0; % 0 = all cells, 1 = find correlations only between red cells
NetworkOps.NoiseCorr.match_RFs = 0; % for noise/signal correlations, only compare stimuli within RF for which both cells are responsive (1 = match RFs) otherwise, 0 = don't match RFs, 
NetworkOps.NoiseCorr.maxlag = 1; % maxlag in seconds for trial noisecorr using xcorr (amount of lag allowed)

% Setup for visualizing networks
NetworkOps.Visualize.selected_cell = [1]; % cell which you'd like to visualize correlations for saved graph
NetworkOps.Visualize.plot_graph = 1;
NetworkOps.Visualize.sig_only = 1; % only plot significant correlations with selected cell
NetworkOps.Visualize.use_matched_cell = 0; % block cells = 0, matched cells =1

% Setup for Summary Plots
NetworkOps.Summary.bins = [-1:0.025:1]; % bins for histograms
NetworkOps.Summary.histogram_type = 'count'; %'count';  % histogram type
NetworkOps.Summary.histogram_line = 0; %0 is standard histogram, 1 is with lines
NetworkOps.Summary.plot_stability = 0;
NetworkOps.Summary.save_figures = 1;
NetworkOps.Summary.nshuffles = 1000;
NetworkOps.Summary.timepoints = 3; % # of timepoints to compare across

% sort by cell type, stim type
NetworkOps.Summary.all_celltypes = {'PYR','NDNF','VIP'}; % NDNF, VIP, PYR or '' if no sort
NetworkOps.Summary.stim = ''; %Spontaneous'; % RF, FM, NoiseITI or '' if no sort



