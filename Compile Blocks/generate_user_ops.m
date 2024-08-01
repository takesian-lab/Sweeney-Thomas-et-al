% List of Suite2p ops to check

% This script saves a user_ops.mat file with a list of Suite2p ops that the 
% user would like to check against their compiled blocks. This check will
% happen during define_suite2p_singleblock.

% Each user should set their own ops file

% Steps:
% 1. Save script as ops_user
% 2. Modify filepath and filename
% 3. Modify list of ops [comment those you don't want to check]
% 4. Run script

% The list of ops is based on the Suite2p GitHub version from May 2020

%% Save filepath

filepath = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
filename = 'Maryse_ops_thy1.mat';

user_ops = struct;

%% File paths
%ops.look_one_level_down
%ops.save_path0
%ops.fast_disk
user_ops.bruker = 0;

%% Main settings
user_ops.nplanes = 1;
user_ops.nchannels = 1;
user_ops.functional_chan = 1;
user_ops.tau = 1.5;
user_ops.fs = 30;
user_ops.do_bidiphase = 1;
user_ops.bidiphase = 0;

%% Output settings
user_ops.preclassify = 0;
%ops.save_mat = 1;
user_ops.combined = 1;
%ops.reg_tif
%ops.reg_tif_chan2
user_ops.aspect = 1;
%ops.delete_bin
%ops.move_bin

%% Registration
user_ops.do_registration = 1;
user_ops.align_by_chan = 1;
user_ops.nimg_init = 300;
user_ops.batch_size = 500;
user_ops.smooth_sigma = 1.15;
user_ops.smooth_sigma_time = 0;
user_ops.maxregshift = 0.1;
user_ops.th_badframes = 1;
%ops.keep_movie_raw
%ops.two_step_registration

%% Nonrigid
user_ops.nonrigid = 1;
user_ops.block_size = [128, 128];
user_ops.snr_thresh = 1.2;
user_ops.maxregshiftNR = 5.0;

%% 1P
%There is no ops record of 1Preg yet (or at least not if it is set to 0)
%Keep an eye on this and make it's not 1 when setting ops.
%ops.spatial_hp = 50;
%ops.pre_smooth = 2;
%ops.spatial_taper = 50;

%% ROI detection
user_ops.roidetect = 1;
%ops.sparse_mode = 0;
%ops.diameter = 10;
user_ops.spatial_scale = 0;
user_ops.connected = 1;
user_ops.threshold_scaling = 1;
user_ops.max_overlap = 0.75;
user_ops.max_iterations = 20;
user_ops.high_pass = 100;

%% Extraction\Neuropil
user_ops.allow_overlap = 0;
user_ops.inner_neuropil_radius = 3;
user_ops.min_neuropil_pixels = 350;

%% Deconvolution
user_ops.spikedetect = 1;
user_ops.win_baseline = 60;
user_ops.sig_baseline = 10;
user_ops.neucoeff = 0.7;

%% Save ops
cd(filepath)
save(filename, 'user_ops', '-mat');

