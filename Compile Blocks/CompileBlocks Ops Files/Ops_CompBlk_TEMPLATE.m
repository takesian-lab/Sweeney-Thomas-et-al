% user ops for compile_blocks:

info_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Info sheets v2';
save_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Final info sheets and compiled blocks\Blocks (Final x3)';
info_filename = 'Info_VxDM092022M3';

recompile = 0; %1 to save over previously compiled blocks, 0 to skip
stim_protocol = []; %leave empty [] to ignore, otherwise only compile blocks with this stim protocol
stim_name = '';  %leave empty '' to ignore, otherwise only compile blocks with this stim name
checkOps = 0; %Check Suite2p Fall.ops against user-specified ops
upsampleTo30 = 1; %Upsample 2p data with framerate < 30 to 30Hz (recommended)
WFtrigChange = 0; %Use for widefield imaging with framerate > 20
compute_whiskerpad = 1; %Run define_whiskerpad for experiments collected with video data
EFR = 30; % for behavior experiments where FP or S2p failed/werent used. Choose 20 if matching FP data, 30 for matching S2P data ('estimated frame rate')
look_for_previous_MarkpointsExperiments = 1; % for markpoints experiments, 1 to load prev. defined experiment from file, 0 to redo, 2 to redo except for ensemble name
substitute_ToscaTimes_for_VR = 0; % See salvage_voltagerecording script. To salvage blocks with missing voltage recording data, compile blocks with this option then run salvage script
redo_loco_check = 0; %Set to 0 the first time you are compiling. Afterwards, 1 to redo loco check to ask user how to treat data with no loco activity, 0 to use saved answer
skip_avi = 0; % sometimes you might want to skip avi analysis in Tosca data (errors/troubleshooting/etc.). 1= skip, 0 = default. Whisker will not run if this is set ==1

% setup for Fiber Photometry analysis
Ops.FibPhot.plot_spectrogram = 0; % plot spectrograms? no = 0, yes = 1
Ops.FibPhot.plot_graphs = 0; % plot graphs? no = 0, yes = 1
Ops.FibPhot.params.remove_blips = 0; 
% if 0, will error on data interuptions, 
% if 1, will interpolate data interuptions and give error under block.FibPhot.error
% if 2, will inject NaNs during data interuptions and give error under block.FibPhot.error
Ops.FibPhot.params.end_analysis = 1; % for Fiber Photometry analyis, 0 = continue with error, 1 = abort analysis if spectrogram fails on all channels (eg. no LED on) 
Ops.FibPhot.params.FitControl = 0; % fit blue channel to green to correct for artifacts

Ops.FibPhot.params.freqStep   = 1; % Step size in Hz for spectrogram
Ops.FibPhot.params.freqStepWidth = 7; % Frequency band around peak for spectrogram
Ops.FibPhot.params.inclFreqWin = 4; % Number of frequency bins to average for signal (on either side of peak freq) 
Ops.FibPhot.params.detrendWindowTime = 60; % in seconds
Ops.FibPhot.params.lowPassCorner = 100; % low pass filter - keep at 100 for no filter
Ops.FibPhot.params.notes = []; % establish notes

% final downsampled freqeuncy, if = [], then program will define it
% automatically based on freqeuncies of LEDs 
Ops.FibPhot.params.finalSampleFreq = 30; % []; % this value can't be above 37Hz!
               
%these values are based on acquisition rates and AIN number
Ops.FibPhot.params.numChannels = 10; % number of AIN channels, defined by hardware 
Ops.FibPhot.params.scanRate = 2000; % acquisition scan rate, determined by program

 
%% set up values for 'align to stim'
% *CHECK THESE EVERY TIME YOU RUN COMPILE BLOCKS*
% Ndnf vs. Vip project: 0.5, 2.5, 1.5, 1.5, 0.8, 0.7
% Peptide sensor (Christine) project: 0.5, **5***, 1.5, 1.5, 0.8, 0.7

% How many seconds of baseline?
constant.baseline_length = 0.5;


% ** Recommended to set this in the Info Sheet.
% How many seconds after stim should we look at?
% Can be overwritten by column in info
constant.after_stim = 2.5; %2.5

% Define (in seconds) where to look for the response peak?
constant.response_window = 1.5;

% define where to look for locomotor responses, in sec?
constant.locowindow = 1.5;

%noise floor of wheel
constant.locoThresh = 0.8;

% Define the neuropil coefficient
% This will be checked against Suite2p value and a warning will appear if they do not match
constant.neucoeff = 0.7;

