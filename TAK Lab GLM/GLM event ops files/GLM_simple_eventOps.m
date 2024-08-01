% set ops for TAK LAB GLM
% V2 by Carolyn Sweeney 1/30/23


%-----------------------------------------
% Set up Paths and set the stim Type
%-----------------------------------------
% data_path = 'Z:\Christine\Analysis_ByProject\8_Behavior_AsscLearning_VIPsensor\extractedData\all_VIPsensor';
% %         blocks_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Final info sheets and compiled blocks\Blocks as of May 2022';
% save_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\GLM\with bootstrap\FM';
% matchFile_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Final info sheets and compiled blocks\Info';
% matchFile = 'Combined Info more mice.xlsx';
% Ops.stimTypes = {'FreqDisc'}; %,'RF','NoiseITI'}; %'SAM','SAMfreq',
% NetworkAnalysisPath = ' Z:\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\Network Analysis_NoiseCorr\Backup';

data_path = 'Z:\Christine\Analysis_ByProject\8_Behavior_AsscLearning_VIPsensor\extractedData\all_VIPsensor';
%         blocks_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Final info sheets and compiled blocks\Blocks as of May 2022';
save_path = 'Z:\Carolyn\Fiber Photometry\Behavior\GLM_VIP mice';
matchFile_path = '';
matchFile = '';
Ops.stimTypes = {'FreqDisc'}; %,'RF','NoiseITI'}; %'SAM','SAMfreq',
NetworkAnalysisPath = '';
%-----------------------------------------
% Settings for the Ridge Regression
%-----------------------------------------
Ops.window = [-0.5 5];          %   shift (in seconds) around events
Ops.dataChunk = 7;              %   time (in sec) to use as a event/trial for training - This is used if Ops.ByTrials = 0
Ops.ByTrials = 0;               %   Use just the trial-based data. 
Ops.testLambda = 2.^(0:15);     %   lambda values to test for fitting
Ops.RepeatFits = 10;            %   number of times to test the lambda fit
Ops.keepContig = 0;             %   for troubleshooting. Keep training data in large, contiguous segments
Ops.numBoots = 10;               %   number of times to boostrap data
Ops.threshFits = 0.25;          %   threshold for visualizing correlations - as of now, this is an arbitrary number. Future iterations may use bootstrap values to set value
Ops.FracTest = 0.2;             %   fraction of data to leave out for testing (and within the validation)
Ops.FibPhotChannel = 2;         %   If running FibPhot data, which channel to run the GLM on (1 = blue, 2 = green, 3 = red);

%-----------------------------------------
% Settings for the Regressors
%-----------------------------------------
Ops.activeCells = 0;            %   Only run on cells deemed active in ExtractData
Ops.MultipleActiveStim = 0;     %   if Ops.activeCells, Only use neurons that are responsive to multiple stim 
Ops.stimnumthesh = 1;           %   if Ops.MultipleActiveStim, the number of stim required to include the cell


Ops.useSounds = 0;              %   use sounds - dont sort subtypes (blanks are removed)
Ops.useV1 = 0;                  %   block.parameters.variable 1, sort out blanks
Ops.useV2 = 1;                  %   block.parameters.variable2, sort out blanks
Ops.useMag = 0;                 %   regressors are by identity not binary response (i.e. analog sound events)
Ops.useAllCells = 0;            %   use all neurons in FOV as regressors, not averaged together.
Ops.useCorrCell = 0;            %   cells with highest correlated activity used as regressors 
Ops.Corrstim = 'ExceptCurrent'  %   if useCorrCell, pick stim type to correlate with current data set. 'ExceptCurrent' will do all stim types except for the one being done in this experiment
Ops.CorrType = 'Noise'          %   if Ops.Corrstim, type of correlation to look at. 
Ops.CorrMean = 0;               %   if Ops.Corrstim, average responses of correlated neurons for single regressor (split into positively and negatively correlated neurons)
Ops.Residuals = 0;              %   Use residuals as regressor - this will use all cells in FOV. If you want to restrict to correlated neurons, choose Ops.useCorrCell && Ops.CorrResid 
Ops.CorrResid = 0;              %   if Ops.Corrstim use only the residuals from neurons that have correlations in different stim block
Ops.useOutcome = 1;             %   for Go/No-go data
Ops.useNeighbor = 0;            %   cells that are located closest to eachother within FOV
Ops.localCellPer = 0;           %   if Ops.useNeighbor, what percent of local cells to use
Ops.useNeighborMean = 0;        %   if Ops.useNeighbor, average neighbor activity

Ops.useLocomotor = 1;           %   use locomotor trace as regressor
Ops.LocoBin = 0;                %   make locomotor data binary instead of a velocity
Ops.NormLoco = 1;               %   normalize the locomotor traces
Ops.LocoOnset = 0;              %   Use locomotor onset rather than whole trace.
Ops.Acceleration = 0;           %   take derivative of loco activity
Ops.ShiftLoco = 1;              %   Apply temporal shift to locomotor activity, similar to what we do with sound (defaults to no shift);

Ops.usePupil = 1;               %   use pupil trace as regressor
Ops.NormPupil = 1;              %   normalize the pupil traces
Ops.ShiftPupil = 1;             %   Apply temporal shift to pupil activity, similar to what we do with sound;

Ops.useLicks = 1;               %   use Licks as regressor (behavior blocks)
Ops.ShiftLicks = 1;             %   Apply temporal shift to lick activity, similar to what we do with sound;

%-----------------------------------------
% Do you want to predict something other than Calicum/neuronal activity?
%-----------------------------------------
Ops.PredictMotor = 0;           % use set regressors to predict motor transient. 


display(Ops)
% end