% Generate synthetic data to test the script Ca_detection_spontaneous
% Contributors: Wisam, Maryse, Carolyn
%
% TODO:
% - Add noise to synthetic data (Wisam script: createRandomComplexTransient)
% - If no noise, synthetic data does not need to be low-pass filtered -> Make option to skip this in Ca_detection_spontaneous

%% Set parameters for synthetic data

nCells = 1; %N synthetic cells to create
fs = 30; %Framerate
T = 600; %Duration in seconds
useRandom = 0; %1 for random spike train, 0 for regular (TODO)
plotFigures = 0;

%% Create fake block with bare minimum data to run through Ca_detection_spontaneous

block = struct;
block.setup = struct;
block.setup.block_supname = 'SyntheticBlock';
block.setup.stim_protocol = 12; %spontaneous
block.setup.framerate = fs;
block.cell_number = [1:nCells]';
block.timestamp = [0:1/fs:T]'; %5 minutes at 30fps
block.F7 = nan(length(block.cell_number), length(block.timestamp)); %to be filled below
block.spks = zeros(length(block.cell_number), length(block.timestamp)); %we won't use this for anything yet
block.loco_activity = zeros(1,length(block.timestamp));
block.locomotion_trace = block.timestamp;
block.setup.constant.locoThresh = 0.7;
block.syntheticData = true;

%% Generate random synthetic spike trains to fill F7 with
% Each cell will have a different spikerate
% We will test how close the code is to guessing the ground truth spike rate

syntheticSpikeRates = 1./(randn(1,nCells) + 10); %[Hz]
measuredSpikeRates = nan(size(syntheticSpikeRates));

for c = 1:size(block.F7,1)
    disp(c)
    spikeRate = syntheticSpikeRates(c);
    
    %% Setup Constant Synthetic Spike Trains
    %spikeRate = 0.1;                             % [Hz]
    duration = T*1000;                            % [ms] coefficient is in seconds 
    timeVector = 0:duration;                      % [ms]

    calciumHalfTimeConstants = [250]; % 750 1000];              % [ms] Half time constants for GCaMP6f/m/s respectively
    calciumTimeConstants = calciumHalfTimeConstants/log(2); % [ms]
    calciumAmplitude = 0.3;

    spikeTrain = zeros(1,duration);                                             % Allocate a vector for the exponential kernel
    spikeTrain(round((1/spikeRate)*1000):round((1/spikeRate)*1000):length(spikeTrain)-1) = 1;	% set up constant spike train (hardcoded numbers s->ms)
    exponentialKernel = zeros(1, duration);                                     % Allocate a vector for the exponential kernel 

    % Create Synthetic Transients 
    % We will convolve the spike trains with exponential kernels
    for tau = calciumTimeConstants
        exponentialKernel = calciumAmplitude*exp(-timeVector(1:end-1)/tau);	% calculate the kernel
        convolvedSignal = conv(spikeTrain, exponentialKernel);              % the convolution
        convolvedSignal = convolvedSignal(1:duration);                      % take the first half, Length(Conv(A, B)) =  Length(A)+ Length(B) 
    end
    
    %% Sub-Sampling A Constant Synthetic GCaMP6s Transient

    msSR = 1000;                % [Hz] Synthetic Sample Rate 
    duration = T*msSR;         % [ms] coefficient is in seconds
    timeVector = 0:duration;    % [ms]

    % Calculate Downsample Rates
    experimentalSRs = [fs];%,15,10,5];                     % [Hz]
    [upsVec,dnsVec] = numden(sym(experimentalSRs/msSR));    % TODO: Currently symbolic needs conversion

    %% Resampling The synthetic transients at 30 Hz

    n = 5;      % n = 5 so that the antialiasing filter is of order 2*n*4 = 40
    beta = 20;	% Filter shape parameter

    ups = upsVec;
    dns = dnsVec;
    ups = 3;	% upsample factor
    dns = 100;	% downsample factor 

    dSR = (ups*msSR)/dns;   % new downsampled rate
    % downsampled time vector - Note: hard coded number (1000) ms->s
    dsTimeVector = 0:(1/dSR)*1000:duration;

    % Resample our transients (b is our filter)
    [resampledConvolvedSignal,b] = resample(convolvedSignal,ups,dns,n,beta);

    resampledConvolvedSignal = [resampledConvolvedSignal, 0]; %Add one extra frame lost during the resampling
    
    
    if plotFigures
        figure
        
        subplot(3,1,1)
        plot(timeVector(1:end-1)/1000, spikeTrain);
        set(gca, 'FontSize', 12)
        title('Random Synthetic Spikes')
        xlabel('Seconds') 
        ylabel('mV','Interpreter','tex') 

        subplot(3,1,2)
        plot(timeVector(1:end-1)/1000, convolvedSignal); 
        title({'Synthetic Transients'});
        xlabel('Seconds') 
        ylabel('\Delta F/F','Interpreter','tex')
        set(gca, 'FontSize', 12)
        legend('6f','Location','northwest','NumColumns',1)

        subplot(3,1,3)
        plot(dsTimeVector,resampledConvolvedSignal)
        title(['Resampled Synthetic GCaMP6s (Sampled Rate = ',num2str(dSR),' Hz)']);
        xlabel('Seconds') 
        ylabel('\Delta F/F','Interpreter','tex')
        set(gca, 'FontSize', 12)
        legend('6f','Location','northwest','NumColumns',1)
    end

    %% Fill F7 with synthetic spike train
    
    %Add some random baseline value (e.g. baseline fluorescence)
    resampledConvolvedSignal = resampledConvolvedSignal + randi([5,100], 1);
    
    block.F7(c,:) = resampledConvolvedSignal;
    
    %% Run through Ca_detection_spontaneous
    [~, data, ~, ~] = Ca_detection_spontaneous(block, c, 0, plotFigures);
    measuredSpikeRates(c) = data(4);
end

%% Figure

%Linear regression
mdl = fitlm(syntheticSpikeRates, measuredSpikeRates);

figure;
plot(mdl)
%scatter(syntheticSpikeRates, measuredSpikeRates)
title(['Synthetic vs. Measured transient rates - r2 = ' num2str(mdl.Rsquared.Ordinary)])
xlabel('Synthetic rates')
ylabel('Measured rates')

