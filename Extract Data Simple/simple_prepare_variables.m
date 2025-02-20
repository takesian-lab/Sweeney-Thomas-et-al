function [v1, v2, v3, stim_length, blank_ind, trialsToIgnore, StimInfo] = simple_prepare_variables(block)
%% This function pulls out stim information from block.parameters
% Blank trials are coded as 0s within the specified parameter list

%% Get variables from parameters

v1 = block.parameters.variable1';
v2 = block.parameters.variable2';
if isfield(block.parameters, 'variable3')
    v3 = block.parameters.variable3';
else
    v3 = [];
end

stim_length = block.parameters.stimLength';

if isfield(block.parameters, 'trialsToIgnore')
    trialsToIgnore = block.parameters.trialsToIgnore;
end

%% Get StimCode (parameters and units) and blank trial type for all stim types

StimInfo = {};
StimInfo.StimProtocol = block.setup.stim_protocol;
 
switch block.setup.stim_protocol
    
    case 1 %Noiseburst
        StimInfo.StimType       = 'Noiseburst';
        StimInfo.Parameters     = {'Intensity', 'Duration'};
        StimInfo.Units          = {'dB', 'ms'};
        blank_ind               = v1 == 0;

    case 2 %RF
        StimInfo.StimType       = 'RF';
        StimInfo.Parameters     = {'Frequency', 'Intensity'};
        StimInfo.Units          = {'kHz', 'dB'};
        v2(v2 == -1000)         = 0;
        v1(v1 == 5.7)           = 5.6;
        blank_ind               = v2 == 0;

    case 3 %FM
        StimInfo.StimType       = 'FM';
        StimInfo.Parameters     = {'Speed', 'Direction'};
        StimInfo.Units          = {'octpersec', 'direction'};
        blank_ind               = v2 == 0;
        
        %Special for FM, use sweep direction as v2 and store intensity in V3
        %v2(v2 == 70) = 65; If we want to combine 65dB and 70dB stim
        v3 = v2;
        v2 = sign(v1);
        v1 = abs(v1); %Remove sign from v1
    
    case 5 %SAM
        StimInfo.StimType       = 'SAM';
        StimInfo.Parameters     = {'Rate', 'Depth'};
        StimInfo.Units          = {'Hz', '%'};
        blank_ind               = isnan(v1);
        
        %Special for SAM, make all stim with 0% AM depth have 0Hz modulation rate
        v2 = v2*100; %Change decimals to percent
        v1(v2 == 0) = 0;
    
    case 6 %SAMfreq
        StimInfo.StimType       = 'SAMfreq';
        StimInfo.Parameters     = {'Frequency', 'Depth'};
        StimInfo.Units          = {'kHz', '%'};
        blank_ind               = isnan(v1);
        
    case 7 %Frequency discrimination
        StimInfo.StimType       = 'FreqDisc';
        StimInfo.Parameters     = {'Frequency', 'Intensity'};
        StimInfo.Units          = {'delta_oct', 'dB'};
        blank_ind               = []; %No blank trials
        
        v3 = v1; %store frequeny values
        % convert frequencies to difference in octaves
        ov_v1 = log2(v1) - log2(block.TargetFreq);
        v1 = abs(ov_v1);
    
        
    case 8 % ABI (Modified 9/27/22 VAD)
        StimInfo.StimType       = 'ABI';
        StimInfo.Parameters     = {'Pulse Amplitude', 'Pulse Rate', 'Pulse-train Duration'};
        StimInfo.Units          = {'dBV', 'Hz', 'ms'};
        Cv1 = unique(v1);
        Cv2 = unique(v2);
        Cv3 = unique(v3);
        
        if length (Cv1) ==1 && length (Cv2) > 1 %Amplitude fixed
            v3_temp = v1; %storing V1 for tracing
            v1 = v2; %V1 is now Pulse Rate
            v2 = v3; %V2 is now Pulse-train Duration
            v3 = v3_temp; %V3 is now Pulse Amplitude
            StimInfo.Parameters     = {'Pulse Rate', 'Pulse-train Duration','Pulse Amplitude'};
            StimInfo.Units          = {'Hz', 'ms','dBV'};
            blank_ind               = []; %No blank trials
        else
            blank_ind               =  v1 == -1000; %[]; stim_v0
        end

    case 16 %ABI SAM
        StimInfo.StimType       = 'ABI_SAM';
        StimInfo.Parameters     = {'Mod rate', 'Modulation Depth'}; %to add pulse Rate
        StimInfo.Units          = {'kHz', 'depth'};
        blank_ind               = isnan(v1); %No blank trials
    case 19 %AMBehavior
        if strcmp(block.setup.stim_name, 'ABI')
        StimInfo.StimType       = 'ABI_AM';
        StimInfo.Parameters     = {'Mod rate', 'Modulation Depth', 'Current'}; %to add pulse Rate
        StimInfo.Units          = {'kHz', 'depth', 'mA'};
        else
        StimInfo.StimType       = 'Ac_AM';  
        StimInfo.Parameters     = {'Mod rate', 'Modulation Depth', 'Level'}; %to add pulse Rate
        StimInfo.Units          = {'kHz', 'depth', 'dBSPL'};
        end
        blank_ind               = isnan(v1); %No blank trials
        
        
    case 9 %H20
        StimInfo.StimType       = 'H2O';
        StimInfo.Parameters     = {'Solenoid', 'Outcome'};
        StimInfo.Units          = {'On/Off', 'Outcome'};
        blank_ind               = (~v1 + v2) == 2; %True shams only
        
        %Special for water: v1 = solenoid on/off, v2 = successful hit/sham trials, outcome = all behavioral outcomes
        %Toggle this option if you would rather look at all possible behavioral outcomes of random water
        useRandomWaterOutcomes = 0;
        if ~useRandomWaterOutcomes
            %Only keep hits and shams and add all other outcomes to trialsToIgnore
            trialsToIgnore(v2 == 0) = 1;
        else
            %Replace v2 with outcomes, trialsToIgnore will already contain uncategorized trials
            v2 = block.Outcome';
            v2(isnan(v2)) = 99; %These will be removed by trialsToIgnore
        end

        %Do not use stim_length, because solenoid was calibrated per experiment

    case 10 %Noiseburst ITI
        StimInfo.StimType       = 'NoiseITI';
        StimInfo.Parameters     = {'Intensity', 'Duration'};
        StimInfo.Units          = {'dB', 'ms'};
        blank_ind               = v1 == 0;
        
        %Special for noise, v2 is empty. Set v2 = stim_length
        v2 = stim_length;
        
    case 11 %Air
        StimInfo.StimType       = 'Air';
        StimInfo.Parameters     = {'Air', 'Air'};
        StimInfo.Units          = {'Air', 'Air'};
        blank_ind               = v1 == 0;
        
        %Special for air, v2 is all zeros. Set v2 = v1
        %Do not use stim_length, because solenoid was calibrated per experiment
        v2 = v1;
        
    case 111 %Air sham
        StimInfo.StimType       = 'Air Sham';
        StimInfo.Parameters     = {'Air Sham', 'Air Sham'};
        StimInfo.Units          = {'Air Sham', 'Air Sham'};
        blank_ind               = v1 == 0;
        
        %Special for air, v2 is all zeros. Set v2 = v1
        %Do not use stim_length, because solenoid was calibrated per experiment
        v2 = v1;
        
    case 12 %Spontaneous
        StimInfo.StimType       = 'Spontaneous';
        StimInfo.Parameters     = {'Sham', 'Sham'};
        StimInfo.Units          = {'Sham', 'Sham'};

        %Special for spontaneous, v1 and v2 are 0, make using stim_length
        [v1, v2] = deal(ones(size(stim_length)));
        
        %Special for spontaneous, label N random trials as "stim" and the
        %rest blank to evaluate false positive rate of simple_extract_data
        use_spont_as_control = 0;
        if use_spont_as_control
            N = 5;
            blank_ind = true(size(stim_length));
            stim_ind = 1:ceil(length(v1)/N):length(v1);
            blank_ind(stim_ind) = 0;
        else
            blank_ind = []; %No blank trials
        end

    case 13 %Maryse behavior stim
        StimInfo.StimType       = 'MaryseBehavior';
        StimInfo.Parameters     = {'Alternating tone', 'Repeating tone'};
        StimInfo.Units          = {'oct', 'oct'};
        blank_ind               = []; %No blank trials
    
        %Special for Maryse behavior stim, transform frequencies to octaves and bin
        bins = [0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1];
        
        v3 = v2; %Store kHz frequencies in v3 for reference
        
        repeating_log = unique(log2(v1)); %Find repeating frequency in octaves
        ov_v1 = log2(v1) - repeating_log;
        ov_v2 = abs(log2(v2) - repeating_log);
        
        %Find nearest bin for the alternating frequencies in octaves
        ov_v2_binned = zeros(size(ov_v2));
        for i = 1:length(ov_v2)
            if ov_v2(i) ~= 0
                [~, idx] = min(abs(bins - ov_v2(i)));
                ov_v2_binned(i) = bins(idx);
            end
        end
        
        %Swap v1 and v2
        v1 = ov_v2_binned;
        v2 = ov_v1;

    case 14 %Markpoints
        StimInfo.StimType       = 'Markpoints';
        StimInfo.Parameters     = {'Power', 'Duration'};
        StimInfo.Units          = {'', 'ms'};
        v1 = block.setup.XML.markPoints.Power;
        v2 = block.setup.XML.markPoints.Duration;
        blank_ind = (block.parameters.uncagingShutter == 0)';
                
    case 15 %Ripple
        StimInfo.StimType       = 'Ripple';
        StimInfo.Parameters     = {'Density', 'Intensity'};
        StimInfo.Units          = {'cyc per oct', 'dB'};
        blank_ind               = v2 == 0;
        
    case 17 %Ephys
        StimInfo.StimType       = 'Ephys';
        StimInfo.Parameters     = {'Rate', 'Pulse'};
        StimInfo.Units          = {'Hz', 'N'};
        blank_ind               = []; %No blank trials
                
    case 18 %Vocalizations
        StimInfo.StimType       = 'Vocalizations';
        StimInfo.Parameters     = {'Vocalization', 'Intensity'};
        StimInfo.Units          = {'#', 'dB'};
        blank_ind               = []; %No blank trials
        
    case 21 %ABCXYZ pattern
        StimInfo.StimType       = 'ABCXYZ_pattern';
        StimInfo.Parameters     = {'StimLabel', 'FMrates', 'Intensity'};
        StimInfo.Units          = {'triplet', 'oct per sec','dB'};
        blank_ind               = []; %No blank trials  
        
    case 211 %ABCXYZ random
        StimInfo.StimType       = 'ABCXYZ_random';
        StimInfo.Parameters     = {'StimLabel', 'FMrates', 'Intensity'};
        StimInfo.Units          = {'triplet', 'oct per sec','dB'};
        blank_ind               = []; %No blank trials
        
    case 50 %Loco (blocks with data aligned to loco bouts)
        StimInfo.StimType       = 'Loco';
        StimInfo.Parameters     = {'NA', 'BoutDuration'};
        StimInfo.Units          = {'NA', 'Seconds'};
        blank_ind               = []; %No blank trials

        %Use binned stim_length as V2 (corresponding to loco bout durations)
        v2_cats = 2:2:1000; %Make an arbitrarily  long set of bins spaced 2s apart
        v2 = nan(size(stim_length));
        for i = 1:length(stim_length)
            v2(i) = v2_cats(find(v2_cats >= stim_length(i),1,'first'));
        end


    case 1018 %Vocalizations_noise
        StimInfo.StimType       = 'VocNormFlip-Noise';
        StimInfo.Parameters     = {'vocORnoise', 'Intensity'};
        StimInfo.Units          = {'#', 'dB'};
        blank_ind               = v1 == 0 ;


    otherwise
        error('Add stim protocol to list.')
        
end
 
%% Standardize coding across all stim types
% Blanks will now be coded for as 0 in all stim parameters

%Make blank indices logical
if isempty(blank_ind)
    blank_ind = false(size(v1));
end

v1(blank_ind) = 0;
v2(blank_ind) = 0;
stim_length(blank_ind) = 0;

if ~isempty(v3)
    v3(blank_ind) = 0;
end


%Catch NaN stimuli to remove
if ~any(StimInfo.StimProtocol == [21 211]) % different structure to V1/V2
    if any(isnan(v1)) || any(isnan(v2)) || any(isnan(v3))
        error('NaNs found in variable list. Fix within simple_prepare_variables.')
    end
end

%% Remove trialsToIgnore from v1, v2, v3, blank_ind and stim_length
%This will be removed from the rest of the trials in makeTrialList

%Duplicate v1 and v2 so that we can pass "clean" trials to add_RF_mapping_to_StimInfo
v1b = v1;
v2b = v2;
    
if exist('trialsToIgnore','var')
    if any(trialsToIgnore)
        T = ~trialsToIgnore;
        v1b = v1b(T);
        v2b = v2b(T);
    end
else
    trialsToIgnore = zeros(size(v1));
end

%% Make RF map for StimInfo (for now, this only includes V1 and V2)

combined_stim = [v1b, v2b];
StimInfo = add_RF_mapping_to_StimInfo(StimInfo, combined_stim);
