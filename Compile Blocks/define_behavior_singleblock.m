function [block] = define_behavior_singleblock(block)
% This function accesses the stim and locomotor data from the Tosca folder
% and stores it in block
%
% Argument(s):
%   block (struct)
%
% Returns:
%   block (struct)
%
% VERSION HISTORY
%  - V1 March 2020
%  - V2 = current version adds pupil data

%% Skip this function if Tosca data is not available

if ismissing(block.setup.Tosca_path)
    disp('Skipping Tosca data...');
    return
end

disp('Pulling out Tosca data...');

%% Go to Tosca folder and pull out files related to setup.Tosca_run

setup = block.setup;
cd(setup.Tosca_path{1})
Tosca_Run_number = num2str(setup.Tosca_run);
Tosca_Session = num2str(setup.Tosca_session);
mouseID = char(setup.mousename);

%Check whether Tosca files match the information given in info
allfiles=dir('*Run*');
desiredFilename = [mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number];
if isempty(find(contains({allfiles(:).name}, desiredFilename),1))
    %Accommodate datasets where Tosca path mousename was saved with/without dashes 
    temp_mouseID = char(add_or_remove_dashes_from_mousename(setup.mousename));
    desiredFilename = [temp_mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number];
    
    %If still no files found, user made a mistake
    if isempty(find(contains({allfiles(:).name}, desiredFilename),1))
        error('Could not find Tosca files. Check that mouseID, Session, and Run are correct.')
    end
end

%% NEW: generate TL from tosca log....
FN = [mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.txt'];
TL = tosca_create_log(FN,'skip_avi',block.setup.skip_avi);
Params = TL.params;
setup.Tosca_date = datestr(Params.Info.Date); %Store the date the Tosca file was created (to compare with Bruker data)

%% Check if locomotor data were recorded

if isfield(Params,'Tosca') == true % newer tosca structure
    if Params.Tosca.Connect.Loco ==1
        recordLoco = 1;
    else
        recordLoco = 0;
        warning('no locomotor data recorded for the run')
    end
    
elseif isfield(Params,'DAQ') == true % older tosca structure
    if Params.DAQ.Loco == 1
        recordLoco = 1;
    else
        recordLoco = 0;
        warning('no locomotor data recorded for the run')
    end

else
    error('Tosca Params have a different configuration to detct loco. Ask 2P code slack for help')
end

%% check if Pupil data were recorded

if isfield(Params,'Tosca') == true % newer tosca structure
%     if Params.Tosca.Connect.Pupil == 1 && info_video == 1
    if Params.Tosca.Connect.Pupil == 1
        recordPupil = 1;
    else
        recordPupil = 0;
    end
    
elseif isfield(Params,'DAQ') == true % older tosca structure
    if strcmp(Params.DAQ.Pupil, 'OFF') % || info_video == 0
        recordPupil = 0;
    else
        recordPupil = 1;
    end
else
    error('Tosca Params have a different configuration to detect Pupil. Ask 2P code slack for help')
end

%% Read data from the run
Data = TL.trials;
ntr = length(Data);% number of trials

if recordLoco ==1
    loco_data = tosca_read_loco([mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt']); %locomotor data
    [loco_data, ignoreLoco] = check_loco_data(block, loco_data, block.setup.redo_loco_check); %Check loco data to make sure wheel was working
    loco_trace_times = [];
    loco_trace_activity =[];
end

%% loop through trials to get all the relevent data

%Preallocate variables: 
nanmat = nan(1,ntr);
emptycell = cell(1,ntr);

targetFreq = nan;
[result, result_orig, group] = deal(strings(size(nanmat)));
[start_time, end_time, rxn_time, StateChange, New_sound_times, New_sound_idx,New_pulse_times, New_pulse_idx,pulseDuration_sec,water_delivery, optoOn, target_state_duration] = deal(nanmat);
[Tosca_times, zero_times, zero_loc, activity_trial, licks, states, b_Outcome, trialType] = deal(emptycell);
concatenated_trial_times = [];
concatenated_lick_times = [];

if any(setup.stim_protocol == [13 19]) %Freq/AM detection behavior
    [holdingPeriod, waitPeriod, waitCondition] = deal(nanmat);
elseif block.setup.stim_protocol == 33 % Fear conditioning; FM sweep
    stim_Reps = Params.Tosca.Flowchart(2).State.Timeouts.Expr;
    % --- CJL to add stim reps calculation in case the 'state3' has
    % additional stim reps (aka the air puff delivered on the last rep) 
    
end


%Loop through trials        
for t=1:ntr
    
    s = tosca_read_trial(Params,Data,t); %read_trial gives us more info than read_run alone
    if isempty(s)
        continue
    end
    
    Tosca_times{t}  = s.Time_s; %pulls out the tosca generated timestamps for each trial
    start_time(t)   = Tosca_times{1,t}(1,1);
    end_time(t)     = Tosca_times{1,t}(1,end);
    zero_times{t}   = Tosca_times{1,t}(1,:)-start_time(t); %set the start of each trial to zero and normalize
    licks{t}        = s.Lickometer;
    states{t}       = [0 (diff(s.State_Change)>0)];
    if isfield(s, 'Rxn_time_ms'); rxn_time(t) = s.Rxn_time_ms; end
    if isfield(s, 'Group'); group{t} = char(string(s.Group)); end
    
    %Each Tosca trial can have several state changes (e.g. Start, Sound, Wait, etc)
    %Figure out which state change corresponds to the sound
    allStateChanges = [1, find(states{t})]; %Index of each state change. Concatenate 1 to represent trial start so that allStateChanges will match size of statenames. MET 9/17/22
    statenames = string({s.states.name});
    result(t) = s.Result; %Result that may be corrected
    result_orig(t) = s.Result; %Save original result in block
    
    %Recover error trials that made it past the target state
    %These will be recorded as Hit;Error, Miss;Error, etc.
    %Remove the ';Error' portion and keep as regular trials
    charResult = char(result(t));
    if length(charResult) > 6 %6 is the number of characters in ;Error
        if isequal(charResult(end-5:end), ';Error')
            result(t) = string(charResult(1:end-6));
        end
    end
    
    if strcmp(result(t), 'suspect di')
        result(t) = 'Error';
    end
    
    if any(setup.stim_protocol == [7])
          if strcmp(result(t), 'No trial')
               result(t) = 'Error';
          end
    end
                
    %Freq/AM detection behavior
    if any(setup.stim_protocol == [13 19]) 
        
        %Catch other types of error trials
        if strcmp(result(t), 'No trial')
           if any(strcmp(s.History, 'Abort'))
               result(t) = 'Error';
           end
        end
    
        if ~strcmp(result(t), 'Error')
            %Find duration of holding period (fixed duration for each trial)
            holdingPeriod(t) = s.Script.output; %should be equivalent to allStateChangesInS(2) - allStateChangesInS(1);

            %Added variable wait condition to paradigm on 9/14/21
            %Figure out how long the holding period was extended if the mouse didn't meet the wait condition
            HoldStateIdx = find(statenames == 'Hold', 1, 'first');
            waitPeriod(t) = 0;
            if ~isempty(HoldStateIdx) %Before I introduced the wait condition
                waitPeriod(t) = s.states(HoldStateIdx).duration;

                %Record wait condition
                if isfield(s, 'Events')
                    waitNames = {'Withold_2_minus7s', 'Withold_4_minus6s', 'Withold_1_minus3s'};
                    for w = 1:length(waitNames)
                        if isfield(s.Events, waitNames{w})
                            waitCondition(t) = s.Events.(waitNames{w}).Lickometer.Low;
                        end
                    end
                elseif isfield(s, 'Withold_2s') || isfield(s, 'Withold')
                    waitCondition(t) = 2;
                else
                    error('Look for wait condition variable')
                end
            end
        end
    end
    
    %Finish finding where the sound happened
    %NEW 10/10/23 allow error trials to go through this loop. If the error
    %happened after the sound, we will still be able to record the SoundTime for it -MET
    
    %if ~strcmp(result(t), 'Error')
        %Find where the 'Sound' or 'cue' is (usually the second state in a trial, after 'Wait/Nothing')
        stateOI = find(statenames == 'cue' | statenames == 'Sound' | statenames == 'Standard' | statenames == 'Hold' | statenames == 'End' | statenames == 'ABI' | statenames == 'Sound1', 1, 'first');

        %Record state change representing the sound
        if ~isempty(stateOI)
            StateChange(t) = allStateChanges(stateOI); %MET: 9/17/22 removed -1, see allStateChanges above            
            New_sound_times(t) = zero_times{1,t}(1,StateChange(t));
            New_sound_idx(t) = StateChange(t);
            target_state_duration(t) = s.states(stateOI).duration;   
            if isfield(Params, 'Tosca')
            stateChannelStruct = Params.Tosca.Flowchart(stateOI).State.SigMan.Channels; % save for later finding the channel delays    
             else
                stateChannelStruct = [];
            end
        elseif strcmp(result(t), 'Error') || length(allStateChanges) == 1
            %This means this trial stopped and did not get to the Sound state (e.g. could be an 'Abort' trial)
            result(t) = 'Error'; %Classify as error if it wasn't already
            %NEW 10/10/23 Record info for first state change, that way we can still have approximate trial times for errors
            StateChange(t) = allStateChanges(1);
            New_sound_times(t) = zero_times{1,t}(1,StateChange(t));
            New_sound_idx(t) = StateChange(t);
            target_state_duration(t) = s.states(1).duration;
            if isfield(Params, 'Tosca')
            stateChannelStruct = Params.Tosca.Flowchart(1).State.SigMan.Channels; % save for later finding the channel delays    
            else
                stateChannelStruct = [];
            end
        else
            %This means the trial had more than one state, but none of them were the sound
            %(might need to add a new state to the list of statenames above)
            error('stateOI not found. Ask 2P slack for help')
        end
    %end
    
    %Check for Opto during cue -- CJL added 4/24/2024

    if isfield(s.states(stateOI),'Pulse') % check if there is pulse during cue period for optogenetic activation
        optoOn(t) = s.states(stateOI).Pulse.Level.Volts>0; %check if opto is on - greater than 0 Volts
    end
    
    % --- check for any delays in the target state --- % CJL added 4/24/2024
    % catch any trials where the start of target state is not when the sound is played 

    if ~isempty(stateOI) % AT added 5/11/2024
        stateChannelStruct = Params.Tosca.Flowchart(stateOI).State.SigMan.Channels; 
    else
        stateChannelStruct = [];
    end
    

    if ~isempty(stateChannelStruct)
        SigChan = find(contains({stateChannelStruct.Name},'Signal')); % channel index for signal
        PulseChan = find(contains({stateChannelStruct.Name},'Pulse')); % channel index for pulse
        if ~isempty(SigChan)
            SigDelay_sec = stateChannelStruct(SigChan).Gate.Delay_ms./1000;
            New_sound_times(t)= zero_times{1,t}(1,StateChange(t)) + SigDelay_sec; % in seconds wrt to start of trail
            realSoundTime = s.Time_s(StateChange(t))+SigDelay_sec;% find the sound time with respect to the Time_s
            % find the index of when this sound occurs
            [~,New_sound_idx(t)] = min(abs(s.Time_s-realSoundTime)); % smallest difference in sec between the real sound time and all the times of this trial
        end
        
        if ~isempty(PulseChan) % check for second channel in target state (pulse) as well
            PulseDelay_sec = stateChannelStruct(PulseChan).Gate.Delay_ms./1000;% in seconds
            New_pulse_times(t)= zero_times{1,t}(1,StateChange(t)) + PulseDelay_sec; % pulse time wrt start of trial
            realPulseTime = s.Time_s(StateChange(t))+PulseDelay_sec;
            [~,New_pulse_idx(t)] = min(abs(s.Time_s-realPulseTime)); % smallest difference in sec between the real pulse time and all the times of this trial
            
            % record pulse length if opto is on
            if optoOn(t) == 1
                pulseDuration_sec(t) = stateChannelStruct(PulseChan).Gate.Duration_ms./1000;
            else
                pulseDuration_sec(t) = 0;
            end
        end
    end
 
    %Get CS+/CS- results
    try
        switch result(t)
            case 'Hit'
                b_Outcome{t}=1;
                trialType{t}=1;
                water_delivery(t) = 1; %Set water delivery to 1 for all hits, then update if block has catch trials
                
                %catch trials
                if isfield(s, 'cue') %Carolyn behavior
                    if isfield(s.cue, 'lick')
                        water_delivery(t) = s.cue.lick.H2O.Duration_ms>0;
                    end
                elseif isfield(s, 'Sound') %Maryse behavior
                    if isfield(s.Sound, 'Lick')
                        water_delivery(t) = s.Sound.Lick.H2O.Duration_ms>0;
                    end
                end
                                
                if setup.stim_protocol == 7
                    if isfield(s.cue.Signal,'Waveform')
                        targetFreq = s.cue.Signal.Waveform.Frequency_kHz; %pull out the target frequency
                    else
                        targetFreq = s.cue.Signal.Tone.Frequency_kHz; %pull out the target frequency
                    end
                end

            case 'Miss'
                b_Outcome{t}=0;
                trialType{t}=1;
                if setup.stim_protocol == 7
                    if isfield(s.cue.Signal,'Waveform')
                        targetFreq = s.cue.Signal.Waveform.Frequency_kHz;%pull out the target frequency
                    else
                        targetFreq = s.cue.Signal.Tone.Frequency_kHz;%pull out the target frequency
                    end
                elseif setup.stim_protocol == 9
                    try
                        targetFreq = s.cue.Signal.FMSweep.Rate_oct_s;
                    catch
                        targetFreq = nan;
                    end
                end

            case 'Withhold'
                b_Outcome{t}=3;
                trialType{t}=0;

            case 'False Alarm'
                b_Outcome{t}=4;
                trialType{t}=0;

            otherwise
                b_Outcome{t}=NaN;
                trialType{t}=NaN;
        end
    catch
        warning('No CS+/CS- trial type found (classifying as Sound Pav)')
        b_Outcome{t} = s.Result;
        trialType{t} = ('SoundPav');
    end
    
    % now that all the loco times are corrected per trial, put them
    % back together to get a loc trace that is on a correct timescale.Also,
    % do the same for licks and trial times.
    if t == 1
        trial_add = zero_times{t};
    else
        trial_add = zero_times{t} + trial_add(end);
    end
    concatenated_trial_times = [concatenated_trial_times; trial_add'];
    concatenated_lick_times = [concatenated_lick_times; licks{t}'];
    
    if length(concatenated_trial_times) ~= length(concatenated_lick_times)
        error('Traces should be the same length')
    end
    
    if recordLoco == 1
        zero_loc{t} = s.loco.t - start_time(t); %each trial's locomotor trials, corrected by zeroing out the start of each trial.
        activity_trial{t} = abs(s.loco.speed); %divide the loco activity by trials to use in define_sound_singleblock
        
        %Replace loco with NaNs if problem was found with loco
        if ignoreLoco
            activity_trial{t}(:) = NaN;
        end
        
        if t == 1
            loc_add = zero_loc{t};
        else
            loc_add = zero_loc{t} + loc_add(end);
        end
        loco_trace_times = [loco_trace_times; loc_add];
        loco_trace_activity = [loco_trace_activity; activity_trial{t}];
        
        if length(loco_trace_times) ~= length(loco_trace_activity)
            error('Traces should be the same length')
        end
    end
      
    % Pupillometry/camera frames/times
    if recordPupil == 1
        FrameData = Data{1,t}.states;
        if isfield(FrameData, 'frames')
            startpupil = FrameData(1).tframe(1);
            [framenum, frametime] = pupil_trial_frames(FrameData,startpupil);
            PupilFrameData{t}.frames = framenum;
            PupilFrameData{t}.frame_times = frametime;
            PupilFrameData{t}.AVIname = Data{t}.aviFile;

        else % condition created for when camera recording turns off spontaneously
            warning(['no frame data for trial ' num2str(t)])
            
            PupilFrameData{t}.frames = NaN;
            PupilFrameData{t}.frame_times = NaN;
            PupilFrameData{t}.AVIname = NaN;
        end
    end
end

%% Extract stimulus-specific variables

[V1, V2, V3, stimLength] = define_stim_parameters(Data, Params, setup);

%% ABI suppl. params

if setup.stim_protocol == 8 %ABI
    pulse_width = zeros(1,length(Data));
    for m = 1:length(Data)
        if isfield(Params, 'Output_States')
            pulse_width(1,m)= Params.Output_States(2).StimChans.Waveform.Pulse_Train.Pulse_width_us(1);
%             pulse_width(1,m)= Params.Output_States(2).StimChans.Stimulus.Waveform.Pulse_Train.Pulse_width_us(1); %inventory purpose (Tosca version 2366)
        elseif isfield(Params, 'Tosca')
            %         pulse_width(1,m)= Params.Tosca.Flowchart(3).State.SigMan.Channels.Channel.Waveform.Width_us(1); %store for inventory purpose %Added Channel. (VAD 12/12/22)
            pulse_width(1,m)= Params.Tosca.Flowchart(3).State.SigMan.Channels.Waveform.Width_us(1); %store for inventory purpose
        end
    end
end

if setup.stim_protocol == 16 %ABI-SAM
    pulse_amplitude = zeros(1,length(Data));
    pulse_width = zeros(1,length(Data));
    for m = 1:length(Data)
        try
            if isfield(Params, 'Output_States')
%             pulse_amplitude(1,m) = Data{m}.cue.Current.Level.dB_Vrms;
%             pulse_width(1,m)= Params.Output_States(2).StimChans.Waveform.Pulse_Train.Pulse_width_us(1);
            pulse_amplitude(1,m) = Data{m}.cue.Current.Level.dB_re_1_Vrms; %tosca version 2366
            pulse_width(1,m)= Params.Output_States(2).StimChans.Waveform.Pulse_Train.Pulse_width_us(1);
%             pulse_width(1,m)= Params.Output_States(2).StimChans.Stimulus.Waveform.Pulse_Train.Pulse_width_us(1);%tosca version 2366
            elseif isfield(Params, 'Tosca')
                %             pulse_amplitude(1,m) = Data{m}.cue.Current.Level.dB_Vrms;
                pulse_amplitude(1,m) = Data{m}.cue.Current.Level.dB_re_1_Vrms; %tosca version 2366
                %             pulse_width(1,m)= Params.Tosca.Flowchart(3).State.SigMan.Channels.Channel.Waveform.Width_us(1); %Added Channel. (VAD 12/12/22)
                pulse_width(1,m)= Params.Tosca.Flowchart(3).State.SigMan.Channels.Waveform.Width_us(1);
            end
        catch
            %Blank Trial
            stimLength(1,m) = 0;
            pulse_amplitude(1,m) = 0;
            pulse_width(1,m)= 0;
        end
    end  
end

%% Freq/AM detection behavior

if any(setup.stim_protocol == [13 19]) 
    targetFreq = unique(V1);
    
    %StimLength:  Total stim duration = holding period + wait condition + response window
    responseWindow = Params.Script.Fixed.ResponseWindow_s; %Response window is the same for all trials
    stimLength = holdingPeriod + waitPeriod + responseWindow;
    
    %Correct reaction time (this will get reaction times up to responseWindow instead of responseWindow - 0.2s):
    rxn_time = tosca_find_reaction_times(Params, Data, 'Lick', responseWindow);
    rxn_time(isnan(rxn_time)) = -1;
end
%% check for complete triplets in rhythmic behavior.
if any(setup.stim_protocol == [21]) 
    result = check_complete_triplet(V1,result);
end
%% Check for tosca trials that are errors, and remove them from the data

% save all of the trials before removing errors
error_trials = find(strcmp(result, 'Error'));
block.errorData.error_trials = error_trials;
block.errorData.result = result;
block.errorData.result_orig = result_orig;
block.errorData.start_time = start_time;
block.errorData.Tosca_times = Tosca_times;
block.errorData.New_sound_times = New_sound_times;
block.errorData.New_sound_idx = New_sound_idx;
block.errorData.New_pulse_times = New_pulse_times;
block.errorData.New_pulse_idx = New_pulse_idx;
block.errorData.target_state_duration = target_state_duration;
block.errorData.Group = group;
block.errorData.lick_time = licks;
if recordLoco
    block.errorData.loc_Trial_times =  zero_loc;
    block.errorData.loc_Trial_activity =  activity_trial;
end
if recordPupil
    block.errorData.PupilFrameData = PupilFrameData;
end

% remove errors
k = error_trials;
k(k > ntr) = []; %Cannot remove trials that don't exist

if ~isempty(k)
    start_time(:,k) = [];
    Tosca_times(:,k)  = [];
    New_sound_times(:,k) = [];
    New_sound_idx(:,k) = [];
    New_pulse_times(:,k) = [];
    New_pulse_idx(:,k) = [];
    pulseDuration_sec(:,k) = [];
    licks(:,k) = [];
    b_Outcome(:,k) = [];
    trialType(:,k) = [];
    rxn_time(:,k) = [];
    stimLength(:,k) = [];
    water_delivery(:,k) = [];
    optoOn(:,k) = [];
    
    if recordLoco
        zero_loc(:,k) = [];
        activity_trial(:,k) = [];
    end
    
    %If not noise or spontaneous
    if ~any(setup.stim_protocol == [1 12])
        V1(:,k)=[];
        V2(:,k)=[];
        if ~isempty(V3)
            V3(:,k)=[];
        end
    end
    
    %Freq/AM detection behavior
    if any(setup.stim_protocol == [13 19]) 
        holdingPeriod(:,k) = [];
        waitPeriod(:,k) = [];
        waitCondition(:,k) = [];
    end
    
    if recordPupil == 1
        flpk = flip(k);
        for d = 1:length(flpk)
            PupilFrameData(flpk(d)) = [];
        end
    end
end

%% Save everything to block

block.start_time = start_time;
block.Tosca_times = Tosca_times;
block.New_sound_times = New_sound_times;
block.New_sound_idx = New_sound_idx;
block.New_pulse_times = New_pulse_times;
block.New_pulse_idx = New_pulse_idx;
block.lick_time = licks;
block.concat_times = concatenated_trial_times;
block.concat_licks = concatenated_lick_times;
block.Outcome = cell2mat(b_Outcome);
block.trialType = cell2mat(trialType);
block.TargetFreq = targetFreq;
block.parameters.variable1 = V1; %index of variable 1 (e.g. frequency)
block.parameters.variable2 = V2; %index of variable 2 (e.g. level)
block.parameters.variable3 = V3; %index of variable 3 (e.g. level)
block.parameters.stimLength = stimLength;
block.parameters.optoPulseLength = pulseDuration_sec;
block.parameters.optoOn = optoOn;
block.rxn_time = rxn_time;
block.setup = setup;
block.water_delivery = water_delivery;


%Stim specific variables
switch setup.stim_protocol
    case 8 %ABI
        block.parameters.pulse_width = pulse_width; 

    case 16 %ABI SAM 
        block.parameters.pulse_width = pulse_width; 
        block.parameters.pulse_amplitude = pulse_amplitude; 
    
    case {13 19} %Freq/AM detection
        block.holdingPeriod = holdingPeriod;
        block.waitPeriod = waitPeriod;
        block.waitCondition = waitCondition;
        block.stim_level = unique(V3);

    case 33 %Fear conditioning
        block.stim_Reps = stim_Reps;
    case 7 % frequency discrimination:
         [block] = FreqDisc_Behavior_singleblock(block);
end

if recordLoco == 1
    block.loco_data_raw = loco_data; %raw loco data
    block.loco_activity = loco_trace_activity; %trace of velocity
    block.loco_times = loco_trace_times; %time-corrected, should match the velocity trace
    block.loc_Trial_times = zero_loc; %timestamps for loco by each trial
    block.loc_Trial_activity = activity_trial; % time-corrected velocity for each trial
end

if recordPupil == 1
    block.PupilFrameData = PupilFrameData;
end
