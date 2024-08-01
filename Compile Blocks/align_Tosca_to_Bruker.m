function [block] = align_Tosca_to_Bruker(block, Bruker_trial_time)
% Align Tosca variables (Loco, Licks) to Bruker timestamp
% 
% Inputs:
% - block
% - Bruker_trial_time - time (s) of Tosca triggers detected on Bruker Voltage Recording
% 
% Updates:
% - Sound_Time
% - locomotion_time
% - activate_trials
% - concat licks

%% Preallocate and check for loco data

lick_timestamp = [];
lick_raster = [];
lick_BrukerTime = cell(1, length(Bruker_trial_time));
if isfield(block,'loco_data_raw')
    recordLoco = 1;
    locomotion_trace = [];
    loco_activity = [];
    loco_times = [];
    Loc_BrukerTime = cell(1, length(Bruker_trial_time));
else
    recordLoco = 0;
    warning('no locomotor data to align to bruker data')
end

%% Check that number of final Bruker trials matches number of Tosca trials
%If they don't match, run function to try to fix

%Record errors before fixing block with check_Bruker_triggers
originalErrors = block.errorData.error_trials;

nBrukerTrialsWithoutErrors = length(Bruker_trial_time) - length(originalErrors);
nToscaTrialsWithoutErrors = length(block.Tosca_times);

if nBrukerTrialsWithoutErrors ~= nToscaTrialsWithoutErrors
    %This function works by removing the extra Tosca triggers that don't match the Bruker voltage recording
    [block, Bruker_trial_time] = check_Bruker_triggers(Bruker_trial_time, block.errorData.New_sound_times, block);
    toscaTrialInd = block.errorData.corrected_ToscaTrialIndex;
else
    toscaTrialInd = 1:length(Bruker_trial_time);
end

%% find bruker times for entire tosca trial (do this for all trials,
% including errors, so that full loco and lick traces will not be missing data)

Bruker_trial_time_entire = cell(1, length(Bruker_trial_time));
for i = 1:length(Bruker_trial_time)
    ind = toscaTrialInd(i); %ind will be different from i when block is fixed only, needed to correctly access errorData
    t0 = block.errorData.Tosca_times{1,ind}(:) - block.errorData.Tosca_times{1,ind}(1,1);
    Bruker_trial_time_entire{1,i}(:) = t0 +  Bruker_trial_time(i);
end

% find the time (normalized from trial start) of sound time 
for i = 1:length(Bruker_trial_time)
    ind = toscaTrialInd(i); %ind will be different from i when block is fixed only, needed to correctly access errorData

    % put locomotor times on Bruker timescale
    if recordLoco == 1
        Loc_BrukerTime{i} = block.errorData.loc_Trial_times{ind}(:) + Bruker_trial_time(i);
        locomotion_trace = [locomotion_trace; Loc_BrukerTime{i}(:)]; % concat loco trials adjusted to Bruker time
        
        %remake concatenated activity in case block had to be fixed 
        if i == 1; loc_add = 0; else; loc_add = loco_times(end); end
        loco_times = [loco_times, block.errorData.loc_Trial_times{ind}' + loc_add];
        loco_activity = [loco_activity, block.errorData.loc_Trial_activity{ind}']; 
    end
    
    % put licks on Bruker timescale
    lick_BrukerTime{i} = (block.errorData.Tosca_times{ind} - block.errorData.start_time(ind)) + Bruker_trial_time(i);
    lick_timestamp = [lick_timestamp, lick_BrukerTime{i}]; % concat lick trials adjusted to Bruker time
    lick_raster = [lick_raster, block.errorData.lick_time{ind}]; %remake concatenated activity in case block had to be fixed
end

%% Go back and remove error trials

Bruker_trial_time(originalErrors,:) = [];
Bruker_trial_time_entire(originalErrors) = [];
if recordLoco == 1
    Loc_BrukerTime(originalErrors) = [];
end

%% Align sounds to Bruker timestamp
%Have to do this after removing errors

Sound_Time = Bruker_trial_time' + block.New_sound_times;

%% Save new variables in block
block.Sound_Time = Sound_Time;
block.BrukerTrialTimes = Bruker_trial_time_entire;
block.BrukerTrialStarts = Bruker_trial_time;
block.concat_times = lick_timestamp; %rewrite lick timestamps after bruker/fibphot alignment
block.concat_licks = lick_raster; %rewrite in case block had to be fixed
if recordLoco == 1
    block.locomotion_trace = locomotion_trace;
    block.loco_activity = loco_activity; %rewrite in case block had to be fixed
    block.loco_times = loco_times; %rewrite in case block had to be fixed
end
