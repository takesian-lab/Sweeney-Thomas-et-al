function [TrialData, StimInfo] = removeUnbalancedStim(TrialData, StimInfo, varargin)
%Unbalanced stim are stim that are not present in all blocks
%This script will remove them from TrialData and update StimInfo
%Alternatively, you can use varargin 'KeepOnlyStim' to specify a list of
%stim you would like to keep instead

%e.g.: For SAM stim keep only 0% AM noise: removeUnbalancedStim(TrialData, StimInfo, 'KeepOnlyStim', [0 0])

%% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% ADDITIONAL PARAMETERS
% Instead of removing unbalanced stim, use this script to keep only the list of stim you provide.
% List should be an N x 2 array with V1 in column 1 and V2 in column 2
addOptional(p, 'KeepOnlyStim', []);

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% Determine V1/V2 stim that are not present in every block

all_stim = StimInfo.Combined; %List of all possible stim in dataset

if isempty(ops.KeepOnlyStim)
    %Loop through all blocks in TrialData and record unique stim per block
    %Figure out which of all_stim are missing
    unique_stim = cell(size(TrialData));
    all_stim_present = zeros(size(all_stim,1),size(TrialData,2));
    for b = 1:size(TrialData,2)
        combined_stim = [TrialData(b).Stim_V1, TrialData(b).Stim_V2];
        unique_stim{b} = unique(combined_stim,'rows');

        for a = 1:size(all_stim,1)
            v1_present = unique_stim{b}(:,1) == all_stim(a,1);
            v2_present = unique_stim{b}(:,2) == all_stim(a,2);
            both_v1_v2 = (v1_present + v2_present) == 2;

            if sum(both_v1_v2)
                all_stim_present(a,b) = 1; 
            elseif sum(both_v1_v2) > 1
                error('This should not be > 1, troubleshoot')
            end
        end
    end

    unbalanced_stim_ind = sum(all_stim_present,2) < size(TrialData,2);
    unbalanced_stim = all_stim(unbalanced_stim_ind,:);

    %Create string result table to record
    results_mat = string(all_stim_present);
    results_mat = [StimInfo.V1_Label, StimInfo.V2_Label, results_mat];
    results_mat = [{'V1', 'V2', TrialData.Block}; results_mat];

    %Record and report results
    if ~isempty(unbalanced_stim)

        %Record information in StimInfo
        StimInfo.Original_StimInfo = StimInfo;
        StimInfo.OriginalStimPerBlock = results_mat;
        StimInfo.UnbalancedStimRemoved = unbalanced_stim;

        %Print unbalanced stim
        disp('Unbalanced stim found and removed:')
        disp(unbalanced_stim)

    else
        disp('No unbalanced stim found')
        return;
    end
end

%% Alternatively, identify all stim as unbalanced to keep only desired stim

if ~isempty(ops.KeepOnlyStim)
    unbalanced_stim = all_stim;
    
    %Remove the stim you want to keep from this list
    stim_to_remove = false(size(unbalanced_stim,1),1);
    
    for a = 1:size(ops.KeepOnlyStim,1)
        v1_present = unbalanced_stim(:,1) == ops.KeepOnlyStim(a,1);
        v2_present = unbalanced_stim(:,2) == ops.KeepOnlyStim(a,2);
        both_v1_v2 = (v1_present + v2_present) == 2;
        stim_to_remove(both_v1_v2) = 1;
    end
    %verify by printing: unbalanced_stim(stim_to_remove,:)

    unbalanced_stim(stim_to_remove,:) = [];
    
    %Record and report results
    if ~isempty(unbalanced_stim)

        %Record information in StimInfo
        StimInfo.Original_StimInfo = StimInfo;
        StimInfo.KeepOnlyStim = ops.KeepOnlyStim;
        StimInfo.UndesiredStimRemoved = unbalanced_stim;

        %Print unbalanced stim
        disp('Undesired stim found and removed:')
        disp(unbalanced_stim)

    else
        disp('No undesired stim found')
        return;
    end
end

%% Remove unbalanced (or undesired) stim from TrialData
% This method does not affect the blank trials

for b = 1:size(TrialData,2)
    %Determine which trials to remove
    combined_stim = [TrialData(b).Stim_V1, TrialData(b).Stim_V2];
    trials_to_remove = false(size(combined_stim,1),1);
    
    for a = 1:size(unbalanced_stim,1)
        v1_present = combined_stim(:,1) == unbalanced_stim(a,1);
        v2_present = combined_stim(:,2) == unbalanced_stim(a,2);
        both_v1_v2 = (v1_present + v2_present) == 2;
        trials_to_remove(both_v1_v2) = 1;
    end
    %verify by printing: combined_stim(trials_to_remove,:)
    
    %Remove from all types of TrialData
    TrialData(b).Stim_V1(trials_to_remove) = [];
    TrialData(b).Stim_V2(trials_to_remove) = [];
    if isfield(TrialData(b), 'Stim_V3')
        if ~isempty(TrialData(b).Stim_V3)
            TrialData(b).Stim_V3(trials_to_remove) = [];
        end
    end
    TrialData(b).Stim_Length(trials_to_remove) = [];
    TrialData(b).IsBlank(trials_to_remove) = [];
    TrialData(b).Sound_Time(trials_to_remove) = [];
    
    %BEHAVIORAL DATA
    TrialData(b).Loco(trials_to_remove,:) = [];
    TrialData(b).IsRunning(trials_to_remove) = [];
    if isfield(TrialData(b), 'Licks'); TrialData(b).Licks(trials_to_remove,:) = []; end
    if isfield(TrialData(b), 'Pupil'); TrialData(b).Pupil(trials_to_remove,:) = []; end
    if isfield(TrialData(b), 'Whisker'); TrialData(b).Whisker(trials_to_remove,:) = []; end
        
    %FUNCTIONAL DATA
    if isfield(TrialData(b), 'Trials'); TrialData(b).Trials(:,trials_to_remove,:) = []; end

    %PROTOCOL SPECIFIC VARIABLES
    switch StimInfo.StimProtocol
        
        case 13 %Maryse behavior stim
%             TrialData(b).StimLevel(trials_to_remove) = [];
            TrialData(b).Outcomes(trials_to_remove) = [];
            TrialData(b).ReactionTime(trials_to_remove) = [];
            TrialData(b).HoldingPeriod(trials_to_remove) = [];
            TrialData(b).WaitPeriod(trials_to_remove) = [];
            TrialData(b).WaitCondition(trials_to_remove) = [];
            TrialData(b).Result(trials_to_remove) = [];
            if isfield(TrialData(b), 'water_delivery'); TrialData(b).WaterDelivery(trials_to_remove) = []; end
    end
end

%% Update StimInfo

fields = {'Stim_V1', 'Stim_V2'};

concat_Combined = [];
for b = 1:size(TrialData,2)
    blank_ind = TrialData(b).IsBlank;
    temp_Combined = [];
    for f = 1:length(fields)
        if isfield(TrialData, fields{f})
            temp_Combined = [temp_Combined, TrialData(b).(fields{f})];
        end
    end
    temp_Combined(blank_ind,:) = []; %Don't include blank trials as part of variable list
    concat_Combined = [concat_Combined; temp_Combined];
end

%Final StimInfo will reflect data for all blocks combined
StimInfo = add_RF_mapping_to_StimInfo(StimInfo, concat_Combined);
TrialData = add_RF_mapping_to_TrialData(StimInfo, TrialData);

%% Check for problems

if isempty(StimInfo.V1) && isempty(StimInfo.V2)
    error('Removing unbalanced stim resulted in no trials left')
end
