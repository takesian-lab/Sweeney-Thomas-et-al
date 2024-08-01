function block = helper_removeTrialsFromBlock(original_block, trialsToRemove)
%Remove specific trials from block

block = original_block;
block.errorData.ManuallyRemovedTrials = trialsToRemove; %Store record of manual trial removal

%Remove trials from all variables

%TOSCA
block.start_time(trialsToRemove) = [];
block.Tosca_times(trialsToRemove) = [];
block.New_sound_times(trialsToRemove) = [];
block.New_sound_idx(trialsToRemove) = [];
block.lick_time(trialsToRemove) = [];
block.Outcome(trialsToRemove) = [];
block.trialType(trialsToRemove) = [];
block.rxn_time(trialsToRemove) = [];
block.water_delivery(trialsToRemove) = [];
if isfield(block,'loc_Trial_times')
    block.loc_Trial_times(trialsToRemove) = [];
    block.loc_Trial_activity(trialsToRemove) = [];
end
paramfields = fields(block.parameters);
for f = 1:length(paramfields)
    if ~isempty(block.parameters.(paramfields{f}))
        if length(block.parameters.(paramfields{f})) > 1
            block.parameters.(paramfields{f})(trialsToRemove) = [];
        end
    end
end
block.Sound_Time(trialsToRemove) = [];
block.active_trials(trialsToRemove) = [];
block.aligned_stim.velocity(trialsToRemove,:) = [];
block.aligned_stim.licks(trialsToRemove,:) = [];
%TODO: add pupil and whisker

%BRUKER
if isfield(block, 'F')
    block.BrukerTrialTimes(trialsToRemove) = [];
    block.BrukerTrialStarts(trialsToRemove) = [];
    block.aligned_stim.F_stim(:,trialsToRemove,:) = [];
    block.aligned_stim.Fneu_stim(:,trialsToRemove,:) = [];
    block.aligned_stim.F7_stim(:,trialsToRemove,:) = [];
    block.aligned_stim.spks_stim(:,trialsToRemove,:) = [];
    block.aligned_stim.df_f(:,trialsToRemove,:) = [];
    block.aligned_stim.zscore(:,trialsToRemove,:) = [];
    
    %Chan2
    if isfield(block.aligned_stim, 'F_stim_chan2')
        block.aligned_stim.F_stim_chan2(:,trialsToRemove,:) = [];
        block.aligned_stim.Fneu_stim_chan2(:,trialsToRemove,:) = [];
        block.aligned_stim.F7_stim_chan2(:,trialsToRemove,:) = [];
        block.aligned_stim.df_f_chan2(:,trialsToRemove,:) = [];
        block.aligned_stim.zscore_chan2(:,trialsToRemove,:) = [];
    end

    %Markpoints
    if isfield(block, 'MarkpointsExperiment')
        block.MarkpointsExperiment.StimTable(trialsToRemove,:) = [];
    end
end

end %function