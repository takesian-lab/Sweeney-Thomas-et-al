function temp_block = setup_Cat_Var
% set up temp_block with empty variables to concatenate
% any variables left empty will be removed at the end of  simple_Concat_Blocks

    temp_block.Outcome = [];
    temp_block.trialType = [];
    temp_block.rxn_time = [];
    temp_block.water_delivery = [];
    temp_block.active_trials = [];
    temp_block.Sound_Time = [];
    temp_block.zcorr = [];

    %parameter variables
    temp_block.parameters.variable1 = [];
    temp_block.parameters.variable2 = [];
    temp_block.parameters.variable3 = [];
    temp_block.parameters.stimLength = [];
    temp_block.parameters.trialsToIgnore = [];
    
    %aligned_stim variables
    temp_block.aligned_stim.F_stim = [];
    temp_block.aligned_stim.F_stim_dff = [];
    temp_block.aligned_stim.F_stim_zscore = [];
    temp_block.aligned_stim.Fneu_stim = [];
    temp_block.aligned_stim.F7_stim = [];
    temp_block.aligned_stim.spks_stim = [];
    temp_block.aligned_stim.df_f = [];
    temp_block.aligned_stim.zscore = [];
    temp_block.aligned_stim.licks = [];
    temp_block.aligned_stim.velocity = [];
    temp_block.aligned_stim.whisker = [];
    temp_block.aligned_stim.pupil = [];
    
    %errorData variables
    temp_block.errorData.result = [];
    temp_block.errorData.result_orig = [];
    temp_block.errorData.error_trials = [];
    
    %Ops variables
    temp_block.ops.xoff = [];
    temp_block.ops.yoff = [];
 
    %Combined multiplane block variables
    temp_block.layer = [];
    temp_block.cell_number_orig = [];
    
    %Freq & AM detection behavior
    temp_block.holdingPeriod = [];
    temp_block.waitPeriod = [];
    temp_block.waitCondition = [];
end        