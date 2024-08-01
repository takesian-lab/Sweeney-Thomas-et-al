%concatenate_blocks
%Run this script to concatenate multiple blocks (e.g. from the same imaging/behavior day)
%For 2P data, the FOVs must have been run together in suite2p such that the cell numbers are identical
%By MET, CGS May 2022

% WORK IN PROGRESS: THIS IS AN AMALGAMATION OF A FEW DIFFERENT SCRIPTS WITH
% THE GOAL OF CONCATENATING BLOCKS. NOT FINISHED WHATSOEVER

%% Define Info sheet and compiled blocks path

PC_name = getenv('computername');

switch PC_name
    case 'RD0366' %Maryse
        info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\NxDxNpC061921F4';
        blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\NxDxNpC061921F4\CompiledBlocks';
        info_filename = 'ConcatTest';
        stimProtocol = 13;

    case 'RD0332' %Carolyn
        info_path = 'Z:\Carolyn\Fiber Photometry\Behavior\5HT sensor\Info sheets';
        blocks_path = 'Z:\Carolyn\Fiber Photometry\Behavior\5HT sensor\Compiled Blocks';
        info_filename = 'Info_VxDK123121M2';
        stimProtocol = 7;
        
    case 'RD-6-TAK2' %Anna's computer
        info_path = 'Z:\Carolyn\Fiber Photometry\Behavior\5HT sensor\Info sheets';
        blocks_path = 'Z:\Carolyn\Fiber Photometry\Behavior\5HT sensor\Compiled Blocks';
        info_filename = 'Info_VxDK123121M2';
        stimProtocol = 7;
        
    case 'TakesianArea51-1' %Christine's computer
        info_path = 'Z:\Carolyn\Fiber Photometry\Behavior\5HT sensor\Info sheets';
        blocks_path = 'Z:\Carolyn\Fiber Photometry\Behavior\5HT sensor\Compiled Blocks';
        info_filename = 'Info_VxDK123121M2';
        stimProtocol = 7;


    otherwise
        disp('Computer does not match known users')
        return
end

%% Load blocks corresponding to stimProtocol

cd(info_path)
Info = importfile(info_filename);
[blockData] = fillSetupFromInfoTable_v3(Info, blocks_path, stimProtocol,0);

% find size of data and preallocate
setup = blockData.setup;
mouse = setup.mousename{1};

D = blockData.([mouse]);
P = D.parameters;
D = rmfield(D,'parameters');
fields = fieldnames(D);

%% Create new block struct using the first block as a template
combined_block = struct;
combined_block.files = true;
combined_block.files = {block1.setup.block_filename, block2.setup.block_filename};
combined_block.setup = block1.setup;
combined_block.setup2 = block2.setup;
filename = strcat('Combined_', block1.setup.block_filename, '.mat');

%% Combine Tosca data
combined_block = hcat(combined_block, block1, block2, 'New_sound_times');
combined_block = hcat(combined_block, block1, block2, 'start_time');
combined_block = vcat(combined_block, block1, block2, 'lick_time');
combined_block = hcat(combined_block, block1, block2, 'Tosca_times');
combined_block = hcat(combined_block, block1, block2, 'Outcome');
combined_block = hcat(combined_block, block1, block2, 'trialType');
combined_block = check_equal(combined_block, block1, block2, 'TargetFreq');
parameters = struct;
combined_block.parameters = hcat(parameters, block1.parameters, block2.parameters, 'variable1');
combined_block.parameters = hcat(combined_block.parameters, block1.parameters, block2.parameters, 'variable2');
combined_block = vcat(combined_block, block1, block2, 'loco_data');
combined_block = vcat(combined_block, block1, block2, 'loco_activity');
combined_block = vcat(combined_block, block1, block2, 'loco_times');
combined_block = vcat(combined_block, block1, block2, 'locomotion_data');
combined_block = vcat(combined_block, block1, block2, 'active_time');
combined_block = hcat(combined_block, block1, block2, 'Sound_Time');
combined_block = vcat(combined_block, block1, block2, 'locTime');

   
combinedBlock = struct;
combinedBlock.setup = block1.setup;
if ~isfield(block1, 'setup2')
    combinedBlock.setup2 = block2.setup;
elseif ~isfield(block1, 'setup3')
    combinedBlock.setup2 = block1.setup2;
    combinedBlock.setup3 = block2.setup;
else
    error('More than 2 blocks found, combine manually')
end
combinedBlock.Combined = true;
combinedBlock.start_time = [block1.start_time, block2.start_time];
combinedBlock.Tosca_times = [block1.Tosca_times, block2.Tosca_times];
combinedBlock.errors = [block1.errors, block2.errors];
combinedBlock.New_sound_times = [block1.New_sound_times, block2.New_sound_times];
combinedBlock.New_sound_idx = [block1.New_sound_idx, block2.New_sound_idx];
combinedBlock.lick_time = [block1.lick_time; block2.lick_time];
combinedBlock.concat_times = [block1.concat_times; block2.concat_times + block1.concat_times(end)];
combinedBlock.concat_licks = [block1.concat_licks; block2.concat_licks];
combinedBlock.Outcome = [block1.Outcome, block2.Outcome];
combinedBlock.trialType = [block1.trialType, block2.trialType];
combinedBlock.TargetFreq = [block1.TargetFreq, block2.TargetFreq];
combinedBlock.parameters.variable1 = [block1.parameters.variable1, block2.parameters.variable1];
combinedBlock.parameters.variable2 = [block1.parameters.variable2, block2.parameters.variable2];
%skip combinedBlock.loco_data
combinedBlock.loco_activity = [block1.loco_activity; block2.loco_activity];
combinedBlock.loco_times = [block1.loco_times; block2.loco_times + block1.loco_times(end)];
combinedBlock.loc_Trial_times = [block1.loc_Trial_times, block2.loc_Trial_times];
combinedBlock.loc_Trial_activity = [block1.loc_Trial_activity, block2.loc_Trial_activity];
combinedBlock.rxn_time = [block1.rxn_time, block2.rxn_time];
combinedBlock.locIDX = [block1.locIDX, block2.locIDX];
combinedBlock.holdingPeriod = [block1.holdingPeriod, block2.holdingPeriod];
combinedBlock.stim_level = unique([block1.stim_level, block2.stim_level]);

%% Combine Bruker data
combined_block = vcat(combined_block, block1, block2, 'timestamp');
combined_block = hcat(combined_block, block1, block2, 'isLoco');

%% Combine Suite2p data
%combined_block.ops = check_equal(combined_block, block1, block2, 'ops'); %Does not need to be equal
combined_block = check_equal_struct(combined_block, block1, block2, 'img');
combined_block = check_equal_struct(combined_block, block1, block2, 'iscell');
combined_block = check_equal_struct(combined_block, block1, block2, 'cell_number');
combined_block = check_equal_struct(combined_block, block1, block2, 'stat');
combined_block = hcat(combined_block, block1, block2, 'F');
combined_block = hcat(combined_block, block1, block2, 'Fneu');
combined_block = hcat(combined_block, block1, block2, 'spks');
combined_block = check_equal_struct(combined_block, block1, block2, 'redcell');

%% Combine aligned to stim data
try
    combined_block = vcat(combined_block, block1, block2, 'aligned_to_stim');
catch
    warning('aligned_to_stim mats do not match on all dimensions. Not included in combined data')
end

      
        mouseID=setup.mousename{a,b};
        Imaging_Block=setup.Imaging_sets{a,b};
 
        for i=1:length(Imaging_Block)
            
            unique_block_name = setup.unique_block_names{a,b}(i);
            block = data.([mouseID]).([unique_block_name]);
            
            if ismissing(block.setup.suite2p_path)
                disp('Skipping Suite2p data for block...');
                disp(unique_block_name);
                return
            end
            
            if setup.run_redcell ==1
              data.([mouseID]).parameters.red_idx  = block.redcell;
              data.([mouseID]).parameters.green_idx = ~block.redcell;
            end
            
            
%             timestamp =  block.timestamp;
%             Sound_Time = block.Sound_Time;
            isLoco = block.isLoco;
            
            
            % correct for the third dimension of data.
            if i == 1
                F_cat = block.aligned_stim.F_stim;
                F7_cat = block.aligned_stim.F7_stim;
                spks_cat = block.aligned_stim.spks_stim;
                neu_cat = block.aligned_stim.Fneu_stim;
                V1_cat = block.parameters.variable1;
                V2_cat = block.parameters.variable2;
                isLoco_cat = block.isLoco;
            else i>1
                % check to make sure that the size (frames) of F, F7,
                % spks,and neu is the same across blocks. The third dim of
                % each of these matricies should be the same, so I only
                % used F_cat to test for size issues:
                if size(F_cat,3) == size(block.aligned_stim.F_stim,3)
                    F_cat = cat(2,F_cat,block.aligned_stim.F_stim);
                    F7_cat = cat(2,F7_cat,block.aligned_stim.F7_stim);
                    spks_cat = cat(2,spks_cat,block.aligned_stim.spks_stim);
                    neu_cat = cat(2,neu_cat,block.aligned_stim.Fneu_stim);
                    
                elseif size(F_cat,3) > size(block.aligned_stim.F_stim)
                    % make new array dim3 same size as F_cat dim3
                    diff_size = size(F_cat,3) - size(block.aligned_stim.F_stim,3);
                    add_size = NaN(size(F_cat,1),size(F_cat,2),size(diff_size));
                    
                    % pad the arrays with Nans
                    F_nan = cat(3,block.aligned_stim.F_stim,add_size);
                    F7_nan =  cat(3,block.aligned_stim.F7_stim,add_size);
                    spks_nan = cat(3,block.aligned_stim.spks_stim,add_size);
                    neu_nan =  cat(3,block.aligned_stim.Fneu_stim,add_size);
                    
                    % concat the arrays
                    F_cat = cat(2, F_cat,F_nan);
                    F7_cat = cat(2,F7_cat,F7_nan);
                    spks_cat = cat(2,spks_cat,spks_nan);
                    neu_cat = cat(2,neu_cat,neu_nan);
                    
                else size(F_cat,3) < size(block.aligned_stim.F_stim)
                    %make dim3 larger for F_cat...
                    diff_size = size(block.aligned_stim.F_stim,3)- size(F_cat,3); 
                    add_size = NaN(size(F_cat,1),size(F_cat,2),size(diff_size));
                    F_cat = cat(3,F_cat,add_size);
                    F7_cat =  cat(3,F7_cat,add_size);
                    spks_cat = cat(3,spks_cat,add_size);
                    neu_cat =  cat(3,neu_cat,add_size);
                    
                    % then concat the arrays...
                    F_cat = cat(2, F_cat,block.aligned_stim.F_stim);
                    F7_cat = cat(2,F7_cat,block.aligned_stim.F7_stim);
                    spks_cat = cat(2,spks_cat,block.aligned_stim.spks_stim);
                    neu_cat = cat(2,neu_cat,block.aligned_stim.Fneu_stim);
                end
                % cat is Loco and variables...
                V1_cat = cat(2,V1_cat,block.parameters.variable1);
                V2_cat = cat(2,V2_cat,block.parameters.variable2);
                isLoco_cat = cat(2,isLoco_cat,block.isLoco);
            end     
        end
    end

    data.([mouseID]).cat.F_cat = F_cat;
    data.([mouseID]).cat.F7_cat = F7_cat;
    data.([mouseID]).cat.spks_cat = spks_cat;
    data.([mouseID]).cat.neu_cat = neu_cat;
    data.([mouseID]).cat.V1_cat = V1_cat;
    data.([mouseID]).cat.V2_cat = V2_cat;
    data.([mouseID]).cat.isLoco_cat = isLoco_cat;
    
    

%% put together blocks
% elements to be concatenated are:
%   - stim types (Var1 and Var2)
%   - stim traces
%   - more to come...



stimTrace1 = [];
stimTrace2 = [];
varstim = [];
dB = [];
if stimProtocol ==7
    outcome = [];
end
LoCat = []; Licks = [];
for i = 1:length(fields)
    
    DM = D.([fields{i}]);
    
    if i ==1
        baseline_length = DM.setup.constant.baseline_length;
        after_stim = DM.setup.constant.after_stim;
        fs = DM.setup.framerate;
        basepoints = baseline_length*fs;
        timelength  = after_stim*fs;
        graphx = linspace((baseline_length.*-1),after_stim,size(DM.aligned_stim.F_stim,3));
    end
    stimTrace1 = cat(2,stimTrace1, DM.aligned_stim.F_stim);
    stimTrace2 = cat(2,stimTrace2, DM.aligned_stim.F_stim_dff);
    Licks      = cat(1,Licks,DM.aligned_stim.licks);
    varstim = [varstim, DM.parameters.variable1];
    dB = [dB, DM.parameters.variable2];
    LoCat = [LoCat, DM.active_trials];
    if stimProtocol ==7
        outcome = [outcome, DM.Outcome];
    end
end
