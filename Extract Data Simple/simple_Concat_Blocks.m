function [Blocks, BlockInfo, CellList, Archive] = simple_Concat_Blocks(Ops, Blocks, BlockInfo, CellList)
% concatenate blocks for simple_Extracted_Data
% Ops.catBlocks - 1 to perform concatenation, 0 to skip
% Ops.catRemoveBlocks - 1 to remove original (uncatenated) blocks from ExtractedData, 0 to keep in addition to concatenated blocks
% Ops.catBy - concatenation method:
%   - Session - concatenate all blocks from the same Tosca session
%   - Day - concatenate all blocks from the same day
%   - FOV - concatenate all blocks from FOVs with the same Suite2p analysis path
%   - TODO: Cell - concatenate cell data based on matched row number

% Version Control:
% v1 - archived
% v2 - current version: includes more functions and creates check for concatenating
% blocks that have different variables

%% Setup

Archive = {}; %Return empty archive if blocks don't need concatenating for any reason

if ~Ops.catBlocks % user must indicate to concatenate blocks
    disp('Skip concatenating blocks....')
    return
end

Hit_threshold = 0.5; %For freq_disc blocks (stim protocol 7)

removeBlock = zeros(height(BlockInfo),1);
removeCellList = zeros(height(CellList),1);
Days = unique(BlockInfo.Date);
Mice = unique(BlockInfo.MouseID);
fields = fieldnames(Blocks); %block names
channels = {'blue', 'green', 'red'};

%% Loop through each mouse to find blocks to concatenate

for i = 1:length(Mice)
    
    %For FOV concatenation, make sure we don't concatenate data from different expt_group
    %Combine analysis path and group into one string:
    orig_APs = BlockInfo.AnalysisPath;
    orig_Groups = BlockInfo.Group;
    combined_AP_with_Group = strcat(orig_APs, '\', orig_Groups);
    orig_APs_Cells = CellList.AnalysisPath;
    orig_Groups_Cells = CellList.Group;
    combined_AP_with_Group_Cells = strcat(orig_APs_Cells, '\', orig_Groups_Cells);
    
    M = find(BlockInfo.MouseID == Mice(i));
    MC = find(CellList.MouseID == Mice(i));
    APs = unique(combined_AP_with_Group(M));
    
    %Define variable to concatenate by
    switch Ops.catBy
        case "Day"
            var1 = Days;
        case "Session"
            var1 = unique(get_Tosca_variable_ExDat(Blocks));
        case "FOV"
            var1 = APs;
        otherwise
            error('variable for concatenation not present')
    end
    
    %Start concatenating
    for j = 1:length(var1)
        switch Ops.catBy
            case "Day"
                D = find(BlockInfo.Date == Days(j));
                DC = find(CellList.Date == Days(j));
            case "Session"
                y = strcat('Session',sprintf('%02d ',var1(j)));
                D = find(contains(BlockInfo.Block,y));
                DC = find(contains(CellList.Block,y));
            case "FOV"
                D = find(combined_AP_with_Group == APs(j));
                DC = find(combined_AP_with_Group_Cells == APs(j));
        end
        
        rows = intersect(D,M); % data to cat
        rowsC = intersect(DC,MC); % data to cat
        
        %If only one block was found, skip loop (except for stim protocol 7)
        if length(rows) <= 1 && ~(Ops.stim_protocol == 7)
            continue
        end
   
        %Found blocks to concatenate
        removeBlock(rows) = true;
        removeCellList(rowsC) = true;
        
        % stim specific elimination of rows--------------
        switch Ops.stim_protocol
            case 7 % frequency discrimination
                keeprows = true(size(rows));
                for k = 1:length(rows)
                    if Blocks.(fields{rows(k)}).prepTrial == 1 || Blocks.(fields{rows(k)}).HitRate<=Hit_threshold
                        keeprows(k) = false;
                    end
                end
                rows = rows(keeprows);
        end
        %  -----------------------------------------------
        
        %Skip loop if no rows left to concatenate
        if isempty(rows)
            continue
        end
        
        %% make new block!
        
        % update cell list:
        [tempCellList] = simple_cat_CellList(rows,BlockInfo,CellList);
        [temp_block] = setup_Cat_Var; % sets up all variables for concat
        F_cat = [];
        F_dff_cat = [];
        timestamp_cat = [];
        loco_cat = [];
        whisker_cat = [];
        pupil_cat = [];
        lick_cat = [];
        pupil_time_cat = [];
        check_stim_level = [];
        archive = {};

        
        %If concat type is FOV, do a quick check to make sure the cell numbers match
        cellNumberCheck = cell(1,length(rows));
        for nb = 1:length(rows) %number of blocks to concatenate
            cellNumberCheck{nb} = Blocks.(fields{rows(nb)}).cell_number;
        end
        for nb = 1:length(rows)
            for nb2 = 1:length(rows)
                if ~isequal(cellNumberCheck{nb}, cellNumberCheck{nb2})
                    error('Cell numbers do not match across FOV, check blocks and recompile')
                end
            end
        end

            
        %Loop through blocks to concatenate    
        for nb = 1:length(rows) %number of blocks to concatenate

            if nb == 1 %If first block, get new row for BlockInfo
                block_info_template = BlockInfo(rows(nb),:);
            end
             
            old_block = Blocks.(fields{rows(nb)});

            %check for the size of the old_block data-- should be the
            %same length as all the other ones
                
            archive.setup_archive{nb} = old_block.setup;
            temp_block.Outcome      = [temp_block.Outcome,old_block.Outcome];
            temp_block.trialType    = [temp_block.trialType,old_block.trialType];
            temp_block.rxn_time     = [temp_block.rxn_time,old_block.rxn_time];
            temp_block.TargetFreq   = old_block.TargetFreq;

            %from error data:
            if nb>1
                temp_block.errorData.error_trials = [temp_block.errorData.error_trials,(old_block.errorData.error_trials)+length(temp_block.errorData.result_orig)];
            else
                temp_block.errorData.error_trials = [temp_block.errorData.error_trials,old_block.errorData.error_trials];
            end
            
            if isfield(old_block.errorData, 'result')
                temp_block.errorData.result = [temp_block.errorData.result,old_block.errorData.result];
                temp_block.errorData.result_orig = [temp_block.errorData.result_orig,old_block.errorData.result_orig];
            end
                
            %optional variables:            
            if isfield(old_block, 'water_delivery')
                temp_block.water_delivery = [temp_block.water_delivery,old_block.water_delivery];
            end
                
            if isfield(old_block,'active_trials')
                temp_block.active_trials = [temp_block.active_trials,old_block.active_trials];
            end
                
            %Freq/AM detection behavior
            if any(Ops.stim_protocol == [13 19]) 
                temp_block.holdingPeriod = [temp_block.holdingPeriod, old_block.holdingPeriod];
                temp_block.waitPeriod    = [temp_block.waitPeriod, old_block.waitPeriod];
                temp_block.waitCondition = [temp_block.waitCondition, old_block.waitCondition];
                check_stim_level         = [check_stim_level,old_block.stim_level];
            end
                
            % parameters
            temp_block.parameters.variable1  = [temp_block.parameters.variable1,old_block.parameters.variable1];
            temp_block.parameters.variable2  = [temp_block.parameters.variable2,old_block.parameters.variable2];
            temp_block.parameters.variable3  = [temp_block.parameters.variable3,old_block.parameters.variable3];
            temp_block.parameters.stimLength = [temp_block.parameters.stimLength,old_block.parameters.stimLength];
            if isfield(old_block.parameters, 'trialsToIgnore')
                temp_block.parameters.trialsToIgnore = [temp_block.parameters.trialsToIgnore, old_block.parameters.trialsToIgnore];
                if Ops.stim_protocol ~= 12 && (size(temp_block.parameters.trialsToIgnore,2) ~= size(temp_block.parameters.variable1,2)); error('Accomodate data where not all blocks have trialsToIgnore'); end
            end
                
            % aligned to stim
            if isfield(old_block, 'FibPhot')
                temp_block.aligned_stim.F_stim        = cat(2,temp_block.aligned_stim.F_stim,old_block.aligned_stim.F_stim);
                temp_block.aligned_stim.F_stim_dff    = cat(2,temp_block.aligned_stim.F_stim_dff,old_block.aligned_stim.F_stim_dff);
                temp_block.aligned_stim.F_stim_zscore = cat(2,temp_block.aligned_stim.F_stim_zscore,old_block.aligned_stim.F_stim_zscore);
            end
            if isfield(old_block, 'cell_number') %S2P data
                temp_block.aligned_stim.F_stim        = cat(2,temp_block.aligned_stim.F_stim,old_block.aligned_stim.F_stim);
                temp_block.aligned_stim.Fneu_stim     = cat(2,temp_block.aligned_stim.Fneu_stim,old_block.aligned_stim.Fneu_stim);
                temp_block.aligned_stim.F7_stim       = cat(2,temp_block.aligned_stim.F7_stim,old_block.aligned_stim.F7_stim);
                temp_block.aligned_stim.spks_stim     = cat(2,temp_block.aligned_stim.spks_stim,old_block.aligned_stim.spks_stim);
                temp_block.aligned_stim.df_f          = cat(2,temp_block.aligned_stim.df_f,old_block.aligned_stim.df_f);
                temp_block.aligned_stim.zscore        = cat(2,temp_block.aligned_stim.zscore,old_block.aligned_stim.zscore);
            end
            temp_block.aligned_stim.licks = [temp_block.aligned_stim.licks;old_block.aligned_stim.licks];
            [temp_block.aligned_stim.velocity] = simple_cat_check_var(temp_block.aligned_stim,old_block.aligned_stim, 'velocity');
            [temp_block.aligned_stim.whisker] = simple_cat_check_var(temp_block.aligned_stim,old_block.aligned_stim, 'whisker');
            [temp_block.aligned_stim.pupil] = simple_cat_check_var(temp_block.aligned_stim,old_block.aligned_stim, 'pupil');

            %Markpoints experiments
            if isfield(old_block, 'MarkpointsExperiment')
                if ~isfield(temp_block, 'MarkpointsExperiment')
                    temp_block.MarkpointsExperiment.StimTable = old_block.MarkpointsExperiment.StimTable;
                else
                    temp_block.MarkpointsExperiment.StimTable = [temp_block.MarkpointsExperiment.StimTable; old_block.MarkpointsExperiment.StimTable];
                end
            end
                
            %Full traces
            if isfield(old_block, 'FibPhot')
                %Stack channels in matrix; if channel doesn't exist placeholder will be NaNs
                F_dff = nan(length(channels),length(old_block.FibPhot.timestamp_ds));
                F_org = nan(length(channels),length(old_block.FibPhot.timestamp_ds));
                for c = 1:length(channels)
                    if isfield(old_block.FibPhot, [channels{c} '_z'])
                        F_dff(c,:) = old_block.FibPhot.([channels{c} '_z']);
                        %elseif isfield(old_block.FibPhot, [channels{c} '_dff']) %to accomodate for older FibPhot blocks that didn't have the new z-scored traces
                        %F_dff(c,:) = old_block.FibPhot.([channels{c} '_dff']);
                    end
                    
                    if isfield(old_block.FibPhot, [channels{c} '_F'])
                        F_org(c,:) = old_block.FibPhot.([channels{c} '_F']);
                    end
                end
                [timestamp, F_new, full_loco] = match_fluor_loco_timestamps(F_org, old_block.loco_activity,old_block.FibPhot.timestamp_ds, old_block.locomotion_trace);
                [~, F_dff_new, ~] = match_fluor_loco_timestamps(F_dff, old_block.loco_activity,old_block.FibPhot.timestamp_ds, old_block.locomotion_trace);
                 [~, ~, full_licks] = match_fluor_loco_timestamps(F_dff_new, old_block.concat_licks,timestamp, old_block.concat_times);
              
                
            elseif isfield(old_block, 'F7') && isfield(old_block,'locomotion_trace') %2P
                [timestamp, F_new, full_loco] = match_fluor_loco_timestamps(old_block.F7, old_block.loco_activity, old_block.timestamp, old_block.locomotion_trace);
                [~, ~, full_licks] = match_fluor_loco_timestamps(F_new, old_block.concat_licks,timestamp, old_block.concat_times);
               
                try
                    [~, xoff, ~] = match_fluor_loco_timestamps(double(old_block.ops.xoff), old_block.loco_activity, old_block.timestamp, old_block.locomotion_trace);
                    [~, yoff, ~] = match_fluor_loco_timestamps(double(old_block.ops.yoff), old_block.loco_activity, old_block.timestamp, old_block.locomotion_trace);
                    temp_block.ops.xoff = [temp_block.ops.xoff, xoff];
                    temp_block.ops.yoff = [temp_block.ops.yoff, yoff];
                catch
                    %The error is that xoff and yoff lengths frame length used to be wrong
                    error('This part will throw an error if you have an older block. Recommended to recompile or rerun define_suite2p_singleblock')
                end
                
                if isfield(old_block, 'zcorr')
                    [~, zcorr, ~] = match_fluor_loco_timestamps(old_block.zcorr, old_block.loco_activity, old_block.timestamp, old_block.locomotion_trace);
                    temp_block.zcorr = [temp_block.zcorr, zcorr];
                end
                
            elseif isfield(old_block,'locomotion_trace') || isfield(old_block,'loco_times') %blocks without functional data
                timestamp = old_block.loco_times';
                full_loco = old_block.loco_activity';
                [~, ~, full_licks] = match_fluor_loco_timestamps(full_loco, old_block.concat_licks,timestamp, old_block.concat_times);
            else
                timestamp = old_block.timestamp;
                F_new = old_block.F7;
                [~, ~, full_licks] = match_fluor_loco_timestamps(F_new, old_block.concat_licks,timestamp, old_block.concat_times);
            end
            
            % licking data
             full_licks(find(isnan(full_licks))) = 0;
             full_licks(full_licks>0) = 1;
             lick_cat = [lick_cat,full_licks];
                
            % pupil/whisker
            if isfield(old_block, 'FibPhot') || isfield(old_block, 'F7')
                if isfield(old_block,'concat_whisker')
                    if isempty(whisker_cat) % if the block is missing pupil from all early experiments
                        
                        whisker_cat = nan(size(timestamp_cat));
                        pupil_time_cat = timestamp_cat;
                    end
                    [trimmed_pupiltime, ~, full_whisker] = match_fluor_loco_timestamps(F_new, old_block.concat_whisker, timestamp, old_block.pupil_timestamp, 1);
                    whisker_cat = [whisker_cat,full_whisker];
                    pupil_time_cat = [pupil_time_cat,trimmed_pupiltime];
                    
                elseif ~isempty(whisker_cat) && ~isfield(old_block, 'concat_whisker') % in case one block in the cat data are missing video
                    whisker_cat = [whisker_cat,nan(size(timestamp))];
                end
                
                if isfield(old_block,'concat_pupil')
                    if isempty(pupil_cat)
                        pupil_cat = nan(size(timestamp_cat));
                    end
                    [~, ~, full_pupil] = match_fluor_loco_timestamps(F_new, old_block.concat_pupil, timestamp, old_block.pupil_timestamp, 1);
                    pupil_cat = [pupil_cat,full_pupil];
                elseif ~isempty(pupil_cat) % in case one block in the cat data are missing video
                    pupil_cat = [pupil_cat,nan(size(timestamp))];
                end
                
            else % pupil/whisker data without functional data - upsample loco to match
                if isfield(old_block,'concat_whisker')
                    if isempty(whisker_cat) % if the block is missing pupil from all early experiments
                        
                        whisker_cat = nan(size(timestamp_cat));
                        pupil_time_cat = timestamp_cat;
                    end
                    [timestamp, full_whisker, full_loco] = match_fluor_loco_timestamps(old_block.concat_whisker, full_loco, old_block.pupil_timestamp, old_block.loco_times);
                    whisker_cat = [whisker_cat,full_whisker];
                    pupil_time_cat = [pupil_time_cat,timestamp];
                    if isfield(old_block,'concat_pupil')
                        if isempty(pupil_cat)
                            pupil_cat = nan(size(timestamp_cat));
                        end
                        [~, ~, full_pupil] = match_fluor_loco_timestamps(old_block.concat_whisker, old_block.concat_pupil, timestamp, old_block.pupil_timestamp,1);
                        pupil_cat = [pupil_cat,full_pupil];
                    elseif ~isempty(pupil_cat) % in case one block in the cat data are missing video
                        pupil_cat = [pupil_cat,nan(size(timestamp))];
                    end
                    
                elseif ~isempty(whisker_cat) && ~isfield(old_block, 'concat_whisker') % in case one block in the cat data are missing video
                    whisker_cat = [whisker_cat,nan(size(timestamp))];
                else
                    timestamp = old_block.loco_times;
                    full_loco = old_block.loco_activity;
                end
            end
                
            fr = old_block.setup.framerate;
            
            if nb == 1
                st = 0;
                tt = timestamp;
            else
                st =  timestamp_cat(end) + (1/fr); % last timestamp + 1 frame
                zerott = timestamp - timestamp(1); % zero out the current timestamps
                tt = zerott + st;
            end
            
            % dff trace for FibPhot
            if isfield(old_block, 'FibPhot')
                F_dff_cat = cat(2, F_dff_cat, F_dff_new);
            end
            
            % Raw trace (FibPhot) or F7 (S2p)
            if isfield(old_block, 'FibPhot') || isfield(old_block, 'F7')
                F_cat = cat(2,F_cat,F_new);
            end
            
            % these variables should be created even if there is no functional data
            tt = check_dim(tt);
            timestamp_cat = check_dim(timestamp_cat);
            timestamp_cat = [timestamp_cat, tt];
            if isfield(old_block,'locomotion_trace') || isfield(old_block,'loco_activity')
                loco_cat = check_dim(loco_cat);
                full_loco = check_dim(full_loco);
                loco_cat = [loco_cat,full_loco];
            end
            
            if isfield(old_block,'Sound_Time')
                S_T = (old_block.Sound_Time - timestamp(1)) + st;
                temp_block.Sound_Time = [temp_block.Sound_Time, S_T];
            else % blocks without functional data wont have a Sound_Time variable
                S_T = get_Sound_Time_behavior(old_block);
                S_T = S_T +st;
                S_T = check_dim(S_T);
                temp_block.Sound_Time = check_dim(temp_block.Sound_Time );
                temp_block.Sound_Time = [temp_block.Sound_Time, S_T];
            end
            
            % create a setup variable that has common values:
            % assumes all recordings done on same day
            % make new concat block name
            if nb == 1
                temp_block.setup = old_block.setup;
                if ~ismissing(temp_block.setup.FOV)
                    try
                        FOVtag = ['FOV' temp_block.setup.FOV{1} '_'];
                    catch
                        FOVtag = ['FOV' num2str(temp_block.setup.FOV(1)) '_'];
                    end
                else
                    FOVtag = '';
                end
                
                %Find diff between APs to make sure concat blocks don't have duplicate names
                if strcmp(Ops.catBy, "FOV")
                    shortAPnames = makeUniqueShortNames(var1);
                    APtag = strcat('_', shortAPnames{j});
                else
                    APtag = '';
                end
                
                new_block_name = convertStringsToChars(matlab.lang.makeValidName(strcat('CatBlock_', temp_block.setup.mousename, '_', FOVtag, temp_block.setup.expt_date, '_', temp_block.setup.stim_name, APtag)));
                temp_block.setup.block_filename = new_block_name;
                temp_block.setup.block_supname = (strcat(temp_block.setup.mousename, '-', temp_block.setup.expt_date, '-', temp_block.setup.stim_name));
                temp_block.setup.Tosca_path = 'cat_block';
                temp_block.setup.block_path = 'cat_block';
            end
        end % rows
            
        %% Finishing touches
        
        if isfield(old_block, 'FibPhot')
            temp_block.FibPhot.timestamp_ds = timestamp_cat;
            temp_block.timestamp = timestamp_cat;
            for c = 1:length(channels)
                if c == 1
                    temp_block.FibPhot.blue_F = F_cat(1,:);
                    temp_block.FibPhot.blue_dff = F_dff_cat(1,:);
                elseif c == 2
                    temp_block.FibPhot.green_F = F_cat(2,:);
                    temp_block.FibPhot.green_dff = F_dff_cat(2,:);
                elseif c ==3
                    temp_block.FibPhot.red_F = F_cat(3,:);
                    temp_block.FibPhot.red_dff = F_dff_cat(3,:);
                else
                    error('fib phot channel missing?')
                end
            end
        elseif isfield(old_block, 'F7') %2P
            temp_block.F7 = F_cat;
            temp_block.timestamp = timestamp_cat;
        else
            disp('no functional data to cat');
            temp_block.timestamp = timestamp_cat;
        end
        
        temp_block.loco_activity = loco_cat;
        temp_block.locomotion_trace = timestamp_cat;
        temp_block.pupil_timestamp =  pupil_time_cat;
        temp_block.concat_licks = lick_cat;
        
        % store pupil and whisker
        if ~isempty(pupil_cat)
            temp_block.concat_pupil = pupil_cat;
        end
        
        if ~isempty(whisker_cat)
            temp_block.concat_whisker = whisker_cat;
        end
        
        %Stim level
        if any(Ops.stim_protocol == [13 19]) && ~isempty(check_stim_level)
            if all(check_stim_level == check_stim_level(1))
                temp_block.stim_level = check_stim_level(1);
            else
                warning('Stim Levels are not the same in cat blocks')
                temp_block.stim_level = mode(check_stim_level);
            end
        end
            
        % remove empty variables from block
        temp_block_fields = fieldnames(temp_block);
        for f = 1:length(temp_block_fields)
            if isempty(temp_block.(temp_block_fields{f}))
                temp_block = rmfield(temp_block,temp_block_fields{f});
            end
        end
        
        % remove empty variables from block.aligned_stim
        temp_aligned_fields = fieldnames(temp_block.aligned_stim);
        for f = 1:length(temp_aligned_fields)
            if isempty(temp_block.aligned_stim.(temp_aligned_fields{f}))
                temp_block.aligned_stim = rmfield(temp_block.aligned_stim,temp_aligned_fields{f});
            end
        end
        
        % remove empty variables from parameters
        if isempty(temp_block.parameters.trialsToIgnore)
            temp_block.parameters = rmfield(temp_block.parameters,'trialsToIgnore');
        end

        block_info_template.Block = new_block_name;
        BlockInfo = [BlockInfo; block_info_template];
        Blocks.(new_block_name) = temp_block;
        Archive.(new_block_name) = archive;
        
        for c = 1:height(tempCellList)
            tempCellList.Block(c) = new_block_name;
        end
        CellList = [CellList;tempCellList];
        
        removeBlock(end+1) = false; % add this new "cat block" to the data set
        removeCellList(end+height(tempCellList)) = false; % add this new "cat block" to the data set
        
     end % days/session/etc
end % mice

%% Remove original blocks (optional)

BlockInfo.keepBlock = removeBlock;
CellList.removeRow = removeCellList;

if Ops.catRemoveBlocks
    removeBlock = logical(removeBlock );
    removeCellList = logical(removeCellList);
    BlockInfo(removeBlock,:) = [];
    blockfields = fieldnames(Blocks); %new fields containing new concatenated blocks
    rmbf = blockfields(removeBlock);
    for i =1:length(rmbf)
        Blocks = rmfield(Blocks,rmbf(i));
    end
    CellList(removeCellList,:) = [];
end

%Catch if concatenated blocks don't have unique names
if length(unique(BlockInfo.Block)) ~= length(BlockInfo.Block)
    error('Blocks do not have unique names')
end

end %function

function array2 = check_dim(array)

[dim1, dim2] = size(array);
if dim2 < dim1
    array2 = array';
else
    array2 = array;
end

end