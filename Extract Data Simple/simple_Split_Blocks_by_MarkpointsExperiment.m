function [BlockInfo, CellList, TrialData] = simple_Split_Blocks_by_MarkpointsExperiment(BlockInfo, CellList, TrialData, StimInfo, varargin)
%Split ExtractedData block variables by MarkpointsExperiment targeted ensemble

% PART 1:
% - Blocks with > 1 Markpoints ensemble (e.g. Sham and Sound Responsive) will be split based on ensemble trial identity
% - The split blocks will be given a new name starting with 'SplitBlock_'
% - As a result, these trials will be processed separately during extract data (e.g. sound-responsiveness will be assessed independently)

% PART 2 (OPTIONAL):
% - Remove trials where no targeted cells were activated

% PART 3 (OPTIONAL):
% - Split trials from SHAM blocks to estimate effects coming from repeated measurements alone
% - Use parse varargin options below to set rules for splitting
% - ***For now, this will DISCARD any non-sham ensembles***

%% Skip function if MarkpointsExperiment is not found in TrialData

if ~any(ismember(fields(TrialData), 'MarkpointsExperiment'))
    return;
end

%% Options

% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% Apply set of paramaters corresponding to the following analyses:
% 'Regular', 'Uncaging Shutter Control', 'Spatial XY', 'Spatial Z', 'Remove trials without activation'
addOptional(p, 'Analysis', 'Regular')

% Split ensembles based on name only (allow XY position to vary), used for  XY controls
addOptional(p, 'UseEnsembleNameOnly', 0);

% Remove trials where no target cells were activated
addOptional(p, 'RemoveTrials', 0);

%Settings for splitting Sham trials to do control analyses where you compare Sham to Sham

% Whether or not to split sham trials. If 1, all non-sham ensembles will be discarded
addOptional(p, 'SplitShamTrials', 0);

% How to perform split
addOptional(p, 'SplitType', 'Random'); %'Alternate' or 'Random'

% How many trials per stim condition required
addOptional(p, 'MinimumTrials', 20);

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% Apply Parameter Set [no need to supply additional varargin]
%USER: DO NOT SUPPLY YOUR OWN VARARGIN IF PROVIDING 'ANALYSIS'

switch ops.Analysis
    case {'Regular', 'Z Spatial'}
        %Do nothing, use default parameters above
        
    case 'Uncaging Shutter Control'
        ops.SplitShamTrials = 1;
        
    case 'XY Spatial'
        ops.UseEnsembleNameOnly = 1; %Use for XY spatial controls where the target position varies
        
    case 'Remove trials without activation'
        ops.RemoveTrials = 1; %Remove individual trials with no activated cells
        
    otherwise
        error('Unexpected analysis name provided')
end

%% PART 1: Loop through all blocks and split up trials by targeted ensemble

New_BlockInfo = table;
New_CellList = table;
New_TrialData = struct;

%Preallocate Ensemble columns
BlockInfo.Ensemble = strings(size(BlockInfo,1),1);
CellList.Ensemble = strings(size(CellList,1),1);

B = 1; %New block count

for b = 1:size(BlockInfo,1)
    cellListInd = strcmp(CellList.Block, BlockInfo.Block(b));
    
    %Figure out how many ensembles there are in MarkpointsExperiment
    MPE = TrialData(b).MarkpointsExperiment;
    newEnsembleList = MPE.EnsembleName; %Temporary naming including location info
    
    %Get info about uncaging shutter (1 = open, 0 = closed)
    UncShutter = unique(MPE.UncShutter);
    
    %DETECT ENSEMBLES WITH DIFFERENT PRAIRIE VIEW PARAMETERS
    %First, figure out if ensembles have the same given name but different PointNames or different number of points
    %This means they were separate ensembles in the Prairie View Markpoints window
    origEnsembles = unique(MPE.EnsembleName);
    if length(UncShutter) > 1
        %User concatenated blocks where uncaging shutter was open and closed
        error('Add code')
    elseif UncShutter == 1
        for E = 1:length(origEnsembles)
            ensembleRows = strcmp(MPE.EnsembleName, origEnsembles(E));
            temp_MPE = MPE(ensembleRows,:);
            temp_EnsembleList = temp_MPE.EnsembleName;

            uniquePointNames = unique(temp_MPE.PointNames);

            if length(uniquePointNames) > 1
                %Ensembles with different number of points found. Append this info to Ensemble name
                for m = 1:height(temp_MPE)
                    temp_EnsembleList(m) = strcat(temp_EnsembleList(m), '-', temp_MPE.PointNames(m));
                end

                %Replace data in MPE with updated ensembles
                newEnsembleList(ensembleRows) = temp_EnsembleList;
            end
        end
    elseif UncShutter == 0
        %If uncaging shutter is closed for the entire block, apply first row data to all points
        for m = 1:height(MPE)
            MPE(m,:) = MPE(1,:);
        end
    end
    
    %DETECT ENSEMBLES WITH THE SAME PRAIRIE VIEW PARAMETERS BUT TARGET POSITIONS MOVED
    %Set ops.UseEnsembleNameOnly to skip position analysis and just use ensemble names
    if ~ops.UseEnsembleNameOnly
        %For each ensemble, evaluate whether markpoint target locations differ
        %If so, do not treat them as the same ensemble unless:
        % - They have the same target cells (e.g. user just readjusted)
        % - They are a Sham ensemble and have moved by less than 20 microns
        % - They are a Sham ensemble with uncaging shutter closed
        
        origEnsembles = unique(newEnsembleList);

        for E = 1:length(origEnsembles)
            ensembleRows = strcmp(newEnsembleList, origEnsembles(E));
            temp_MPE = MPE(ensembleRows,:);
            temp_EnsembleList = newEnsembleList(ensembleRows);
            EnsembleUncShutter = unique(temp_MPE.UncShutter);

            X = cell2mat(temp_MPE.X);
            Y = cell2mat(temp_MPE.Y);
            XY = [X, Y];
            unique_XY = unique(XY,'rows');


            if length(EnsembleUncShutter) > 1
                %This means a Sham uncaging shutter closed ensemble was concatenated with a Sham uncaging shutter open ensemble
                %Map open shutter values onto closed shutter so that we can have more Sham trials in the same ensemble
                %Closed shutter will ALWAYS be Sham trials (automatically assigned during define_Markpoints_Experiment)
                UncShutterValue = temp_MPE.UncShutter;
                first1 = find(UncShutterValue == 1, 1, 'first');
                for m = 1:height(temp_MPE)
                    if UncShutterValue(m) == 0
                        temp_MPE.X{m} = temp_MPE.X{first1};
                        temp_MPE.Y{m} = temp_MPE.Y{first1};
                    end
                end
                EnsembleUncShutter = 1;

                %Assess again
                X = cell2mat(temp_MPE.X);
                Y = cell2mat(temp_MPE.Y);
                XY = [X, Y];
                unique_XY = unique(XY,'rows');
            end

            if size(unique_XY,1) > 1 && EnsembleUncShutter == 1
                %Target locations differ: determine whether to treat as the same or different 
                nLocations = size(unique_XY,1);
                shared_identity = nan(nLocations,1); %determine if ensembles are actually the same

                if ~strcmp(origEnsembles(E),'Sham')
                    %For non-sham trials, get Target cells for each location
                    correspondingTargets = cell(nLocations,1);
                    for i = 1:nLocations
                        for ii = 1:height(temp_MPE)
                            current_XY = XY(ii,:);
                            if isequal(unique_XY(i,:),current_XY)
                                break
                            end
                        end
                        correspondingTargets{i} = temp_MPE.TargetedCells{ii}';
                    end

                    %If target cells are the same, count as the same ensemble even if XY position differed
                    uniqueCells = unique(cell2mat(correspondingTargets),'rows');
                    for i = 1:length(shared_identity)
                        for ii = 1:height(uniqueCells)
                            if isequal(correspondingTargets{i}(:)', uniqueCells(ii,:))
                                shared_identity(i) = ii;
                            end
                        end
                    end
                else
                    %If this was a sham stimulation, allow the ensembles to differ by 20 microns
                    for i = 1:length(shared_identity)
                        for ii = 1:length(shared_identity)
                            if ii == i
                                continue;
                            end
                            if isnan(shared_identity(i)) && isnan(shared_identity(ii))
                                diff_XY = unique_XY(i,:) - unique_XY(ii,:);
                                if all(diff_XY < 20)
                                    shared_identity(i) = i;
                                    shared_identity(ii) = i;
                                end
                            elseif isnan(shared_identity(ii))
                                diff_XY = unique_XY(i,:) - unique_XY(ii,:);
                                if all(diff_XY < 20)
                                    shared_identity(ii) = shared_identity(i);
                                end
                            end
                        end
                    end
                end

                %For rows with no match, give unique identity
                count = max(shared_identity) +1;
                for i = 1:length(shared_identity)
                    if isnan(shared_identity(i))
                        shared_identity(i) = count;
                        count = count + 1;
                    end
                end

                %If different ensembles found, append this info to Ensemble name
                identity_list = nan(height(temp_MPE),1);
                for m = 1:height(temp_MPE)
                    for n = 1:nLocations
                        if isequal(XY(m,:),unique_XY(n,:))
                            identity_list(m) = shared_identity(n);
                        end
                    end
                end
                for m = 1:height(temp_MPE)
                    temp_EnsembleList(m) = strcat(temp_EnsembleList(m), '-', num2str(identity_list(m)));
                end

            elseif size(unique_XY,1) > 1 && EnsembleUncShutter == 0
                %FOR SHAM TRIALS ONLY
                %Target locations differ but uncaging shutter was closed
                %Use first row of data for X/Y positions and do not change Ensemble name
                for m = 1:height(temp_MPE)
                    temp_MPE.X{m} = temp_MPE.X{1};
                    temp_MPE.Y{m} = temp_MPE.Y{1};
                end
            end

            %Replace data in MPE with updated ensembles
            newEnsembleList(ensembleRows) = temp_EnsembleList;
            MPE(ensembleRows,:) = temp_MPE;
        end

        %Save MPE with changes
        TrialData(b).MarkpointsExperiment = MPE;
    end
    
    %Get unique ensembles
    Ensembles = unique(newEnsembleList);

    %SPLIT ENSEMBLES
    if length(Ensembles) == 1 %Only one ensemble found
        
        %If there is only one ensemble, we just need to save the ensemble name and copy over the block data
        BlockInfo.Ensemble(b) = Ensembles;
        CellList.Ensemble(cellListInd) = Ensembles;
        
        New_BlockInfo(B,:) = BlockInfo(b,:);
        New_CellList = [New_CellList; CellList(cellListInd,:)];
        if B == 1
            New_TrialData = TrialData(b);
        else
            New_TrialData(B) = TrialData(b);
        end
        
        B = B + 1;
        
    else
        for e = 1:length(Ensembles) %More than one ensemble found
            
            %If there is more than one ensemble, loop through ensembles and save each as a separate block
            ensembleRows = strcmp(newEnsembleList, Ensembles{e});

            %Make new unique block name for split blocks
            orig_block_name = BlockInfo.Block(b);
            new_block_name = char(strcat('SplitBlock', num2str(e), {'_'}, BlockInfo.Block(b)));
                        
            %Get original ensemble name
            origEnsemble = unique(MPE.EnsembleName(ensembleRows)); %We never changed it in MPE
            
            %Temporarily assign new block name and ensemble (easier than finding what position we are at in CellList)
            BlockInfo.Ensemble(b) = origEnsemble;
            CellList.Ensemble(cellListInd) = origEnsemble;
            BlockInfo.Block(b) = new_block_name;
            CellList.Block(cellListInd) = new_block_name;
            TrialData(b).Block = new_block_name;
            
            %Save new BlockInfo and CellList
            New_BlockInfo(B,:) = BlockInfo(b,:);
            New_CellList = [New_CellList; CellList(cellListInd,:)];
            
            %Copy over TrialData as is...
            if B == 1
                New_TrialData = TrialData(b);
            else
                New_TrialData(B) = TrialData(b);
            end
            
            %Replace original block name (next loop will use this to generate new_block_name)
            BlockInfo.Block(b) = orig_block_name;
            
            %Then loop through variables and only keep trials corresponding to current Ensemble
            variablesToUpdate = {'Stim_V1', 'Stim_V2', 'Stim_V3', 'StimID', 'Stim_Length', 'IsBlank', 'Sound_Time',...
                'Loco', 'Licks', 'Pupil', 'Whisker', 'IsRunning', 'Trials', 'Spikes', 'MarkpointsExperiment'...
                'Outcomes', 'ReactionTime', 'HoldingPeriod', 'WaitPeriod', 'WaitCondition', 'Result', 'Result_orig'};
            
            %Remove variables not in TrialData (e.g. Stim_V3 might not be included)
            variablesToUpdate(~ismember(variablesToUpdate,fields(TrialData))) = [];
            
            for v = 1:length(variablesToUpdate)
                New_TrialData(B).(variablesToUpdate{v}) = keepRows(New_TrialData(B).(variablesToUpdate{v}), ensembleRows);
            end
            
            B = B + 1;
        end
    end
end

BlockInfo = New_BlockInfo;
CellList = New_CellList;
TrialData = New_TrialData;

%% PART 2: Remove trials where no targets were activated (Optional)
%FYI: This code is mostly identical to the code in evaluate_markpoints_ExtractedData

if ~ops.SplitShamTrials && ops.RemoveTrials

    %% Parameters
    Z_level = 1; %For single trials
    AUC_level = 1;
    smTime = 0;
    LatePeak = 0.5; %For target cell activation, be more strict on the LatePeak criteria for activated vs. prolonged responses 
    
    New_BlockInfo = table;
    New_CellList = table;
    New_TrialData = TrialData(1); %Start variable using first row that will be replaced
    B = 1; %New block count

    for b = 1:height(BlockInfo)

        blockname = BlockInfo.Block(b);
        blockind = strcmp(CellList.Block, blockname);
        block_CellList = CellList(blockind,:);
        
        %If this is a sham block, skip
        if strcmp(BlockInfo.Ensemble(b),"Sham")
            New_BlockInfo(B,:) = BlockInfo(b,:);
            New_CellList = [New_CellList; block_CellList];
            New_TrialData(B) = TrialData(b);
            B = B+1;
            continue;
        end

        %Get target indices
        MPE = TrialData(b).MarkpointsExperiment;
        targetedCells = MPE.TargetedCells{1};
        targetedCellIndices = nan(size(targetedCells));
        for t = 1:length(targetedCells)
            targetedCellIndices(t) = block_CellList.CellOrder(block_CellList.CellNumber == targetedCells(t));
        end

        Trials = TrialData(b).Trials(targetedCellIndices,:,:);
        StimID = TrialData(b).StimID;

        %Find corresponding sham block

        %First look for matching block index (e.g. interleaved blocks)
        origBlockName = BlockInfo.ShortBlockName(b); %Original block name before split by markpoints
        shamBlockInd = strcmp(BlockInfo.ShortBlockName, origBlockName);
        shamBlockInd(b) = 0; %Remove current block
        shamBlockInd(~strcmp(BlockInfo.Ensemble,"Sham")) = 0; %Limit to sham blocks

        %If there are no matching blocks, use analysis path
        if ~any(shamBlockInd)
            analysisPath = BlockInfo.AnalysisPath(b);
            shamBlockInd = strcmp(BlockInfo.AnalysisPath, analysisPath);
            shamBlockInd(b) = 0; %Remove current block
            shamBlockInd(~strcmp(BlockInfo.Ensemble,"Sham")) = 0; %Limit to sham blocks
        end

        %Prepare to store data
        unique_stim = unique(StimID);
        nStim = length(unique_stim);
        nCells = length(targetedCellIndices);
        averageShamResponse = nan(nCells,nStim,StimInfo.nFrames);

        if sum(shamBlockInd) > 1
            %Confirm there's only one matching sham block for this AP
            %This should be the case if concatenate by FOV was selected
            error('Expecting only one sham block per ACT block')

        elseif sum(shamBlockInd) == 1
            %Subtract sham trials from act trials to find 'true' activation
            sham_Trials = TrialData(shamBlockInd).Trials(targetedCellIndices,:,:);
            sham_StimID = TrialData(shamBlockInd).StimID;

            %Go through each stimulus to compute SHAM average and subtract from individual ACT trials
            %This is important to do by stimulus because the cell could be tuned differently to each stim
            delta_Trials = nan(size(Trials));
            for c = 1:nCells
                for s = 1:nStim
                    averageShamResponse(c,s,:) = squeeze(mean(sham_Trials(c,sham_StimID == unique_stim(s),:),2,'omitnan'));
                    delta_Trials(c,StimID == unique_stim(s),:) = Trials(c,StimID == unique_stim(s),:) - averageShamResponse(c,s,:);
                end
            end

        else
            %If no sham trials available, use raw trials
            delta_Trials = Trials;
        end

        %% Compute peak responsiveness on delta_Trials for each trial

        %Store in MPE
        for i = 1:height(MPE)
            current_trials = squeeze(delta_Trials(:,i,:));
            [PeakData, ~] = simple_check_if_responsive(current_trials, StimInfo.nBaselineFrames, StimInfo.fs, Z_level, AUC_level, smTime, 0, 'LatePeak', LatePeak);

            MPE.Peak{i} = PeakData.Peak';
            MPE.Peak_AUC{i} = PeakData.Peak_AUC';
            MPE.IsActivated{i} = strcmp(PeakData.ResponseType','activated');
            MPE.nActivated(i) = sum(MPE.IsActivated{i});
        end
        TrialData(b).MarkpointsExperiment = MPE;
            
        %Then loop through variables and only keep trials with activated targets
        keepTrial = sum(cell2mat(MPE.IsActivated),2) > 0;

        if any(keepTrial) %If no trials found to keep, block will not be included
            variablesToUpdate = {'Stim_V1', 'Stim_V2', 'Stim_V3', 'StimID', 'Stim_Length', 'IsBlank', 'Sound_Time',...
                'Loco', 'Licks', 'Pupil', 'Whisker', 'IsRunning', 'Trials', 'Spikes', 'MarkpointsExperiment'...
                'Outcomes', 'ReactionTime', 'HoldingPeriod', 'WaitPeriod', 'WaitCondition', 'Result', 'Result_orig'};

            %Remove variables not in TrialData (e.g. Stim_V3 might not be included)
            variablesToUpdate(~ismember(variablesToUpdate,fields(TrialData))) = [];

            New_TrialData(B) = TrialData(B);
            for v = 1:length(variablesToUpdate)
                New_TrialData(B).(variablesToUpdate{v}) = keepRows(New_TrialData(b).(variablesToUpdate{v}), keepTrial);
            end
            New_BlockInfo(B,:) = BlockInfo(b,:);
            New_CellList = [New_CellList; block_CellList];
            
            B = B+1;
        else
            warning(strcat({'No activated targets found in '}, blockname, '. Removing from experiment.'))
        end
    end

    BlockInfo = New_BlockInfo;
    CellList = New_CellList;
    TrialData = New_TrialData;
 
end

%% PART 3: Split sham trials (Optional)

if ~ops.SplitShamTrials
    return;
end

%Remove all blocks that are not Sham
isShamBlock = strcmp(BlockInfo.Ensemble, 'Sham');
BlockInfo = BlockInfo(isShamBlock,:);
TrialData = TrialData(isShamBlock);
CellList = CellList(strcmp(CellList.Ensemble, 'Sham'),:);

New_BlockInfo = table;
New_CellList = table;
New_TrialData = struct;

B = 1; %New block count

for b = 1:size(BlockInfo,1)
    
    cellListInd = strcmp(CellList.Block, BlockInfo.Block(b));
    StimID = TrialData(b).StimID;
    
    %Figure out how many sets of Sham trials we can make from this block
    tabulate_stim = tabulate(StimID);
    stimCount = tabulate_stim(:,2); %Number of trials per stim
    
    %To split Sham trials, we need at least 2*ops.MinimumTrials fo each stim
    if ~all(stimCount >= 2*ops.MinimumTrials)
        continue; %Move onto next block if we don't have this
    end
        
    %Find max number of sets of Sham trials we can make
    nSets = 2; nSetsPossible = 1;
    while nSetsPossible == 1
        if ~all(stimCount/nSets >= ops.MinimumTrials)
            nSetsPossible = 0;
        else
            nSets = nSets + 1;
        end
    end
    nSets = nSets - 1;
    
    %Determine which rows to assign to each set based on ops.SplitType
    uniqueStim = unique(StimID);
    nStim = length(uniqueStim);
    SetRows = cell(1,nSets); %Preallocate
        
    switch ops.SplitType

        case 'Alternate'
            %Take every Nth trial (per stim)
            %The benefit of this option is the result will be the same every time
            %Stimuli are also already randomized within Tosca, so Alternate is already pseudo-random
            for n = 1:nSets
                rows = [];
                for s = 1:nStim
                    stimRows = find(StimID == uniqueStim(s));
                    ind = (1:nSets:length(stimRows)) + (n-1);
                    ind(ind > length(stimRows)) = [];
                    rows = [rows; stimRows(ind)];
                end
                SetRows{n} = sort(rows,'ascend');
            end
            
        case 'Random'
            %Randomize set allocation within stimulus type
            setInd = nan(size(StimID));
            for s = 1:nStim
                %Get all rows corresponding to each stimulus
                stimRows = find(StimID == uniqueStim(s));
                %Randomly permute row order
                randInd = randperm(length(stimRows));
                randStimRows = stimRows(randInd);
                
                temp_setInd = repmat(1:nSets, [1,floor(length(randStimRows)/nSets)]);
                for n = 1:length(temp_setInd)
                    setInd(randStimRows(n)) = temp_setInd(n);
                end
            end
            for n = 1:nSets
                SetRows{n} = sort(find(setInd == n),'ascend');
            end

        otherwise
            error('SplitType not recognized')
    end

    %Check that trials are unique and meet criteria
    for n = 1:nSets
        checkRows = SetRows{n};
        checkStim = StimID(checkRows);
        tabulate_stim = tabulate(checkStim);
        if ~all(tabulate_stim(:,2) >= ops.MinimumTrials)
            error('Expected number of trials not reached')
        end
        otherSets = 1:nSets;
        otherSets(otherSets == n) = [];
        for n2 = 1:length(otherSets)
            intersect_sets = intersect(checkRows, SetRows{otherSets(n2)});
            if ~isempty(intersect_sets)
                error('Duplicate trials found between sets')
            end
        end
    end
                    
    %Make a new block for each set
    for n = 1:nSets
        if n == 1%nSets
            ensembleName = 'Sham'; %Keep one set as Sham
        else
            ensembleName = strcat('Control', num2str(n-1));
        end
        
        %Make new unique block name for split blocks
        orig_block_name = BlockInfo.Block(b);
        new_block_name = char(strcat('ControlBlock', num2str(n), {'_'}, BlockInfo.Block(b)));
     
        %Determine what rows to assign to each set
        ensembleRows = SetRows{n};

        %Temporarily assign new block name and ensemble (easier than finding what position we are at in CellList)
        BlockInfo.Ensemble(b) = ensembleName;
        CellList.Ensemble(cellListInd) = ensembleName;
        BlockInfo.Block(b) = new_block_name;
        CellList.Block(cellListInd) = new_block_name;
        TrialData(b).Block = new_block_name;

        %Save new BlockInfo and CellList
        New_BlockInfo(B,:) = BlockInfo(b,:);
        New_CellList = [New_CellList; CellList(cellListInd,:)];
            
            %Copy over TrialData as is...
            if B == 1
                New_TrialData = TrialData(b);
            else
                New_TrialData(B) = TrialData(b);
            end
            
            %Replace original block name (next loop will use this to generate new_block_name)
            BlockInfo.Block(b) = orig_block_name;
            
            %Then loop through variables and only keep trials corresponding to current Ensemble
            variablesToUpdate = {'Stim_V1', 'Stim_V2', 'Stim_V3', 'StimID', 'Stim_Length', 'IsBlank', 'Sound_Time',...
                'Loco', 'Licks', 'Pupil', 'Whisker', 'IsRunning', 'Trials', 'Spikes', 'MarkpointsExperiment'...
                'Outcomes', 'ReactionTime', 'HoldingPeriod', 'WaitPeriod', 'WaitCondition', 'Result', 'Result_orig'};
            
            %Remove variables not in TrialData (e.g. Stim_V3 might not be included)
            variablesToUpdate(~ismember(variablesToUpdate,fields(TrialData))) = [];
            
            for v = 1:length(variablesToUpdate)
                New_TrialData(B).(variablesToUpdate{v}) = keepRows(New_TrialData(B).(variablesToUpdate{v}), ensembleRows);
            end
            
            %Update ensemble name in TrialData
            New_TrialData(B).MarkpointsExperiment.EnsembleName(:) = ensembleName;
            B = B + 1;
    end
end

BlockInfo = New_BlockInfo;
CellList = New_CellList;
TrialData = New_TrialData;
        
end %function

%% check dimensions of TrialData variable and subset by ensembleRows

function subsetVariable = keepRows(variable, ensembleRows)

if numel(size(variable)) == 2 %1 or 2D
    subsetVariable = variable(ensembleRows,:);  
elseif numel(size(variable)) == 3 %3D
    subsetVariable = variable(:,ensembleRows,:);
end

end