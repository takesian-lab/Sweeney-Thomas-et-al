
function  [D , params, ztraces,ResidTraces] = gSetup(bl,Data,s,blist,Ops,bsfr,t)
% gSetup will generate regressors and stimuli indicies for the GLM
% generated from simple Extracted Data.
% Version 2; updated Carolyn Sweeney 10/24/2023


% Input:
% bl        :index of block within Data.StimType.BlockInfo
% Data      : ExtractedData files used in modeling anlysis
% s         : stim type (ie 'FM')
% blist     : list of blocks
% Ops       : Ops settings 
% bsfr      : frames of baseline data
% t         : timestamp for block


%
%
% Output:
% D         : Stim variables used to generate regressors
% params    : parameters for regression
% ztraces   : zscored/normalized calcium traces for the block, if
%             Ops.ByTrials this variable is concatenated traces 
% ResidTraces: traces for residuals

%-------------------------------------------------------------------------
% determine if data are from FibPhot block or not:

 % sort out blanks trials and find stim V1/V2
D.blankIDX = find(Data.([s]).TrialData(bl).IsBlank==1);
D.stimIDX = find(Data.([s]).TrialData(bl).IsBlank==0);
D.uniqueV1 = Data.([s]).StimInfo.V1;
D.uniqueV2 = Data.([s]).StimInfo.V2;
D.Sound_Time =  Data.([s]).TrialData(bl).Sound_Time(~Data.([s]).TrialData(bl).IsBlank);
D.locoTrace = Data.([s]).TrialData(bl).Full_Loco'; % motor activity
D.rows = find(Data.([s]).CellList.Block == blist{bl}); % rows in Data.S that correspond to current block
%--------------------------------------------------------------------------
% store V1 
if Ops.useV1 
    for v = 1:length(D.uniqueV1)
        if Ops.ByTrials
            v1StimTemp = Data.([s]).TrialData(bl).Stim_V1(D.stimIDX);
            D.V1_IDX{v} = find(v1StimTemp == D.uniqueV1(v));
        else
            D.V1_IDX{v} = find(Data.([s]).TrialData(bl).Stim_V1 == D.uniqueV1(v));
        end
    end
end
%--------------------------------------------------------------------------
% store V2 
if Ops.useV2
    % FreqDisc data occationally have outcome of NaN - remove here. 
    nanix = find(isnan(D.uniqueV2));
  D.uniqueV2(nanix) = [];
    for v = 1:length(D.uniqueV2)
        if Ops.ByTrials
            v2StimTemp = Data.([s]).TrialData(bl).Stim_V2(D.stimIDX);
            D.V2_IDX{v} = find(v2StimTemp == D.uniqueV2(v));
        else
            D.V2_IDX{v} = find(Data.([s]).TrialData(bl).Stim_V2 == D.uniqueV2(v));
        end
    end
end
%-------------------------------------------------------------------------
if Ops.ByTrials 
    triNum = numel(find(Data.([s]).TrialData(bl).IsBlank ==0));
    for i = 1:size(Data.([s]).TrialData(bl).Trials,1); % loop through cells
        CatTrial = []; % concatenate trials
        ResidTrial = []; % concatenate residuals from trials
        for j = 1:size(Data.([s]).TrialData(bl).Trials,2) % loop through trials
            if Data.([s]).TrialData(bl).IsBlank(j) ==0
                temptrial = squeeze(Data.([s]).TrialData(bl).Trials(i,j,:)); % trace from single trial
                CatTrial = [CatTrial;temptrial]; % add single trace to concatenated experiment
                
                if Ops.Residuals
                    % find all trials that have same stim combo:
                    V1 = Data.([s]).TrialData(bl).Stim_V1(j);
                    V2 = Data.([s]).TrialData(bl).Stim_V2(j);
                    V1ix = find(Data.([s]).TrialData(bl).Stim_V1 == V1);
                    V2ix = find(Data.([s]).TrialData(bl).Stim_V2 == V2);
                    sameStim = intersect(V1ix,V2ix);
                    % average the traces from trials with same stim combo:
                    meanStim = mean(squeeze(Data.([s]).TrialData(bl).Trials(i,sameStim,:)),1);
                    % subtract mean from the single trial trace to generate
                    % residual trace
                    resid = temptrial - meanStim';
                    % concat residual traces
                    ResidTrial = [ResidTrial;resid];
                end
                
            end
        end
        
        ztraces(i,:) = CatTrial; 
        
        if Ops.Residuals
            ResidTraces(i,:) = ResidTrial;
        end
    end
else % non-trial based data start here::::::::::::::::::::::::::::::::::::::::::
   
    % -------------use motor activity as data to model: ------------
    if Ops.PredictMotor
        ztraces =  D.locoTrace ;
        if Ops.LocoBin
            ztraces = smoothdata(ztraces,'gaussian',20); % smooth to remove small blips within large bouts
            ztraces(ztraces>0) = 1;
        end
        
        % ------------------------------------------------------------------------
    else
        % non-trial based neuronal data start here :::::::::::::::::::
        if Ops.Residuals | Ops.CorrResid % non-trial based residuals
            % for non-trial based residuals, we create an array of zeros that
            % is the length of the timestamp of the experiment. We then find
            % the frames/timepoints that correspond to trials, and we input the
            % corresponding residual trace into those points in the array.
            for i = 1:size(Data.([s]).TrialData(bl).Trials,1); % loop through cells
                ResidTrace = zeros(size(t)); % create array of zeros
                count = 0;
                for j = 1:size(Data.([s]).TrialData(bl).Trials,2) % loop through trials
                    if Data.([s]).TrialData(bl).IsBlank(j) ==0 % skip over blank trials
                        count = count +1;
                        temptrial = squeeze(Data.([s]).TrialData(bl).Trials(i,j,:)); % single trial trace
                        % find matching stim and generate residual trace:
                        V1 = Data.([s]).TrialData(bl).Stim_V1(j);
                        V2 = Data.([s]).TrialData(bl).Stim_V2(j);
                        V1ix = find(Data.([s]).TrialData(bl).Stim_V1 == V1);
                        V2ix = find(Data.([s]).TrialData(bl).Stim_V2 == V2);
                        
                        sameStim = intersect(V1ix,V2ix);
                        meanStim = mean(squeeze(Data.([s]).TrialData(bl).Trials(i,sameStim,:)),1);
                        resid = temptrial - meanStim';
                        
                        % figure out where, in time, that residual trace should
                        % fall:
                        stime = D.Sound_Time(count);
                        [~, closest_frame_sound] = min(abs(t(:)-stime));
                        startFrame = closest_frame_sound - bsfr;
                        if startFrame<1
                            startFrame = 1; % edge case for rounding the first trial
                        end
                        
                        ResidTrace(startFrame:startFrame+length(resid)-1) = resid;
                    end
                end
                ResidTraces(i,:) = ResidTrace;
            end
        end
        
        % Get neuronal responses
        
    
            for i = 1:size(Data.([s]).TrialData(bl).F7,1)
                if length(Data.([s]).TrialData(bl).F7(i,:))>8000 % edge case - short data sets cant be detrended with window size
                    ftemp = locdetrend(Data.([s]).TrialData(bl).F7(i,:),30,[300 10]);
                else
                    ftemp = locdetrend(Data.([s]).TrialData(bl).F7(i,:),30,[100 10]);
                end
                F7 = zscore(ftemp);
                ztraces(i,:) = F7(:,:);
            end
       
    end
end

if ~Ops.Residuals & ~Ops.CorrResid
    ResidTraces = NaN;
end

% ----------SET UP PARAMS HERE!!!-------------------
params = GLM_simple_setRegs(Data.([s]),D,bl,ztraces,Ops,bsfr);


end