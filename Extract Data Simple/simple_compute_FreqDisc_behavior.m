function [FreqDisc] = simple_compute_FreqDisc_behavior(LocoState,CellNum, TrialData, StimInfo,CellList)

% TODO: Set Default variables:
% Check number of inputs.



% CellList = ExtractedData.CellList;
% TrialData = ExtractedData.TrialData;
% LocoState = 'All';
% CellNum = 2;

minnumtrials = 5;
numshuff = 1000;
plotfig = 0;
%% loop through each mouse and calculate AUC
% FreqDisc = table;

% Find trial data that correspond to the cell/channel of interest:
for tb = 1:length(TrialData)
    stringBlocks(tb) = string(TrialData(tb).Block);
end

% row in TrialData that contains cell of interest
rowOI = find(strcmp(stringBlocks,CellList.Block(CellNum)));
cellIDX = CellList.CellOrder(CellNum); % index of cell/channel in trial data


% calculate trial discriminability
trials = squeeze(TrialData(rowOI).Trials(cellIDX,:,:));

if all(isnan(trials(:))) % if there's no data in this channel, every field in FreqDisc is NaN
    
    FreqDisc.z_test_blank_HitFA = NaN;
    FreqDisc.z_P_blank_HitFA = NaN;
    FreqDisc.z_CI_blank_HitFA = NaN;
    FreqDisc.z_R_blank_HitFA = NaN;
    FreqDisc.AUC_HitFA= NaN;
    FreqDisc.DataForFig_HitFA= NaN;
    FreqDisc.z_test_blank_CatchFA = NaN;
    FreqDisc.z_P_blank_CatchFA = NaN;
    FreqDisc.z_CI_blank_CatchFA = NaN;
    FreqDisc.z_R_blank_CatchFA = NaN;
    FreqDisc.AUC_CatchFA = NaN;
    FreqDisc.DataForFig_CatchFA = NaN;
    FreqDisc.z_test_blank_HitCatch = NaN;
    FreqDisc.z_P_blank_HitCatch = NaN;
    FreqDisc.z_CI_blank_HitCatch = NaN;
    FreqDisc.z_R_blank_HitCatch = NaN;
    FreqDisc.AUC_HitCatch = NaN;
    FreqDisc.DataForFig_HitCatch = NaN;
    
else
    % index trials by loco condition:
    if contains(LocoState , 'All')
        locIX = ones(size(TrialData(rowOI).IsRunning));
    elseif contains(LocoState , 'Running')
        locIX = TrialData(rowOI).IsRunning;
    else
        locIX = ~(TrialData(rowOI).IsRunning);
    end
    
    
    % find the index for the outcome numbers
    H = TrialData(rowOI).Stim_V2 == 1;
    F = TrialData(rowOI).Stim_V2 == 4;
    C = TrialData(rowOI).Stim_V2 == 5;
    
    Hh = false(size(H));
    Hh(H(:) == 1 & locIX(:) ==1) = true;
    Hit = find(Hh);
    
    Ff = false(size(F));
    Ff(F(:) ==1 & locIX(:) ==1) = true;
    FA = find(Ff);
    
    Cc = false(size(C));
    Cc(C(:) ==1 & locIX(:) ==1) = true;
    Catch = find(Cc);
    
    
    % ROC
    for r = 1:3
        if r == 1 % Hit - FA
            if sum(Hit) > minnumtrials  && sum(FA) > minnumtrials
                [z_test_blank,z_P_blank,z_CI_blank,z_R_blank, AUC,DataForFig] = trial_discrimination(Hit, FA, trials, numshuff,plotfig);
                FreqDisc.z_test_blank_HitFA = z_test_blank;
                FreqDisc.z_P_blank_HitFA = z_P_blank;
                FreqDisc.z_CI_blank_HitFA = z_CI_blank;
                FreqDisc.z_R_blank_HitFA = z_R_blank;
                FreqDisc.AUC_HitFA= AUC;
                FreqDisc.DataForFig_HitFA= DataForFig;
            else
                FreqDisc.z_test_blank_HitFA = NaN;
                FreqDisc.z_P_blank_HitFA = NaN;
                FreqDisc.z_CI_blank_HitFA = NaN;
                FreqDisc.z_R_blank_HitFA = NaN;
                FreqDisc.AUC_HitFA= NaN;
                FreqDisc.DataForFig_HitFA= NaN;
            end
            
        elseif r == 2 % Catch - FA
            if sum(Catch) > minnumtrials  && sum(FA) > minnumtrials
                [z_test_blank,z_P_blank,z_CI_blank,z_R_blank, AUC,DataForFig] = trial_discrimination(Catch, FA, trials, numshuff,plotfig);
                FreqDisc.z_test_blank_CatchFA = z_test_blank;
                FreqDisc.z_P_blank_CatchFA = z_P_blank;
                FreqDisc.z_CI_blank_CatchFA = z_CI_blank;
                FreqDisc.z_R_blank_CatchFA = z_R_blank;
                FreqDisc.AUC_CatchFA= AUC;
                FreqDisc.DataForFig_CatchFA= DataForFig;
                
            else
                FreqDisc.z_test_blank_CatchFA = NaN;
                FreqDisc.z_P_blank_CatchFA = NaN;
                FreqDisc.z_CI_blank_CatchFA = NaN;
                FreqDisc.z_R_blank_CatchFA = NaN;
                FreqDisc.AUC_CatchFA = NaN;
                FreqDisc.DataForFig_CatchFA = NaN;
            end
            
        else  % Hit - Catch
            if sum(Hit) > minnumtrials  && sum(Catch) > minnumtrials
                [z_test_blank,z_P_blank,z_CI_blank,z_R_blank, AUC,DataForFig] = trial_discrimination(Hit, Catch, trials, numshuff,plotfig);
                FreqDisc.z_test_blank_HitCatch = z_test_blank;
                FreqDisc.z_P_blank_HitCatch = z_P_blank;
                FreqDisc.z_CI_blank_HitCatch = z_CI_blank;
                FreqDisc.z_R_blank_HitCatch = z_R_blank;
                FreqDisc.AUC_HitCatch= AUC;
                FreqDisc.DataForFig_HitCatch= DataForFig;
                
            else
                FreqDisc.z_test_blank_HitCatch = NaN;
                FreqDisc.z_P_blank_HitCatch = NaN;
                FreqDisc.z_CI_blank_HitCatch = NaN;
                FreqDisc.z_R_blank_HitCatch = NaN;
                FreqDisc.AUC_HitCatch = NaN;
                FreqDisc.DataForFig_HitCatch = NaN;
            end % end of check for Hit -Catch
        end % r == 1 or 2 or 3
    end % for r = 1:4
end % end if trial were not all nans
end %end of function

