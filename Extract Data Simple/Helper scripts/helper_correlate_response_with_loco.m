%% helper_correlate_response_with_loco
% Are cell sound-responses modulated by speed of locomotor activity?

%% Load all extracted data

savepath = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\New stim analysis\September 19th corrected BW and reliability\LocoTrials';

loadpath = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\New stim analysis\September 19th corrected BW and reliability';
cd(loadpath);

StimNames = {'FM', 'RF', 'SAM', 'SAMfreq', 'NoiseITI'};
ExtractedDataFiles = {'ExtractedData_FM_20230412-010119_newStimAnalysis.mat',...
    'ExtractedData_RF_20230412-031712_newStimAnalysis.mat',...
    'ExtractedData_SAM_20230412-031314_newStimAnalysis.mat',...
    'ExtractedData_SAMfreq_20230412-041442_newStimAnalysis.mat',...
    'ExtractedData_NoiseITI_20230412-001517.mat'};

S = cell(size(ExtractedDataFiles));
for f = 2:length(ExtractedDataFiles)
    load(ExtractedDataFiles{f});
    S{f} = ExtractedData;
end

%% Stack all trials from TrialData into one big sheet

%Columns to stack
Trials = table;

for f = 1:length(S)   
    
    %Subset data to get activated cells only
    SubsetOps = struct;
    SubsetOps.Loco                = 'All'; %All, Running, or NotRunning
    SubsetOps.ResponseType        = 'activated'; %No filtering
    SubsetOps.RF_Type             = ''; %No filtering
    SubsetOps.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
    SubsetOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1)
    SubsetOps.Group               = ''; %'' -> no filtering, or add group name here to only keep that group
    SubsetOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
    SubsetOps.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    SubsetOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    SubsetOps.SuppressOutput      = 1; %0 to print Ops to command line, 1 to suppress

    ExtractedData = simple_subset_ExtractedData(S{f}, SubsetOps);
    
    StimInfo = ExtractedData.StimInfo;
    Summary = ExtractedData.Summary;
    CellDataByStim = ExtractedData.CellDataByStim;
    TrialData = ExtractedData.TrialData;
    TrialData = add_RF_mapping_to_TrialData(StimInfo, TrialData); %Store StimID for each block

    %Go through each cell
    IsBF = [];
    MaxLoco = [];
    MeanLoco = [];
    WinLoco = [];
    PeakWin = [];
    PeakAUC = [];
    CellID = [];
    MatchedRow = [];
    Group = strings(0,0);
    Stim = strings(0,0);
    Mouse = strings(0,0);
    FOV = strings(0,0);
    
    for c = 1:size(Summary,1)
        grp = Summary.Group(c);
        mouse = Summary.MouseID(c);
        fov = Summary.FOV(c);
        matchrow = Summary.MatchedRow(c);
        cellnumber = Summary.CellNumber(c);
        cellorder = Summary.CellOrder(c);
        block = Summary.Block(c);
        isRF = find(CellDataByStim(c).RF);
        isBF = CellDataByStim(c).PeakData.Peak_AUC(isRF) == max(CellDataByStim(c).PeakData.Peak_AUC(isRF));
        
        %Go to trial data and only save trials for isRF stim
        block_ind = strcmp({TrialData.Block},block);
        stimID = TrialData(block_ind).StimID;        
        keeptrial = ismember(stimID,isRF);
        
        temploco = TrialData(block_ind).Loco(keeptrial,:);
        tempmaxloco = max(temploco,[],2);
        tempmeanloco = mean(temploco,2,'omitnan');
        tempwinloco = mean(temploco(:,16:46),2,'omitnan');
        tempPeak_win = TrialData(block_ind).Peak_win(cellorder,keeptrial)';
        tempPeak_AUC = TrialData(block_ind).Peak_AUC(cellorder,keeptrial)';
        tempisBF = stimID(keeptrial) == isRF(isBF);
        [tempgroup, tempstim, tempmouse, tempfov] = deal(strings(size(tempmaxloco)));
        tempgroup(:) = deal(grp);
        tempstim(:) = deal(StimNames{f});
        tempmouse(:) = deal(mouse);
        tempfov(:) = deal(fov);
            
        MaxLoco = [MaxLoco; tempmaxloco];
        MeanLoco = [MeanLoco; tempmeanloco];
        WinLoco = [WinLoco; tempwinloco];
        PeakWin = [PeakWin; tempPeak_win];
        PeakAUC = [PeakAUC; tempPeak_AUC];
        IsBF = [IsBF; tempisBF];
        CellID = [CellID; zeros(size(tempmaxloco))+cellnumber];
        MatchedRow = [MatchedRow; zeros(size(tempmaxloco))+matchrow];
        Group = [Group; tempgroup];
        Stim = [Stim; tempstim];
        Mouse = [Mouse; tempmouse];
        FOV = [FOV; tempfov];
    end
    
    stimtable = table(Stim, Group, Mouse, FOV, CellID, MatchedRow, IsBF, MaxLoco, MeanLoco, WinLoco, PeakWin, PeakAUC);
    Trials = [Trials; stimtable];
end

%% Save spreadsheet and Trials

cd(savepath)

save('Trials.mat', 'Trials');
writetable(Trials,'Trials.csv','Delimiter',',');
