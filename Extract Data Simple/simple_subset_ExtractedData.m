
function SubsetData = simple_subset_ExtractedData(ExtractedData, UserOps)

% This function subsets your data based on the options at the top of the script
% You can use this to filter your ExtractedData and then pass it on to other plotting functions

%Arguments:
% - ExtractedData file from output of simple_extract_data
% - UserOps: supply any of the Ops set at the top of this script in advance

%Version control:
%V1: current version created 2/1/2023

%% Options: Supplied by user in UserOps
%Follow the Ops template below to create your own and use as input variable
%DO NOT MODIFY THE OPTIONS DIRECTLY HERE, OR YOU COULD AFFECT OTHERS USING THIS SCRIPT

Ops = struct;
Ops.Loco                = 'All'; %All, Running, or NotRunning
Ops.ResponseType        = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
Ops.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
Ops.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
Ops.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1), put ~ in front to keep everything excep that FOV
Ops.Group               = ''; %'' -> no filtering, or add group name here to only keep that group, put ~ in front to keep everything except that Group
Ops.Ensemble            = ''; %'' -> no filtering, or add ensemble name for Markpoints Experiments, put ~ in front to keep everything except that ensemble
Ops.Layer               = ''; %'' -> no filtering, or add layer name to only keep that layer, put ~ in front to keep everything except that layer
Ops.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
Ops.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
Ops.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
Ops.plotCorrelations	= 0; %0 = don't plot dependent variable correlations, 1 = plot
Ops.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
Ops.plotShadedErrorBars = 0; %0 = don't plot, 1 = plot shaded error bars
Ops.FibPhotChannel      = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
Ops.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress

%Check UserOps against Ops and fill missing parameters
if nargin == 2
    allfields = fields(Ops);
    for f = 1:length(allfields)
        if ~isfield(UserOps, allfields{f})
            UserOps.(allfields{f}) = Ops.(allfields{f});
        end
    end
    Ops = UserOps;
end

%Make sure that chosen Loco type is present in the data and switch if not found:
if ~ismember(Ops.Loco, ExtractedData.Ops.loco_filter)
    Ops.Loco = ExtractedData.Ops.loco_filter{1};
end

if ~Ops.SuppressOutput
    disp('-----PLOT SETTINGS-----')
    disp(['Loco: ' Ops.Loco])
    disp(['ResponseType: ' Ops.ResponseType])
    disp(['RF_Type: ' Ops.RF_Type])
    disp(['IsResponsive: ' num2str(Ops.IsResponsive)])
    disp(['FOV: ' Ops.FOV])
    disp(['Group: ' Ops.Group])
    disp(['Ensemble: ' Ops.Ensemble])
    disp(['Layer: ' Ops.Layer])
    disp(['Sort by GCaMP: ' num2str(Ops.sortbyGCAMP)])
    disp(['Sort by Condition: ' num2str(Ops.sortbyCondition)])
    disp(['Sort by red cell: ' num2str(Ops.sortbyRedCell)])
    disp(['Plot correlations: ' num2str(Ops.plotCorrelations)])
    disp(['Smooth rasters: ' num2str(Ops.smooth_rasters)])
    disp(['Plot shaded error bars: ' num2str(Ops.plotShadedErrorBars)])
    disp(['FibPhot Channel: ' num2str(Ops.FibPhotChannel)])
    disp('---------------------------------')
end

%% Get data from ExtractedData

%Add Summary table to ExtractedData if not there
if ~isfield(ExtractedData,'Summary')
    ExtractedData.Summary = make_simple_extracted_data_summary(ExtractedData);
else
    if isempty(fields(ExtractedData.Summary))
        ExtractedData.Summary = make_simple_extracted_data_summary(ExtractedData);
    end
end

BlockInfo = ExtractedData.BlockInfo;
TrialData = ExtractedData.TrialData;
CellList = ExtractedData.CellList;

%Figure out if ExtractedData has Loco structure (All, Running, NotRunng)
%Data that has already been subset once will not have it
testfields = fields(ExtractedData.CellData);
if any(strcmp(testfields,'All')) || any(strcmp(testfields,'NotRunning')) || any(strcmp(testfields,'Running'))
    CellData = ExtractedData.CellData.(Ops.Loco);
    CellDataByStim = ExtractedData.CellDataByStim.(Ops.Loco);
    FinalTraces = ExtractedData.FinalTraces.(Ops.Loco);
    Summary = ExtractedData.Summary.(Ops.Loco);
else
    CellData = ExtractedData.CellData;
    CellDataByStim = ExtractedData.CellDataByStim;
    FinalTraces = ExtractedData.FinalTraces;
    Summary = ExtractedData.Summary;
end

%In some old datasets, the ResponseType column did not match the RF column. Fix that here:
Summary.RF_Type(Summary.RF == 0) = "none";
Summary.ResponseType(Summary.RF == 0) = "none";
   
%NoiseITI, Air, H20, etc won't have StimAnalysis
if ~isempty(fields(ExtractedData.StimAnalysis))
    if isfield(ExtractedData.StimAnalysis,Ops.Loco)
        StimAnalysis = ExtractedData.StimAnalysis.(Ops.Loco);
    else
        StimAnalysis = ExtractedData.StimAnalysis;
    end
end

%% Sort by Group vs. GCaMP vs. Condition

if Ops.sortbyGCAMP == 0
    GroupList = Summary.Group;
    GroupListB = BlockInfo.Group;
    if Ops.sortbyCondition > 0
        GroupList = split(GroupList,'_');
        GroupList = GroupList(:,Ops.sortbyCondition);
        GroupListB = split(GroupListB,'_');
        GroupListB = GroupListB(:,Ops.sortbyCondition);
    end 
    
elseif Ops.sortbyGCAMP == 1
    GroupList = Summary.GCaMP;
    GroupListB = BlockInfo.GCaMP;
    
elseif Ops.sortbyGCAMP == 2
    GroupList = Summary.Group;
    GroupListB = BlockInfo.Group;
    if Ops.sortbyCondition > 0
        GroupList = split(GroupList,'_');
        GroupList = GroupList(:,Ops.sortbyCondition);
        GroupListB = split(GroupListB,'_');
        GroupListB = GroupListB(:,Ops.sortbyCondition);
    end 
    GroupList = strcat(GroupList, '_', Summary.GCaMP);
    GroupListB = strcat(GroupListB, '_', BlockInfo.GCaMP);
end

%Replace group column with new group classifications
Summary.Group = GroupList;
BlockInfo.Group = GroupListB;
CellList.Group = GroupList;

%Account for data with no Group
if any(ismissing(Summary.Group))
    Summary.Group(ismissing(Summary.Group)) = "none";
end

if any(ismissing(BlockInfo.Group))
    BlockInfo.Group(ismissing(BlockInfo.Group)) = "none";
end

if any(ismissing(CellList.Group))
    CellList.Group(ismissing(Summary.Group)) = "none";
end

%% Filter by Group [only keep specified group]

Ops.Group = char(Ops.Group);

if isempty(Ops.Group)
    keep = true(size(Summary.Group));
    keepB = true(size(BlockInfo.Group));
elseif Ops.Group(1) == '~' %Keep everything except this Group
    keep = ~strcmp(Summary.Group, Ops.Group(2:end));
    keepB = ~strcmp(BlockInfo.Group, Ops.Group(2:end));
else
    %Check if Group exists in data
    GroupList = Summary.Group;
    GroupListB = BlockInfo.Group;
    
    if ~any(strcmp(GroupList, Ops.Group))
        error('Check Ops.Group')
    end
    
    keep = strcmp(GroupList, Ops.Group);
    keepB = strcmp(GroupListB, Ops.Group);
end
    
if ~isempty(fields(ExtractedData.StimAnalysis))
    StimAnalysis = StimAnalysis(keep,:);
end
FinalTraces = FinalTraces(keep,:);
Summary = Summary(keep,:);
CellList = CellList(keep,:);
CellDataByStim = CellDataByStim(keep);
CellData = CellData(keep,:);
BlockInfo = BlockInfo(keepB,:);
TrialData = TrialData(keepB);

%% Filter by FOV [only keep specified FOV]

Ops.FOV = char(Ops.FOV);

if isempty(Ops.FOV)
    keep = true(size(Summary.FOV));
    keepB = true(size(BlockInfo.FOV));
elseif Ops.FOV(1) == '~' %Keep everything except this FOV
    keep = ~strcmp(Summary.FOV, Ops.FOV(2:end));
    keepB = ~strcmp(BlockInfo.FOV, Ops.FOV(2:end));
else
    %Check if FOV exists in data
    FOVList = Summary.FOV;
    FOVListB = BlockInfo.FOV;
    
    if ~any(strcmp(FOVList, Ops.FOV))
        error('Check Ops.FOV')
    end
    
    keep = strcmp(FOVList, Ops.FOV);
    keepB = strcmp(FOVListB, Ops.FOV);
end
    
if ~isempty(fields(ExtractedData.StimAnalysis))
    StimAnalysis = StimAnalysis(keep,:);
end
FinalTraces = FinalTraces(keep,:);
Summary = Summary(keep,:);
CellList = CellList(keep,:);
CellDataByStim = CellDataByStim(keep);
CellData = CellData(keep,:);
BlockInfo = BlockInfo(keepB,:);
TrialData = TrialData(keepB);

%% Filter by Markpoints Ensemble [only keep specified Ensemble]
% Only found in Markpoits datasets

Ops.Ensemble = char(Ops.Ensemble);

if any(ismember(BlockInfo.Properties.VariableNames,'Ensemble')) %Not all data will have Ensemble
    if isempty(Ops.Ensemble)
        keep = true(size(Summary.Ensemble));
        keepB = true(size(BlockInfo.Ensemble));
    elseif Ops.Ensemble(1) == '~' %Keep everything except this Ensemble
        keep = ~strcmp(Summary.Ensemble, Ops.Ensemble(2:end));
        keepB = ~strcmp(BlockInfo.Ensemble, Ops.Ensemble(2:end));
    else
        %Check if Ensemble exists in data
        EnsembleList = Summary.Ensemble;
        EnsembleListB = BlockInfo.Ensemble;

        if ~any(strcmp(EnsembleList, Ops.Ensemble))
            error('Check Ops.Ensemble')
        end

        keep = strcmp(EnsembleList, Ops.Ensemble);
        keepB = strcmp(EnsembleListB, Ops.Ensemble);
    end

    if ~isempty(fields(ExtractedData.StimAnalysis))
        StimAnalysis = StimAnalysis(keep,:);
    end
    FinalTraces = FinalTraces(keep,:);
    Summary = Summary(keep,:);
    CellList = CellList(keep,:);
    CellDataByStim = CellDataByStim(keep);
    CellData = CellData(keep,:);
    BlockInfo = BlockInfo(keepB,:);
    TrialData = TrialData(keepB);
end

%% Filter by Layer [only keep specified Layer]

Ops.Layer = char(Ops.Layer);

if any(ismember(CellList.Properties.VariableNames,'Layer')) %Not all data will have Layer
    if isempty(Ops.Layer)
        keep = true(size(Summary.Layer));
    elseif Ops.Layer(1) == '~' %Keep everything except this FOV
        keep = ~strcmp(Summary.Layer, Ops.Layer(2:end));
    else
        %Check if FOV exists in data
        LayerList = string(Summary.Layer);

        if ~any(strcmp(LayerList, Ops.Layer))
            error('Check Ops.Layer')
        end

        keep = strcmp(LayerList, Ops.Layer);
    end

    if ~isempty(fields(ExtractedData.StimAnalysis))
        StimAnalysis = StimAnalysis(keep,:);
    end
    FinalTraces = FinalTraces(keep,:);
    Summary = Summary(keep,:);
    CellList = CellList(keep,:);
    CellDataByStim = CellDataByStim(keep);
    CellData = CellData(keep,:);
end

%% Sort by FibPhot Channel or Red Cell

%Determine if FibPhot data
hasFibPhot = any(strcmp('FibPhot', Summary.AnalysisPath));
FibPhotChannels = ["blue", "green", "red", "all"];
FibPhotChannelToAnalyze = FibPhotChannels{Ops.FibPhotChannel};

if hasFibPhot
    if strcmp(FibPhotChannelToAnalyze, "all")
        keep = true(size(Summary.Channel));
    else
        keep = strcmp(Summary.Channel, FibPhotChannelToAnalyze);
    end
else
    redCell = Summary.RedCell;

    if Ops.sortbyRedCell == 0
        keep = true(size(redCell));
    elseif Ops.sortbyRedCell == 1
        keep = logical(redCell);
    elseif Ops.sortbyRedCell == 2
        keep = logical(~redCell);
    end
end

if ~isempty(fields(ExtractedData.StimAnalysis))
    StimAnalysis = StimAnalysis(keep,:);
end
FinalTraces = FinalTraces(keep,:);
Summary = Summary(keep,:);
CellList = CellList(keep,:);
CellDataByStim = CellDataByStim(keep);
CellData = CellData(keep,:);

%% Sort by responsive cells

if Ops.IsResponsive == 1
    isRF = logical(Summary.RF);
elseif Ops.IsResponsive == 0
    isRF = ~logical(Summary.RF);
elseif Ops.IsResponsive == 2
    isRF = true(size(Summary.RF));
end

if ~isempty(fields(ExtractedData.StimAnalysis))
    StimAnalysis = StimAnalysis(isRF,:);
end
FinalTraces = FinalTraces(isRF,:);
Summary = Summary(isRF,:);
CellList = CellList(isRF,:);
CellDataByStim = CellDataByStim(isRF);
CellData = CellData(isRF,:);

%% Sort by response type and/or RF type

[isResponseType, isRF_Type] = deal(ones(size(Summary,1),1));
peaktype = 'None';

%Response type
if ~isempty(Ops.ResponseType)
    if ismember(Ops.ResponseType, {'activated', 'prolonged'})
        isResponseType = strcmp(Ops.ResponseType, Summary.ResponseType);
        peaktype = 'Peak';
    elseif strcmp(Ops.ResponseType, 'excitatory')
        isResponseType = (strcmp('activated', Summary.ResponseType) + strcmp('prolonged', Summary.ResponseType)) == 1;
        peaktype = 'Peak';
    elseif strcmp(Ops.ResponseType, 'suppressed')
        isResponseType = strcmp(Ops.ResponseType, Summary.ResponseType);
        peaktype = 'Trough';
    elseif strcmp(Ops.ResponseType, 'none')
        isResponseType = strcmp(Ops.ResponseType, Summary.ResponseType);
        peaktype = 'None';
    else
        error('ResponseType not found')
    end
end

%RF type
if ~isempty(Ops.RF_Type)
    isRF_Type = strcmp(Ops.RF_Type, Summary.RF_Type);
    if strcmp(Summary.RF_Type, 'inhibitory')
        peaktype = 'Trough';
    elseif strcmp(Summary.RF_Type, 'excitatory')
        peaktype = 'Peak';
    end
end

isType = (isResponseType + isRF_Type) == 2;

if ~isempty(fields(ExtractedData.StimAnalysis))
    StimAnalysis = StimAnalysis(isType,:);
end

FinalTraces = FinalTraces(isType,:);
Summary = Summary(isType,:);
CellList = CellList(isType,:);
CellDataByStim = CellDataByStim(isType);
CellData = CellData(isType,:);

%% Make SubsetData structure to pass to functions
SubsetData = struct;
SubsetData.SubsetOps = Ops;
SubsetData.Date = ExtractedData.Date;
SubsetData.Ops = ExtractedData.Ops;
SubsetData.StimInfo = ExtractedData.StimInfo;
SubsetData.BlockInfo = BlockInfo;
SubsetData.TrialData = TrialData;
SubsetData.CellList = CellList;
SubsetData.CellData = CellData;
SubsetData.CellDataByStim = CellDataByStim;
if ~isempty(fields(ExtractedData.StimAnalysis))
    SubsetData.StimAnalysis = StimAnalysis;
else
    SubsetData.StimAnalysis = [];
end
SubsetData.FinalTraces = FinalTraces;
SubsetData.Summary = Summary;
SubsetData.Filename = ExtractedData.Filename;
SubsetData.peaktype = peaktype;
SubsetData.HasFibPhot = hasFibPhot;

