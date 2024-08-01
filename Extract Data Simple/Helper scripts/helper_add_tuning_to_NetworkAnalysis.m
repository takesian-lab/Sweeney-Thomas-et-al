%helper_add_tuning_to_NetworkAnalysis

%% Load files

%Load ExtractedData & Correlations table
%ExtractedData refers to the RF file
%Correlations table refers to any stim table with cells corresponding to the RF file
%This relies on matched rows. If you do not have matched rows saved in correlations table, run helper_add_matches_to_ExtractedData

%Network correlations file
network_stim = 'Spontaneous'; %This can be any type of stim table
network_datapath = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\NxDxNpC061921F4\Analysis December 2023\Combined blocks\NetworkAnalysis';

%ExtractedData file
RF_filename = 'ExtractedData_RF_20231203-210054.mat'; %RF ExtractedData file
RF_datapath = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\NxDxNpC061921F4\Analysis December 2023\Combined blocks';

%Load correlations,  
cd(network_datapath)
network_filename = ['Correlations_Table_' network_stim '.mat'];
load(network_filename);

%Load ExtractedData files
cd(RF_datapath);
load(RF_filename)

%% First add sig/sign columns to table

% Add columns to table indicating + and - correlations
for r = 1:6
    if r == 1
        T = Correlations_Matrix_All.TraceMaxCorrelations;
        SigT = Correlations_Matrix_All.TraceMaxZTest;
    elseif r == 2
        T = Correlations_Matrix_All.TraceZeroCorrelations;
        SigT = Correlations_Matrix_All.TraceZeroZTest;
    elseif r == 3
        T = Correlations_Matrix_All.NoiseCorrelations;
        SigT = Correlations_Matrix_All.NoiseZTest;
    elseif r == 4
        T = Correlations_Matrix_All.SignalCorrelations;
        SigT = Correlations_Matrix_All.SignalZTest;
    elseif r == 5
        T = Correlations_Matrix_All.NoiseXCorrelations;
        SigT = Correlations_Matrix_All.NoiseXZTest;
    elseif r == 6
        T = Correlations_Matrix_All.SignalXCorrelations;
        SigT = Correlations_Matrix_All.SignalXZTest;
    end

    % columns for noise corr
    pos_neg = nan(length(T),1);
    neg_index = find(T<0);
    pos_neg(neg_index) = -1;
    pos_index = find(T>0);
    pos_neg(pos_index) = 1;

    sig_index = find(SigT==1);

    sig = nan(length(T),1);
    sig(intersect(neg_index,sig_index))=-1;
    sig(intersect(pos_index,sig_index))=1;

    if r == 1
        Correlations_Matrix_All.PosOrNegTraceMaxCorr = pos_neg;
        Correlations_Matrix_All.SigPosOrNegTraceMaxCorr = sig;
    elseif r == 2
        Correlations_Matrix_All.PosOrNegTraceZeroCorr = pos_neg;
        Correlations_Matrix_All.SigPosOrNegTraceZeroCorr = sig;
    elseif r == 3
        Correlations_Matrix_All.PosOrNegNoiseCorr = pos_neg;
        Correlations_Matrix_All.SigPosOrNegNoiseCorr = sig;
    elseif r == 4
        Correlations_Matrix_All.PosOrNegSignalCorr = pos_neg;
        Correlations_Matrix_All.SigPosOrNegSignalCorr = sig;
    elseif r == 5
        Correlations_Matrix_All.PosOrNegNoiseXCorr = pos_neg;
        Correlations_Matrix_All.SigPosOrNegNoiseXCorr = sig;
    elseif r == 6
        Correlations_Matrix_All.PosOrNegSignalXCorr = pos_neg;
        Correlations_Matrix_All.SigPosOrNegSignalXCorr = sig;
    end
end

%% Do we need match variable? If so, create unique CombineMatch column

if strcmp(network_stim, 'RF')
    MatchRequired = 0;
else
    MatchRequired = 1;
    if all(isnan(Correlations_Matrix_All.Cell1_Match))
        error('Go back and add matches to correlations table with code helper_add_matches_to_ExtractedData')
    end
end

if MatchRequired
    %Remove NaNs (from blocks where there was only 1 sound-resp cell)
    nanind = (isnan(Correlations_Matrix_All.Cell1_Match) + isnan(Correlations_Matrix_All.Cell2_Match)) > 0;
    Correlations_Matrix_All(nanind,:) = [];

    %Create unique matched pair labels
    Cell1_Match_Sorted = Correlations_Matrix_All.Cell1_Match;
    Cell2_Match_Sorted = Correlations_Matrix_All.Cell2_Match;
    C2first = Cell1_Match_Sorted > Cell2_Match_Sorted; %wherever cel11 ID > cell2 ID, swap labels
    Cell1_Match_Sorted(C2first) = Correlations_Matrix_All.Cell2_Match(C2first);
    Cell2_Match_Sorted(C2first) = Correlations_Matrix_All.Cell1_Match(C2first);

    %Combine match labels into a unique string and store
    CombineMatch = strings(size(Cell1_Match_Sorted)); 
    for i = 1:length(CombineMatch)
        CombineMatch(i) = strcat(num2str(Cell1_Match_Sorted(i)), {'_'}, num2str(Cell2_Match_Sorted(i)));
    end
    Correlations_Matrix_All.CombineMatch = CombineMatch;
end

%% Loop through correlations table and find cell 1 and cell 2 BF in ExtractedData

Correlations_Matrix_All.Cell1_BF = nan(size(Correlations_Matrix_All,1),1);
Correlations_Matrix_All.Cell2_BF = nan(size(Correlations_Matrix_All,1),1);
Correlations_Matrix_All.DeltaBF = nan(size(Correlations_Matrix_All,1),1);
Correlations_Matrix_All.Cell1_CF = nan(size(Correlations_Matrix_All,1),1);
Correlations_Matrix_All.Cell2_CF = nan(size(Correlations_Matrix_All,1),1);
Correlations_Matrix_All.DeltaCF = nan(size(Correlations_Matrix_All,1),1);
Correlations_Matrix_All.RFOverlap = nan(size(Correlations_Matrix_All,1),1);

Loco = unique(Correlations_Matrix_All.LocoType);

for L = 1:length(Loco)
    locoInd = strcmp(Correlations_Matrix_All.LocoType,Loco{L});
    tempCorr = Correlations_Matrix_All(locoInd,:);
    
    Summary = ExtractedData.Summary.(Loco{L});
    CellDataByStim = ExtractedData.CellDataByStim.(Loco{L});
    
    nPairs = size(tempCorr,1);
    tempCorr.RFOverlap = nan(nPairs,1);

    for p = 1:nPairs
        block = tempCorr.BlockName(p);
        mouse = tempCorr.MouseName(p);
        fov = tempCorr.FieldofView(p);
        cell1 = tempCorr.Cell1_s2P(p);
        cell2 = tempCorr.Cell2_s2P(p);
        match1 = tempCorr.Cell1_Match(p);
        match2 = tempCorr.Cell2_Match(p);

        if isnan(cell1) || isnan(cell2)
            continue;
        end
        
        if MatchRequired
            %mouseind zz
            cell1ind = Summary.MatchedRow == match1;
            cell2ind = Summary.MatchedRow == match2;
            cell1row = find(cell1ind);
            cell2row = find(cell2ind);
        else
            blockind = strcmp(Summary.Block, block);
            cell1ind = Summary.CellNumber == cell1;
            cell2ind = Summary.CellNumber == cell2;
            cell1row = find((blockind+cell1ind)==2);
            cell2row = find((blockind+cell2ind)==2);
        end
        
        if isempty(cell1row) || isempty(cell2row)
            if MatchRequired == 1
                %Matched rows present in the correlations table might not exist in the RF data
                continue;
            else
                error('troubleshoot')
            end
        end
        
        if length(cell1row) > 1 || length(cell2row) > 1
            error('troubleshoot')
        end
        
        tempCorr.Cell1_BF(p) = Summary.BF(cell1row);
        tempCorr.Cell2_BF(p) = Summary.BF(cell2row);
        tempCorr.Cell1_CF(p) = Summary.CF(cell1row);
        tempCorr.Cell2_CF(p) = Summary.CF(cell2row);
        
        %Compute percent RF overlap
        cell1RF = CellDataByStim(cell1row).RF;
        cell2RF = CellDataByStim(cell2row).RF;
        
        %Remove nans from RF (important for Running and NotRunning conditions)
        nanind = (isnan(cell1RF) + isnan(cell2RF)) > 0;
        cell1RF(nanind) = [];
        cell2RF(nanind) = [];
        
        tempCorr.RFOverlap(p) = sum((cell1RF+cell2RF)==2)/length(cell1RF);
    end
    
    %Convert BFs to kHz and then take log2 to convert to octaves
    BF1 = log2(tempCorr.Cell1_BF*1000); 
    BF2 = log2(tempCorr.Cell2_BF*1000);
    tempCorr.DeltaBF = abs(BF1-BF2);
    
    CF1 = log2(tempCorr.Cell1_CF*1000); 
    CF2 = log2(tempCorr.Cell2_CF*1000);
    tempCorr.DeltaCF = abs(CF1-CF2);
    
    Correlations_Matrix_All.Cell1_BF(locoInd) = tempCorr.Cell1_BF;
    Correlations_Matrix_All.Cell2_BF(locoInd) = tempCorr.Cell2_BF;
    Correlations_Matrix_All.DeltaBF(locoInd) = tempCorr.DeltaBF;
    Correlations_Matrix_All.Cell1_CF(locoInd) = tempCorr.Cell1_CF;
    Correlations_Matrix_All.Cell2_CF(locoInd) = tempCorr.Cell2_CF;
    Correlations_Matrix_All.DeltaCF(locoInd) = tempCorr.DeltaCF;
    Correlations_Matrix_All.RFOverlap(locoInd) = tempCorr.RFOverlap;
end
    
%% Optional: Add layer from summary (Maryse multiplane data)

if any(strcmp(Summary.Properties.VariableNames,'Layer'))

    Correlations_Matrix_All.Cell1_Layer = nan(size(Correlations_Matrix_All,1),1);
    Correlations_Matrix_All.Cell2_Layer = nan(size(Correlations_Matrix_All,1),1);
    Correlations_Matrix_All.CombineLayer = strings(size(Correlations_Matrix_All,1),1);

    Loco = unique(Correlations_Matrix_All.LocoType);

    for L = 1:length(Loco)
        locoInd = strcmp(Correlations_Matrix_All.LocoType,Loco{L});
        tempCorr = Correlations_Matrix_All(locoInd,:);

        Summary = ExtractedData.Summary.(Loco{L});

        nPairs = size(tempCorr,1);
        for p = 1:nPairs
            block = tempCorr.BlockName(p);
            cell1 = tempCorr.Cell1_s2P(p);
            cell2 = tempCorr.Cell2_s2P(p);
            match1 = tempCorr.Cell1_Match(p);
            match2 = tempCorr.Cell2_Match(p);

            if isnan(cell1) || isnan(cell2)
                continue;
            end

            if MatchRequired
                cell1ind = Summary.MatchedRow == match1;
                cell2ind = Summary.MatchedRow == match2;
                cell1row = find(cell1ind);
                cell2row = find(cell2ind);
            else
                blockind = strcmp(Summary.Block, block);
                cell1ind = Summary.CellNumber == cell1;
                cell2ind = Summary.CellNumber == cell2;
                cell1row = find((blockind+cell1ind)==2);
                cell2row = find((blockind+cell2ind)==2);
            end

            if isempty(cell1row) || isempty(cell2row)
                error('troubleshoot')
            end

            if length(cell1row) > 1 || length(cell2row) > 1
                error('troubleshoot')
            end

            tempCorr.Cell1_Layer(p) = Summary.Layer(cell1row);
            tempCorr.Cell2_Layer(p) = Summary.Layer(cell2row);
            
            %Order and store combined layer variable
            layers = sort([Summary.Layer(cell1row), Summary.Layer(cell2row)],'ascend');
            tempCorr.CombineLayer(p) = strcat(num2str(layers(1)), {'_'}, num2str(layers(2)));

        end
        
        Correlations_Matrix_All.Cell1_Layer(locoInd) = tempCorr.Cell1_Layer;
        Correlations_Matrix_All.Cell2_Layer(locoInd) = tempCorr.Cell2_Layer;
        Correlations_Matrix_All.CombineLayer(locoInd) = tempCorr.CombineLayer;
    end
end

%% Save

return %make user intentionally use this part

writetable(Correlations_Matrix_All,char(strcat('Correlations_Table_withTuning_', network_stim, '.csv')),'Delimiter',',');
save(char(strcat('Correlations_Table_withTuning_', network_stim, '.mat')), 'Correlations_Matrix_All', '-v7.3');
