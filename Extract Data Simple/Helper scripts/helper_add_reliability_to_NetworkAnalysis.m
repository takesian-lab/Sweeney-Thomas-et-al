%helper_add_reliability_to_NetworkAnalysis

%% Load files

%Load ExtractedData & Correlations table

stim = {'FM', 'NoiseITI', 'RF', 'SAMfreq', 'SAMNoise'};
files = {'ExtractedData_FM_20230412-010119_newStimAnalysis.mat',...
    'ExtractedData_NoiseITI_20230412-001517.mat',...
    'ExtractedData_RF_20230412-031712_newStimAnalysis.mat',...
    'ExtractedData_SAM_20230412-031314_newStimAnalysis.mat',...
    'ExtractedData_SAMfreq_20230412-041442_newStimAnalysis.mat'};
    
tables = {'Correlations_Table_FM.mat',...
    'Correlations_Table_NoiseITI.mat',...
    'Correlations_Table_RF.mat',...
    'Correlations_Table_SAMNoise.mat',...
    'Correlations_Table_SAMfreq.mat'};  

extracted_datapath = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\New stim analysis\September 19th corrected BW and reliability';
network_datapath = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\Network Analysis_NoiseCorr\With_Response_Type_and_Reliability'; 

storeTables = {};

for i = 1:length(files)
    
    %Load extracted data
    cd(extracted_datapath)
    disp(files{i})
    load(files{i})
    
    %Load correlations
    cd(network_datapath)
    disp(tables{i})
    load(tables{i})
    
    %% Loop through correlations table and find cell 1 and cell 2 reliability in ExtractedData

    Correlations_Matrix_All.Cell1_R = nan(size(Correlations_Matrix_All,1),1);
    Correlations_Matrix_All.Cell2_R = nan(size(Correlations_Matrix_All,1),1);
    Correlations_Matrix_All.Average_R = nan(size(Correlations_Matrix_All,1),1);
    Correlations_Matrix_All.Cell1_R_matchRF = nan(size(Correlations_Matrix_All,1),1);
    Correlations_Matrix_All.Cell2_R_matchRF = nan(size(Correlations_Matrix_All,1),1);
    Correlations_Matrix_All.Average_R_matchRF = nan(size(Correlations_Matrix_All,1),1);

    Loco = unique(Correlations_Matrix_All.LocoType);

    for L = 1:length(Loco)
        locoInd = strcmp(Correlations_Matrix_All.LocoType,Loco{L});
        tempCorr = Correlations_Matrix_All(locoInd,:);

        Summary = ExtractedData.Summary.(Loco{L});
        CellDataByStim = ExtractedData.CellDataByStim.(Loco{L});

        nPairs = size(tempCorr,1);

        for p = 1:nPairs
            block = tempCorr.BlockName(p);
            cell1 = tempCorr.Cell1_s2P(p);
            cell2 = tempCorr.Cell2_s2P(p);

            if isnan(cell1) || isnan(cell2)
                continue;
            end

            blockind = strcmp(Summary.Block, block);
            cell1ind = Summary.CellNumber == cell1;
            cell2ind = Summary.CellNumber == cell2;

            cell1row = find((blockind+cell1ind)==2);
            cell2row = find((blockind+cell2ind)==2);

            if isempty(cell1row) || isempty(cell2row)
                error('troubleshoot')
            end

            if length(cell1row) > 1 || length(cell2row) > 1
                error('troubleshoot')
            end

            %Find matched RF (make sure not to include nan)
            cell1RF = CellDataByStim(cell1row).RF;
            cell2RF = CellDataByStim(cell2row).RF;
            matchRF = (cell1RF+cell2RF)>0;
            
            tempCorr.Cell1_R(p) = Summary.R2(cell1row); %R2 is reliability of RF=1 stimuli
            tempCorr.Cell2_R(p) = Summary.R2(cell2row);
            tempCorr.Average_R(p) = mean([tempCorr.Cell1_R(p),tempCorr.Cell2_R(p)],'omitnan');
            tempCorr.Cell1_R_matchRF(p) = mean([CellDataByStim(cell1row).ReliabilityData(matchRF).R],'omitnan');
            tempCorr.Cell2_R_matchRF(p) =  mean([CellDataByStim(cell2row).ReliabilityData(matchRF).R],'omitnan');
            tempCorr.Average_R_matchRF(p) = mean([tempCorr.Cell1_R_matchRF(p),tempCorr.Cell2_R_matchRF(p)],'omitnan');
        end

        Correlations_Matrix_All.Cell1_R(locoInd) = tempCorr.Cell1_R;
        Correlations_Matrix_All.Cell2_R(locoInd) = tempCorr.Cell2_R;
        Correlations_Matrix_All.Average_R(locoInd) = tempCorr.Average_R;
        Correlations_Matrix_All.Cell1_R_matchRF(locoInd) = tempCorr.Cell1_R_matchRF;
        Correlations_Matrix_All.Cell2_R_matchRF(locoInd) = tempCorr.Cell2_R_matchRF;
        Correlations_Matrix_All.Average_R_matchRF(locoInd) = tempCorr.Average_R_matchRF;
    end

    %SAVE
    storeTables{i} = Correlations_Matrix_All;
    save(tables{i}, 'Correlations_Matrix_All', '-v7.3');
    writetable(Correlations_Matrix_All, char(strcat(tables{i}(1:end-4), '.csv')),'Delimiter',',');

end

%% Make compiled table

%vertically combine
Compiled_Table = table;
for i=1:length(storeTables)
    Compiled_Table = vertcat(Compiled_Table,storeTables{1,i});
end

% Add columns to table indicating + and - correlations
for r = 1:6
    if r == 1
        T = Compiled_Table.TraceMaxCorrelations;
        SigT = Compiled_Table.TraceMaxZTest;
    elseif r == 2
        T = Compiled_Table.TraceZeroCorrelations;
        SigT = Compiled_Table.TraceZeroZTest;
    elseif r == 3
        T = Compiled_Table.NoiseCorrelations;
        SigT = Compiled_Table.NoiseZTest;
    elseif r == 4
        T = Compiled_Table.SignalCorrelations;
        SigT = Compiled_Table.SignalZTest;
    elseif r == 5
        T = Compiled_Table.NoiseXCorrelations;
        SigT = Compiled_Table.NoiseXZTest;
    elseif r == 6
        T = Compiled_Table.SignalXCorrelations;
        SigT = Compiled_Table.SignalXZTest;
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
        Compiled_Table.PosOrNegTraceMaxCorr = pos_neg;
        Compiled_Table.SigPosOrNegTraceMaxCorr = sig;
    elseif r == 2
        Compiled_Table.PosOrNegTraceZeroCorr = pos_neg;
        Compiled_Table.SigPosOrNegTraceZeroCorr = sig;
    elseif r == 3
        Compiled_Table.PosOrNegNoiseCorr = pos_neg;
        Compiled_Table.SigPosOrNegNoiseCorr = sig;
    elseif r == 4
        Compiled_Table.PosOrNegSignalCorr = pos_neg;
        Compiled_Table.SigPosOrNegSignalCorr = sig;
    elseif r == 5
        Compiled_Table.PosOrNegNoiseXCorr = pos_neg;
        Compiled_Table.SigPosOrNegNoiseXCorr = sig;
    elseif r == 6
        Compiled_Table.PosOrNegSignalXCorr = pos_neg;
        Compiled_Table.SigPosOrNegSignalXCorr = sig;
    end

end

%Remove NaNs (from blocks where there was only 1 sound-resp cell)
nanind = (isnan(Compiled_Table.Cell1_Match) + isnan(Compiled_Table.Cell2_Match)) > 0;
Compiled_Table(nanind,:) = [];

%Create unique matched pair labels
Cell1_Match_Sorted = Compiled_Table.Cell1_Match;
Cell2_Match_Sorted = Compiled_Table.Cell2_Match;
C2first = Cell1_Match_Sorted > Cell2_Match_Sorted; %wherever cel11 ID > cell2 ID, swap labels
Cell1_Match_Sorted(C2first) = Compiled_Table.Cell2_Match(C2first);
Cell2_Match_Sorted(C2first) = Compiled_Table.Cell1_Match(C2first);

%Combine match labels into a unique string and store
CombineMatch = strings(size(Cell1_Match_Sorted)); 
for i = 1:length(CombineMatch)
    CombineMatch(i) = strcat(num2str(Cell1_Match_Sorted(i)), {'_'}, num2str(Cell2_Match_Sorted(i)));
end
Compiled_Table.CombineMatch = CombineMatch;

save('CompiledTable.mat', 'Compiled_Table', '-v7.3');
writetable(Compiled_Table,'CompiledTable.csv','Delimiter',',');

%% Average all stim together
%Load CompiledTable.mat if you want to start from here

tic

%Preallocate table
Loco = unique(Compiled_Table.LocoType); %Sort by Locomotion
pairs_by_loco = cell(1,length(Loco)); %find all matched numbers per loco condition
nRows = 0;
for L = 1:length(Loco)
    pairs_by_loco{L} = unique(Compiled_Table.CombineMatch(strcmp(Compiled_Table.LocoType,Loco{L})));
    nRows = nRows + length(pairs_by_loco{L});
end
columnnames = Compiled_Table.Properties.VariableNames;
columntypes = varfun(@class,Compiled_Table,'OutputFormat','cell');
sz = [nRows length(columnnames)];
T = table('Size',sz,'VariableTypes',columntypes,'VariableNames',columnnames);

count = 1; %combined pair count

for L = 1:length(Loco)
    
    currentPairs = pairs_by_loco{L};
    
    for p = 1:length(currentPairs)
        
        ind_pair = intersect(find(strcmp(Compiled_Table.CombineMatch, currentPairs(p))),find(strcmp(Compiled_Table.LocoType,Loco{L})));
        
        if isempty(ind_pair)
            continue;
        end
        
        % if pair exists
        T.MouseName(count)              = Compiled_Table.MouseName(ind_pair(1));
        T.BlockName(count)              = Compiled_Table.BlockName(ind_pair(1));
        T.FieldofView(count)            = Compiled_Table.FieldofView(ind_pair(1));
        T.StimName(count)               = "AverageAllStim";
        T.CellType(count)               = Compiled_Table.CellType(ind_pair(1));
        T.LocoType(count)               = Compiled_Table.LocoType(ind_pair(1));
        T.Cell1_s2P(count)              = Compiled_Table.Cell1_s2P(ind_pair(1)); %C1/C2 correspond to first found pair
        T.Cell2_s2P(count)              = Compiled_Table.Cell2_s2P(ind_pair(1)); %C1/C2 correspond to first found pair
        T.Cell1_Match(count)            = Compiled_Table.Cell1_Match(ind_pair(1)); %C1/C2 correspond to first found pair
        T.Cell2_Match(count)            = Compiled_Table.Cell2_Match(ind_pair(1)); %C1/C2 correspond to first found pair
        T.Cell1_RedCell(count)          = Compiled_Table.Cell1_RedCell(ind_pair(1)); %C1/C2 correspond to first found pair
        T.Cell2_RedCell(count)          = Compiled_Table.Cell2_RedCell(ind_pair(1)); %C1/C2 correspond to first found pair
        T.Distance(count)               = Compiled_Table.Distance(ind_pair(1));
        T.TraceMaxCorrelations(count)   = mean(Compiled_Table.TraceMaxCorrelations(ind_pair),'omitnan');
        T.TraceZeroCorrelations(count)  = mean(Compiled_Table.TraceZeroCorrelations(ind_pair),'omitnan');
        T.TraceMaxZTest(count)          = max(Compiled_Table.TraceMaxZTest(ind_pair));
        T.TraceZeroZTest(count)         = max(Compiled_Table.TraceZeroZTest(ind_pair));
        T.NoiseCorrelations(count)      = mean(Compiled_Table.NoiseCorrelations(ind_pair),'omitnan');
        T.SignalCorrelations(count)     = mean(Compiled_Table.SignalCorrelations(ind_pair),'omitnan');
        T.NoiseZTest(count)             = max(Compiled_Table.NoiseZTest(ind_pair));
        T.SignalZTest(count)            = max(Compiled_Table.SignalZTest(ind_pair));
        T.NoiseXCorrelations(count)     = mean(Compiled_Table.NoiseXCorrelations(ind_pair),'omitnan');
        T.SignalXCorrelations(count)    = mean(Compiled_Table.SignalXCorrelations(ind_pair),'omitnan');
        T.NoiseXZTest(count)            = max(Compiled_Table.NoiseXZTest(ind_pair));
        T.SignalXZTest(count)           = max(Compiled_Table.SignalXZTest(ind_pair));
        
        %SPECIFIC TO THIS SCRIPT
        T.Cell1_ResponseType(count)     = "NA";
        T.Cell2_ResponseType(count)     = "NA";
        T.Cell1_R(count)                = nan; %C1/C2 order could be flipped so it is meaningless to average single cell values
        T.Cell2_R(count)                = nan; %C1/C2 order could be flipped so it is meaningless to average single cell values
        T.Average_R(count)              = mean(Compiled_Table.Average_R(ind_pair),'omitnan');
        T.Cell1_R_matchRF(count)        = nan; %C1/C2 order could be flipped so it is meaningless to average single cell values
        T.Cell2_R_matchRF(count)        = nan; %C1/C2 order could be flipped so it is meaningless to average single cell values
        T.Average_R_matchRF(count)      = mean(Compiled_Table.Average_R_matchRF(ind_pair),'omitnan');
        
        %Significant pos/negative
        T.PosOrNegTraceMaxCorr(count)       = sign(T.TraceMaxCorrelations(count));
        T.SigPosOrNegTraceMaxCorr(count)    = nan;
        T.PosOrNegTraceZeroCorr(count)      = sign(T.TraceZeroCorrelations(count));
        T.SigPosOrNegTraceZeroCorr(count)   = nan;
        T.PosOrNegNoiseCorr(count)          = sign(T.NoiseCorrelations(count));
        T.SigPosOrNegNoiseCorr(count)       = nan;
        T.PosOrNegSignalCorr(count)         = sign(T.SignalCorrelations(count));
        T.SigPosOrNegSignalCorr(count)      = nan;
        T.PosOrNegNoiseXCorr(count)         = sign(T.NoiseXCorrelations(count));
        T.SigPosOrNegNoiseXCorr(count)      = nan;
        T.PosOrNegSignalXCorr(count)        = sign(T.SignalXCorrelations(count));
        T.SigPosOrNegSignalXCorr(count)     = nan;
        
        if T.TraceMaxZTest(count) == 1;     T.SigPosOrNegTraceMaxCorr(count)    = T.PosOrNegTraceMaxCorr(count);    end
        if T.TraceZeroZTest(count) == 1;    T.SigPosOrNegTraceZeroCorr(count)   = T.PosOrNegTraceZeroCorr(count);   end
        if T.NoiseZTest(count) == 1;        T.SigPosOrNegNoiseCorr(count)       = T.PosOrNegNoiseCorr(count);       end
        if T.SignalZTest(count) == 1;       T.SigPosOrNegSignalCorr(count)      = T.PosOrNegSignalCorr(count);      end
        if T.NoiseXZTest(count) == 1;       T.SigPosOrNegNoiseXCorr(count)      = T.PosOrNegNoiseXCorr(count);      end
        if T.SignalXZTest(count) == 1;      T.SigPosOrNegSignalXCorr(count)     = T.PosOrNegSignalXCorr(count);     end

        %Match label
        T.CombineMatch(count)           = Compiled_Table.CombineMatch(ind_pair(1));
        
        count = count+1;
    end
end

Compiled_Table_AllStim = T;

toc

%% Save

save('CompiledTable_AverageStim.mat', 'Compiled_Table_AllStim') ;
writetable(Compiled_Table_AllStim,'CompiledTable_AverageStim.csv','Delimiter',',');
