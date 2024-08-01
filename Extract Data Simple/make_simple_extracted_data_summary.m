function Summary = make_simple_extracted_data_summary(ExtractedData)
%The point of this function is to make one final summary table with all of the information
%that a user would want for easy copy-pasting into their statistical software of choice
%Use writetable(T) to write any table variable to csv

Loco = ExtractedData.Ops.loco_filter;
StimInfo = ExtractedData.StimInfo;
CellList = ExtractedData.CellList;
CellData = ExtractedData.CellData;
StimAnalysis = ExtractedData.StimAnalysis;

%Remove x and y coords from CellList
if ismember('xcirc',CellList.Properties.VariableNames)
    CellList.xcirc = [];
    CellList.ycirc = [];
end
if ismember('xpix',CellList.Properties.VariableNames)
    CellList.xpix = [];
    CellList.ypix = [];
end

%Combine CellList, CellData, and StimAnalysis
Summary = struct;
for L = 1:length(Loco)
    
    %Some data will have an empty StimAnalysis field (e.g. noiseburst)
    if ~isempty(fieldnames(StimAnalysis)) 
        
        %Add sparseness by BestResponse for some protocols only
        switch ExtractedData.Ops.stim_protocol
            case 2 %RF
                BestCol = StimAnalysis.(Loco{L}).BF;
                BestRow = StimAnalysis.(Loco{L}).BF_I;
                compute_sparseness = 1;

            case 5 %SAM
                BestCol = StimAnalysis.(Loco{L}).BestRate;
                BestRow = StimAnalysis.(Loco{L}).BestDepth;
                compute_sparseness = 1;

            case 6 %SAMfreq
                BestCol = StimAnalysis.(Loco{L}).BF;
                BestRow = StimAnalysis.(Loco{L}).BestDepth;
                compute_sparseness = 1;

            otherwise
                compute_sparseness = 0;
        end

        %Check if V1_Sparseness and V2_Sparseness exist
        if ~any(strcmp(StimAnalysis.(Loco{L}).Properties.VariableNames,'V1_Sparesness'))
            compute_sparseness = 0;
        end
        
        if compute_sparseness
            StimAnalysis.(Loco{L}).V1_BestSparseness = nan(height(StimAnalysis.(Loco{L})),1);
            StimAnalysis.(Loco{L}).V2_BestSparseness = nan(height(StimAnalysis.(Loco{L})),1);
            StimAnalysis.(Loco{L}).Binary_V1_BestSparseness = nan(height(StimAnalysis.(Loco{L})),1);
            StimAnalysis.(Loco{L}).Binary_V2_BestSparseness = nan(height(StimAnalysis.(Loco{L})),1);

            %get index corresponding to BestRow and BestCol
            for c = 1:height(StimAnalysis.(Loco{L}))
                BestCol_ind = find(StimInfo.V1 == BestCol(c));
                BestRow_ind = find(StimInfo.V2 == BestRow(c));

                if ~isempty(BestCol_ind)
                    StimAnalysis.(Loco{L}).V1_BestSparseness(c) = StimAnalysis.(Loco{L}).V1_Sparseness{c}(BestRow_ind);
                    StimAnalysis.(Loco{L}).V2_BestSparseness(c) = StimAnalysis.(Loco{L}).V2_Sparseness{c}(BestCol_ind);
                    StimAnalysis.(Loco{L}).Binary_V1_BestSparseness(c) = StimAnalysis.(Loco{L}).Binary_V1_Sparseness{c}(BestRow_ind);
                    StimAnalysis.(Loco{L}).Binary_V2_BestSparseness(c) = StimAnalysis.(Loco{L}).Binary_V2_Sparseness{c}(BestCol_ind);
                end
            end
        end
    
        %Columns to remove from StimAnalysis
        columnsToRemove = {'Response', 'Abs_Response', 'V1_Sparseness', 'V2_Sparseness', 'Binary_V1_Sparseness', 'Binary_V2_Sparseness'};
        for c = 1:length(columnsToRemove)
            if strmatch(columnsToRemove{c}, StimAnalysis.(Loco{L}).Properties.VariableNames)
                StimAnalysis.(Loco{L}).(columnsToRemove{c}) = [];
            end
        end
        
        Summary.(Loco{L}) = [CellList, CellData.(Loco{L}), StimAnalysis.(Loco{L})];
    else
        Summary.(Loco{L}) = [CellList, CellData.(Loco{L})];
    end
end
