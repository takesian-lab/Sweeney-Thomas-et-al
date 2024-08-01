function [MotorCorrelations_Matrix_All] = add_spontaneous_data_to_network_analysis(ExtractedData, MotorCorrelations_Matrix_All)
% Add spontaneous data from ExtractedData.StimAnalysis to NetworkAnalysis tables

%Check that stim protocol == 12 for spontaneous data
if ExtractedData.Ops.stim_protocol ~= 12
    return
end

%Network analysis data may have less cells than ExtractedData because user
%can choose whether to include only sound-responsive cells or not

%Use data from loco type ALL unless it doesn't exist
if isfield(ExtractedData.StimAnalysis, 'All')
    StimAnalysis = ExtractedData.StimAnalysis.All;
else
    Loco = fields(ExtractedData.StimAnalysis);
    StimAnalysis = ExtractedData.StimAnalysis.(Loco{1});
end

%Go through each cell from Network Analysis and find corresponding cell in ExtractedData
nCells = height(MotorCorrelations_Matrix_All);

%Preallocate columns
columnsToAdd = {'N', 'tr_rate', 'tr_med', 'tr_mean', 'tr_STD', 'tr_duration_med', 'tr_duration_mean', 'tr_duration_STD'};

for c = 1:length(columnsToAdd)
    MotorCorrelations_Matrix_All.(columnsToAdd{c}) = nan(nCells, 1);
end

for i = 1:nCells
    
    MouseName = MotorCorrelations_Matrix_All.MotorMouseName(i);
    FOV = MotorCorrelations_Matrix_All.MotorFieldofView(i);
    CellNumber = MotorCorrelations_Matrix_All.Cell_s2P_list(i);
    MatchNumber = MotorCorrelations_Matrix_All.Matched_list(i);
    
    ind_mouse = strcmp(ExtractedData.CellList.MouseID, MouseName);
    ind_FOV = strcmp(ExtractedData.CellList.FOV, FOV);
    ind_cell = ExtractedData.CellList.CellNumber == CellNumber;
    ind_match = ExtractedData.CellList.MatchedRow == MatchNumber;
    ind_mat = [ind_mouse, ind_FOV, ind_cell, ind_match];
    
    ind = find(sum(ind_mat,2) == 4);
    
    if isempty(ind) || length(ind) > 1
        error('Could not find matching cell in ExtractedData')
    end
    
    %Store spont data in table
    for c = 1:length(columnsToAdd)
        MotorCorrelations_Matrix_All.(columnsToAdd{c})(i) = StimAnalysis.(columnsToAdd{c})(ind);
    end
end
    
end %end function