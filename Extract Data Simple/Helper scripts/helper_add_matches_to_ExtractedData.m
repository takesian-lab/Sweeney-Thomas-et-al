%helper_add_matches_to_ExtractedData
%add matches from match spreadsheet to ExtractedData after the fact

%% Load ExtractedData first

%%%%%%%%%%%%%% REPLACE MATCHING PATH %%%%%%%%%%%%%%%%
ExtractedData.Ops.matching_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\CompiledBlocks_Activation_PYR_2023\Recompile_CombinedBlocks\Matching_CombinedBlocks.xlsx';
ExtractedData.Ops.include_matches = true;

ExtractedData.CellList = add_matches_to_CellList(ExtractedData.CellList, ExtractedData.Ops.matching_path);
ExtractedData.Summary = make_simple_extracted_data_summary(ExtractedData);

%Save ExtractedData
cd(ExtractedData.Ops.save_path)
save([ExtractedData.Filename '.mat'], 'ExtractedData', '-v7.3');
Loco = ExtractedData.Ops.loco_filter;
for L = 1:length(Loco)
    writetable(ExtractedData.Summary.(Loco{L}),[ExtractedData.Filename '_' Loco{L} '_Summary.csv'],'Delimiter',',');
end

%% Add matches to NetworkAnalysis
%Go to network analysis path

filename = 'Correlations_Table_Loco.mat'; %%%%%%%%%%%%%% REPLACE FILENAME %%%%%%%%%%%%
load(filename)

%Find rows corresponding to each cell in corr mat
CellList = ExtractedData.CellList;
for i = 1:height(CellList)
    mouse = CellList.MouseID(i);
    fov = CellList.FOV(i);
    cell = CellList.CellNumber(i);
    match = CellList.MatchedRow(i);
    mouse_ind = strcmp(Correlations_Matrix_All.MouseName,mouse);
    fov_ind = strcmp(Correlations_Matrix_All.FieldofView,fov); 
    cell1_ind = Correlations_Matrix_All.Cell1_s2P == cell;
    cell2_ind = Correlations_Matrix_All.Cell2_s2P == cell;
    
    %Combine block and cell ind
    mouse_and_FOV_and_cell1_ind = (mouse_ind + fov_ind + cell1_ind) == 3;
    mouse_and_FOV_and_cell2_ind = (mouse_ind + fov_ind + cell2_ind) == 3;
    
    %Fill in matched row
    Correlations_Matrix_All.Cell1_Match(mouse_and_FOV_and_cell1_ind) = match;
    Correlations_Matrix_All.Cell2_Match(mouse_and_FOV_and_cell2_ind) = match;
end    

save(filename,'Correlations_Matrix_All');
writetable(Correlations_Matrix_All,[filename(1:end-4) '.csv'],'Delimiter',',');