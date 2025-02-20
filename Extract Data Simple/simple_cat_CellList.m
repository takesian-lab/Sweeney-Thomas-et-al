% update cell list during simple_Concat_Blocks
function    [tempCellList] = simple_cat_CellList(rows,BlockInfo,CellList);
catcellrow = [];
for i = 1:length(rows)
    blockname = BlockInfo.Block(rows(i));
    cellrows = find(CellList.Block == blockname);
    catcellrow = [catcellrow;cellrows];
    temprows{i} = cellrows;
    AnalysisPath(i) = CellList.AnalysisPath(temprows{i}(1));
end



% fibphot:
hasfibphot = any(AnalysisPath == 'FibPhot');
noFunc = ismissing(AnalysisPath);
if  hasfibphot && sum(strcmp(AnalysisPath,'FibPhot')) == length(rows) % all data are FP
    tempCellList = CellList(temprows{1},:);
elseif ~hasfibphot && sum(noFunc) ==0 %s2p
    
    if all(AnalysisPath == AnalysisPath(1))
        tempCellList = CellList(temprows{1},:);
    else
        
        question_ops = struct;
        question_ops.Default = 'Data are not from online analysis';
        question_ops.Interpreter = 'none';
        
        answer = questdlg('Analysis paths do not match, are these data from quick, online pipeline', 'Possible problem with concatenating data',...
            'Data are not from online analysis', 'Continue, these are online, quick analysis data',question_ops);
        
        switch answer
            case 'Data are not from online analysis'
                error('matching data required to concat S2P data')
            case 'Continue, these are online, quick analysis data'
                tempCellList = CellList(temprows{1},:);   
        end
    end
    
elseif sum(noFunc) == length(rows) % no functinal data
    tempCellList = CellList(temprows{1}(1),:);
else
    error('mixed data, maybe FP and no Functional trying to cat together')
end
end