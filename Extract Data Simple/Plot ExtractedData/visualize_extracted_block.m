function visualize_extracted_block(ExtractedData, loco_filter)
%% Visualize the average cell responses from each block in ExtractedData

plot_responsive_cells_only = 0;

%% Loco filter

if nargin < 2
    Loco = ExtractedData.Ops.loco_filter{1};
else
    switch loco_filter
        case 0
            Loco = 'All';
        case 1
            Loco = 'NotRunning';
        case 2
            Loco = 'Running';
    end
end

%%
StimInfo = ExtractedData.StimInfo;
Ops = ExtractedData.Ops;
Blocks = ExtractedData.BlockInfo.Block;
CellList = ExtractedData.CellList;
CellData = ExtractedData.CellData.(Loco);
IsRF = CellData.RF;
Traces = ExtractedData.FinalTraces.(Loco);

for bb = 1:length(Blocks)
    currentCellRows = find(strcmp(CellList.Block, Blocks{bb}));
    if plot_responsive_cells_only
        responsiveCellRows = intersect(currentCellRows,find(IsRF));
    else
        responsiveCellRows = currentCellRows;
    end
    responsiveCellNums = CellList.CellNumber(responsiveCellRows);
    responsiveCellData = CellData(responsiveCellRows,:);
    responsiveTraces = Traces(responsiveCellRows,:);
    
    if isempty(responsiveCellRows)
        continue
    end

    %Plot figure
    subtitles = strings(length(responsiveCellNums),1);
    for i = 1:length(subtitles)
        subtitles(i) = strcat('Cell' , {' '}, num2str(responsiveCellNums(i)), {' / Row '}, num2str(responsiveCellRows(i)));
    end
    binsize = 64;
    bins = ceil(size(responsiveTraces,1)/binsize);
    for b = 1:bins
        b1 = (b-1)*binsize + 1;
        b2 = b1 + binsize - 1;
        if b2 > size(responsiveTraces,1)
            b2 = size(responsiveTraces,1);
        end
        %Plot pre-existing PeakData
        PeakData = responsiveCellData(b1:b2,:);
        [~, ~] = simple_check_if_responsive(responsiveTraces(b1:b2,:), StimInfo.nBaselineFrames, StimInfo.fs, NaN, NaN, Ops.smTime, 1, 'Subtitles', subtitles(b1:b2), 'Data', PeakData);
        tak_suptitle(regexprep(Blocks{bb},'_',' ','emptymatch')) %Replace underscore with a space in plot title
        drawnow
    end
end

%For ephys data, each block corresponds to data from one particular
%stimulus; simple_extract_data is run for different cells and the code
%below ensures all blocks corresponding to one cell end up in the same subplot:  
if ExtractedData.StimInfo.StimProtocol == 17
    figlist = get(groot,'Children');
    newfig = figure;
    tcl = tiledlayout(newfig,'flow');
    
    for i = numel(figlist):-1:1
        figure(figlist(i));
        ax = gca;
        axcopy = copyobj(ax, tcl);
        axcopy.Layout.Tile = numel(figlist)-i+1;
        close(figlist(i)); % close the original figure
        
    end
end