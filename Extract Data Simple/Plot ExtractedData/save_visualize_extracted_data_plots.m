%helper_save_visualize_extracted_data_plots

%% Options

save_folder = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice\RF examples\RF\NEW';
Loco = 'All';

%% Get variables from ExtractedData


cd(save_folder)

Ops = ExtractedData.Ops;
Ops.plot_figures = 1;
CellData = ExtractedData.CellData;
CellList = ExtractedData.CellList;
TrialData = ExtractedData.TrialData;
StimInfo = ExtractedData.StimInfo;
CellDataByStim = ExtractedData.CellDataByStim;
Traces = ExtractedData.FinalTraces;
Summary = ExtractedData.Summary;

%Generate subtitles for figures
if isfield(StimInfo, 'Labels')
    subtitles = StimInfo.Labels;
else
    subtitles = strings(size(StimInfo.Combined,1),1);
    for s = 1:length(subtitles)
        subtitles(s) = strcat(num2str(StimInfo.Combined(s,1)), StimInfo.Units{1}, {' '}, num2str(StimInfo.Combined(s,2)), StimInfo.Units{2});
    end
end
    
%Update reliability_type for older data
if ~isfield(Ops, 'reliability_type')
    Ops.reliability_type = "Pearsons"; %Accommodate old ExtractedData files
end

if ~isfield(Ops, 'adjust_neighbor_p_value')
    Ops.adjust_neighbor_p_value = 1; %Accommodate old ExtractedData files
end
%% Loop through selected cells
   
for c = 1:height(Summary.(Loco))

    if ~Summary.(Loco).RF(c) || Summary.(Loco).Binary_Sparseness(c) == 1 || ~strcmp(Summary.(Loco).ResponseType(c), 'activated')
        continue
    end
    
    disp(['Plotting ExtractedData row #' num2str(c)]);
    figure_names = {}; %Temporarily store figure titles
    figures = {}; %Temporarily store figures
    
    %% Plot average traces
    figure_names{1} = {strcat('ResponseByStim_', Loco)};
    PeakData = CellDataByStim.(Loco)(c).PeakData;
    [~, figures{1}] = simple_check_if_responsive(CellDataByStim.(Loco)(c).StimTracesAveraged, StimInfo.nBaselineFrames, StimInfo.fs, NaN, NaN, Ops.smTime, 1, 'SmoothFigure', 1, 'Data', PeakData, 'FigureDimensions', [length(StimInfo.V2),length(StimInfo.V1)], 'Subtitles', subtitles);  

    %% Plot receptive field heat maps using peak data
    figure_names{end+1} = strcat('RF_', Loco);
    [figures{end+1}] = simple_plot_frequency_tuning(CellDataByStim.(Loco)(c).RF, CellDataByStim.(Loco)(c).PeakData, StimInfo);

   %% Add title to all figures
    for f = 1:size(figures,2)
        if isempty(figures{1,f})
            continue
        end

        set(0, 'currentfigure', figures{1,f});
        blockName = char(CellList.Block(c));
        cellNumber = CellList.CellNumber(c);
        sgtitle(strcat(regexprep(blockName(1,10:end),'_',' ','emptymatch'), ' Cell ', num2str(cellNumber))) %Replace _ with a space
    end
   
   %Figure out how many figures there were not including replicates for Loco
   nFigures = size(figures,2);
   figure_numbers = sort(repmat(1:nFigures, 1, size(figures,2)/nFigures));

   for f = 1:size(figures,2)
       %Skip if figure is empty or has been closed
        if isempty(figures{1,f}) || ~isvalid(figures{1,f})
            continue
        end
        set(0, 'currentfigure', figures{1,f});
        h_name = strcat(CellList.MouseID(c), '_CellRow', num2str(c), '_ID', num2str(cellNumber), '_Fig', num2str(figure_numbers(f)), '_', figure_names{f});
        saveas(figures{1,f}, h_name, 'fig');
        saveas(figures{1,f}, h_name, 'jpg');
   end

   close all
   
end

disp('All finished!')