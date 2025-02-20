function visualize_extracted_data(ExtractedData, cellnum, loco_filter)
% Visualize the data on a cell by cell basic from ExtractedData variable
%
% Version history:
% - V1 = current version

%% Default inputs

if nargin < 2
    loco_filter = 0;
    cellnum = [];
elseif nargin < 3
    loco_filter = 0;
end

%% Apply activity filter

onlyPlotResponsive = 1;
plotRandomOrder = 1;
onlyPlotGroup = ''; %''
onlyPlotResponseType = ''; %''

if isempty(cellnum) %Options aren't relevant if cellnum is specified
    disp('------PARAMETERS-----')
    disp(['Only plot responsive cells: ' num2str(onlyPlotResponsive)])
    disp(['Plot in random order: ' num2str(plotRandomOrder)])
    disp(['Only plot group: ' onlyPlotGroup])
    disp(['Only plot response type: ' onlyPlotResponseType])
    disp('---------------------')
end

%% Apply loco filter

switch loco_filter
    case 0
        Loco = {'All'};
    case 1
        Loco = {'Running'};
    case 2
        Loco = {'NotRunning'};
    case 3
        Loco = {'All', 'Running', 'NotRunning'};
    otherwise
        error('Check loco_filter')
end

%% Make list of cells
%If user did not specify any cells, generate random list

if isempty(cellnum)
    nCells = size(ExtractedData.CellList,1);
    
    if plotRandomOrder
        cellnum = randperm(nCells);
    else
        cellnum = 1:nCells;
    end
else
    %If user specified cell(s), allow to plot even if not responsive
    onlyPlotResponsive = 0;
end

%% Get variables from ExtractedData

plot_figures = 1;
Ops = ExtractedData.Ops;
Ops.plot_figures = 1;
CellData = ExtractedData.CellData;
CellList = ExtractedData.CellList;
TrialData = ExtractedData.TrialData;
StimInfo = ExtractedData.StimInfo;
CellDataByStim = ExtractedData.CellDataByStim;
Traces = ExtractedData.FinalTraces;
Summary = ExtractedData.Summary;

%Accommodate for old ExtractedData formats (before GUI)
if isempty(Ops.smTime)
    Ops.smTime = 0;
end

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
   
for cc = 1:length(cellnum)
    c = cellnum(cc);

    if onlyPlotResponsive
        responsive = 0;
        for L = 1:length(Loco)
            if CellData.(Loco{L}).RF(c)
                responsive = 1;
            end
        end
        if ~responsive
            continue
        end
    end
    
    if ~isempty(onlyPlotGroup)
        if ~strcmp(CellList.Group(c), onlyPlotGroup)
            continue;
        end
    end
    
    if ~isempty(onlyPlotResponseType)
        if ~strcmp(Summary.(Loco{L}).ResponseType(c), onlyPlotResponseType)
            continue;
        end
    end
    
    disp(['Plotting ExtractedData row #' num2str(c)]);
    figure_names = {}; %Temporarily store figure titles
    figures = {}; %Temporarily store figures
    
    %% Detect peaks and troughs
    for L = 1:length(Loco)
        figure_names{L} = {strcat('ResponseByStim_', Loco{L})};
        PeakData = CellDataByStim.(Loco{L})(c).PeakData;
        [~, figures{L}] = simple_check_if_responsive(CellDataByStim.(Loco{L})(c).StimTracesAveraged, StimInfo.nBaselineFrames, StimInfo.fs, NaN, NaN, Ops.smTime, plot_figures, 'Data', PeakData, 'FigureDimensions', [length(StimInfo.V2),length(StimInfo.V1)], 'Subtitles', subtitles);
    end    
    
    %% Reliability
    if Ops.use_reliability && Ops.reliability_type == "Pearsons"             
        for L = 1:length(Loco)   
            if plot_figures
                figure_names{end+1} = strcat('Reliability_Traces_', Loco{L});
                figure_names{end+1} = strcat('Reliability_Blanks_', Loco{L});
                figure_names{end+1} = strcat('Reliability_Histograms_', Loco{L});
                [figures{end+1}, figures{end+2}, figures{end+3}] = plot_simple_reliability(CellDataByStim.(Loco{L})(c), StimInfo, Ops);
            end
        end
    end
    
    %% Xcorr
    if Ops.use_reliability && Ops.reliability_type == "XCorr"
        for L = 1:length(Loco)
            if plot_figures
                figure_names{end+1} = strcat('Xcorr_Traces_', Loco{L});
                figure_names{end+1} = strcat('Xcorr_Blanks_', Loco{L});
                figure_names{end+1} = strcat('Xcorr_Max_Histograms_', Loco{L});
                figure_names{end+1} = strcat('Xcorr_Zero_Histograms_', Loco{L});
                [figures{end+1}, figures{end+2}, figures{end+3}, figures{end+4}] = plot_simple_xcorr(CellDataByStim.(Loco{L})(c), StimInfo, Ops);
            end
        end
    end
    
    %% Reliability neighbors
    %This one might not make as much sense to visualize, because the "old" RF
    %will have already had neighbors incorporated during extract_data
    
    if Ops.use_neighbors
        for L = 1:length(Loco) 
            figure_names{end+1} = strcat('ReliabilityNeighbors_RF_', Loco{L});
            figure_names{end+1} = strcat('ReliabilityNeighbors_Traces_', Loco{L});
            
          
            if Ops.reliability_type == "Pearsons"  
                %Regenerate CellDataByStim.RF from RF and Reliability before neighbors
                TempForNeighbors = CellDataByStim.(Loco{L})(c);
                IsResponsive = [TempForNeighbors.PeakData.IsResponsive] == 1; %Convert NaNs to zeros
                
                if Ops.use_reliability
                    IsReliable = [TempForNeighbors.ReliabilityData.IsReliable]' == 1; %Convert NaNs to zeros
                    TempForNeighbors.RF = (IsResponsive + IsReliable) == 2; %Find overlap
                else
                    TempForNeighbors.RF = IsResponsive;
                end

                [New_RF, NeighborData, figures{end+1}, figures{end+2}] = simple_reliability_neighbors(TempForNeighbors, StimInfo, plot_figures, Ops);
            
            elseif Ops.reliability_type == "XCorr"
                %Regenerate CellDataByStim.RF from RF and Reliability before neighbors
                TempForNeighbors = CellDataByStim.(Loco{L})(c);
                IsResponsive = [TempForNeighbors.PeakData.IsResponsive] == 1; %Convert NaNs to zeros
                
                if Ops.use_reliability
                    IsReliable = TempForNeighbors.XcorrData.cc_max_z == 1; %Convert NaNs to zeros
                    TempForNeighbors.RF = (IsResponsive + IsReliable) == 2; %Find overlap
                else
                    TempForNeighbors.RF = IsResponsive;
                end
                
                [New_RF, NeighborData, figures{end+1}, figures{end+2}] = simple_xcorr_neighbors(TempForNeighbors, StimInfo, plot_figures, Ops);
            end
        end 
    end
                 
    %% Final trace
    for L = 1:length(Loco)
        figure_names{end+1} = strcat('FinalTrace_', Loco{L});
        PeakData = CellData.(Loco{L})(c,:);
        [~, figures{end+1}] = simple_check_if_responsive(Traces.(Loco{L})(c,:), StimInfo.nBaselineFrames, StimInfo.fs, NaN, NaN, Ops.smTime, plot_figures, 'Data', PeakData);
    end
    
    %% Plot receptive field heat maps using peak data
    if plot_figures
        for L = 1:length(Loco)
            figure_names{end+1} = strcat('RF_', Loco{L});
            [figures{end+1}] = simple_plot_PeakData(StimInfo, CellDataByStim.(Loco{L})(c).PeakData, 'RF', CellDataByStim.(Loco{L})(c).RF);
        end
    end

    %% Sparseness
    for L = 1:length(Loco)
        switch Ops.stim_protocol
            case {2 3 5 6}
                %Compute tuning sparseness
                figure_names{end+1} = strcat('Sparseness_', Loco{L});
                [SparsenessData, figures{end+1}] = simple_compute_tuning_sparseness(CellDataByStim.(Loco{L})(c).PeakData, CellDataByStim.(Loco{L})(c).RF, StimInfo, plot_figures);

                %Compute again using just zeros and ones for RF
                figure_names{end+1} = strcat('BinarySparseness_', Loco{L});
                [BinarySparsenessData, figures{end+1}] = simple_compute_tuning_sparseness(CellDataByStim.(Loco{L})(c).PeakData, CellDataByStim.(Loco{L})(c).RF, StimInfo, plot_figures, 'RF', CellDataByStim.(Loco{L})(c).RF);
        end
    end
    
    %% Stim-specific receptive field analyses
    for L = 1:length(Loco)
        figure_names{end+1} = strcat(StimInfo.StimType, Loco{L});

        switch Ops.stim_protocol
            
            case 2 %RF
                [TuningData, figures{end+1}] = simple_compute_frequency_tuning(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, plot_figures);
            
            case 3 %FM
                R = nan(size(CellDataByStim.(Loco{L})(c).RF));
                if Ops.use_reliability && Ops.reliability_type == "Pearsons"
                    R = [CellDataByStim.(Loco{L})(c).ReliabilityData.R]';
                end
                [FMData, figures{end+1}] = simple_compute_FM(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, R, StimInfo, Ops.plot_figures);
            
            case 5 %SAM
                [SAMData, figures{end+1}] = simple_compute_SAM(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, plot_figures);

            case 6 %SAMfreq
                [SAMFreqData, figures{end+1}] = simple_compute_SAMfreq(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, plot_figures);
                
            case 12 %SPONTANEOUS
                %Ca_detection pulls out spontaneous calcium transients and calculates their amplitude, frequency, and duration
                blockName = CellList.Block(c);
                b = find(strcmp([TrialData.Block], blockName));
                cellIndex = CellList.CellOrder(c);
                
                [~, SpontData, figures{end+1}] = simple_Ca_detection_spontaneous(TrialData(b).Timestamp, TrialData(b).F7, TrialData(b).Full_Loco, StimInfo.fs, cellIndex, Loco, plot_figures);

            case 13 %MARYSE BEHAVIOR STIM
                [BehaviorData, figures{end+1}] = simple_compute_tuning_maryse_behavior(CellDataByStim.(Loco{L})(c).RF, CellDataByStim.(Loco{L})(c).PeakData, StimInfo, plot_figures);
        end
    end
    
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
   
   %% Pause and wait for user input before plotting next cell
   uinput = input('Press enter to close figures and move on to the next cell, 5 to save figures and move on to next cell, 9 to quit: ');
   if uinput == 5
       %Figure out how many figures there were not including replicates for Loco
       nFigures = size(figures,2)/length(Loco);
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
   elseif uinput == 9
       return
   end
   close all
   
end
