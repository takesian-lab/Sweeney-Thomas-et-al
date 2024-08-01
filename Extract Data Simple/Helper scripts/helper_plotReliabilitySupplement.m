% helper_plotReliabilitySupplement

%% Apply activity filter
%FM Activated: 348, 3040, 2927*
%FM Suppressed: 385, 170, 420

cellnum = []; 
Loco = {'All'}; %{'All', 'Running', 'Not Running'}

%% Apply activity filter

onlyPlotResponsive = 1;
plotRandomOrder = 1;
onlyPlotGroup = 'PYR'; %''
onlyPlotResponseType = ''; %''

if isempty(cellnum) %Options aren't relevant if cellnum is specified
    disp('------PARAMETERS-----')
    disp(['Only plot responsive cells: ' num2str(onlyPlotResponsive)])
    disp(['Plot in random order: ' num2str(plotRandomOrder)])
    disp(['Only plot group: ' onlyPlotGroup])
    disp(['Only plot response type: ' onlyPlotResponseType])
    disp('---------------------')
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
            
%             for v = 1:size(StimInfo.Combined,1)
%                 %Only compute reliability on stim conditions that pass Z-score threshold (use low Z_level for this reason)
%                 ReliabilityOps = Ops; %Duplicate Ops for use in simple_reliability function
%                 ReliabilityOps.ComputeReliabilityControls = 0;
%                 ReliabilityOps.PlotPearsonFigure = 1;
%                 
%                 %Compute reliability by comparing trials to either blank trials or shifted traces
%                 trials = CellDataByStim.(Loco{L})(c).StimTraces{v};
%                 blanks = CellDataByStim.(Loco{L})(c).BlankTraces; %If empty, script will automatically use shifted traces
%                 ReliabilityData = simple_reliability(trials, blanks, StimInfo, ReliabilityOps);
%             end
            
            if plot_figures
                figure_names{end+1} = strcat('Reliability_Traces_', Loco{L});
                figure_names{end+1} = strcat('Reliability_Blanks_', Loco{L});
                figure_names{end+1} = strcat('Reliability_Histograms_', Loco{L});
                [figures{end+1}, figures{end+2}, figures{end+3}] = plot_simple_reliability(CellDataByStim.(Loco{L})(c), StimInfo, Ops);
            end
        end
    end
    
    %% Full trace
    figure_names{end+1} = {strcat('FullTrace_', Loco{L})};
    blockname = CellList.Block(c);
    blockInd = strcmp({TrialData.Block},blockname);
    cellOrder = CellList.CellOrder(c);
    F7 = TrialData(blockInd).F7(cellOrder,:);
    timestamp = TrialData(blockInd).Timestamp;
    soundtime = TrialData(blockInd).Sound_Time;
    IsBlank = TrialData(blockInd).IsBlank;
    dff = (F7 - mean(F7, 2, 'omitnan'))./mean(F7, 2, 'omitnan'); %DF/F = (total-mean)/mean

    figures{end+1} = figure; hold on
    plot(timestamp,dff);
    vline(soundtime(IsBlank==1),'r')
    vline(soundtime(IsBlank==0),'g')
    xlim([0, timestamp(end)])
    ylabel('DF/F')
    xlabel('Seconds')
    


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
