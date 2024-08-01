function [block] = define_markpoints_experiment(block, varargin)
% Confirm details of markpoints experiment and targeted cells with user input
%
% Definitions
% - Target/Point - targeted position during markpoints experiment
% - Targeted cell - cell intended to be stimulated
% - Ensemble - group of cell(s) targeted during a single markpoints stimultation

% Description:
% - For each targeted ensemble, a figure will be plotted showing the target
%   markpoints positions, detected target cells, and their traces
% - The user will be prompted to identify sham ensembles (called 'Sham') if
%   they cannot be determined automatically, and to provide a name for each
%   non-sham ensemble. NOTE: Be careful to use identical names for corresponding
%   ensembles across blocks. This ensures they are processed together in
%   later stages of the analysis pipeline. If you make a mistake, you can
%   rerun this function on the block to redo it
% - If your experiment contains only one ensemble (other than Sham), OR two ensembles
%   that you want to be named the same think (ex. two 'None' ensembles), you can
%   record this in the Info spreadsheet in a column called "markpoints_experiment"
%   for it to be detected automaticaly. If the experiment only has a sham ensemble,
%   write 'Sham'. If it has one sham and one non-sham ensemble, write the non-sham
%   ensemble's name (e.g. 'Sound Responsive'). If you have more than one
%   non-sham ensemble, leave blank (unless you want them to both be named the same thing).
% - When complete, the script will save a new structure called 'MarkpointsExperiment' in 
%   the block. It will also save the figures and struct directly in the
%   block folder, where it can be uploaded again automatically the next
%   time you run this function.
% - It is recommended that YOU MUST CHECK THE SAVED FIGURES to confirm that
%   all the cells you intended to target have been found by Suite2p and
%   identified correctly as such. This can also help you determine whether your
%   activation worked, and catch any other problems right away.
%
% Argument(s): 
%   block (struct)
%   varargin:
%   'LoadPreviousMarkpointsExperiment' (0, 1, or 2)
%   - 0: (default) Don't look for previous MarkpointsExperiment. Use for compiling
%        block from scratch or writing over previous data
%   - 1: Load MarkpointsExperiment.mat from block folder (saves the user
%        from having to provide any input). Use when you know the previous data
%        is correct, but you have to recompile the block for some other reason
%   - 2: If MarkpointsExperiment already exists in the block, use that data
%        for determining ensemble names (so user doesn't have to retype
%        them). All other data will be recomputed from scratch (for use if
%        you need to add a new element to MarkpointsExperiments or correct the code)
%   'SaveFigPath'
%   - 'default': (default) figures will be saved in each block's folder
%   - path: save all figures in this designated path
%   'MaxTargetDistance'
%   - 20: (default) max distance in microns from each point to be considered a target cell
%   'UserOffset'
%   - [0 0]: (default) Do not perform x/y offset correction
%   - [   ]: Perform x/y offset correction by computing average Suite2p xoff/yoff values and iterating through each combo of +/- xy to find the
%            one that minimizes the distance between target locations and target cells within 20 microns
%   - [x y]: Input your own x/y offset correction
% 
% Returns:
%   block (struct)
%   block.MarkpointsExperiment (struct)
%   - block_filename
%   - RecordingType - Trial-Based (T-Series recording) or Continuous (BOT recording)
%   - ControlType - NotInThisBlock (Sham not present), ShutterClosed (Sham stim. with uncaging shutter closed), or
%     ControlPointsWithShutterOpen (Sham points in FOV with uncaging shutter open)
%   - UncagingShutter - Open, Closed, or Variable (what was the position
%     of the uncaging shutter for the entire experiment?)
%   - Multiplane - Was data multiplane
%   - HasTosca - Is there a corresponding Tosca recording?
%   - HasVoltageRecording - Is there a corresponding VoltageRecording?
%   - StartsWithDummyRow - Whether the first PraireView markpoints iteration
%     was skipped (This happens with BOT recordings and the code corrects for it)
%   - Ensemble - Ensemble name if only one non-Sham ensemble
%   - MarkPoints - Markpoints data from XML file
%   - StimInfo - from simple_prepare_variables
%   - StimTable - Markpoints information for every trial in the block
%     -- EnsembleName - Sham or user-provided name
%     -- TargetedCells - Suite2p cell numbers of detected target cells
%        NOTE: There can be fewer targeted cells than points if not all
%        cells were detected by suite2p, or if there was no cell within
%        MaxTargetDistance (15 microns) of all points
%     -- Distance - Distance in microns from each target cell to the nearest target point
%     -- PointIndex - Which target point (in X and Y) corresponds to each cell 
%     -- nCells - Number of targeted cells detected in ensemble
%     -- TrialType - 0 (sham), 1 (stim), 2 (sound + sham), 3 (sound + stim), 9 (uncaging shutter opening)
%     -- Stim_V1 - corresponding to StimInfo
%     -- Stim_V2 - corresponding to StimInfo
%     -- ToscaGroup - Tosca program Group name (used for interleaved experiments)
%     -- UncShutter - 1 for open, 0 for closed
%     -- IgnoreTrial - 1 for trials on which uncaging shutter was opening during trial-based experiments (discard from analysis)
%     -- Duration - Duration in ms of stimulation
%     -- Power - PrairieView power of 1040 laser during  stimulation
%     -- Revs - # spiral revolutions
%     -- PointNames - Markpoints ensemble name specified in PV Markpoints window
%     -- nPoints - Number of points targeted simultaneously
%     -- Size - Size of markpoints spiral stimulation in microns
%     -- X - X position in microns of each target point
%     -- Y - Y position in microns of each target point
%
% VERSION HISTORY:
% - v1: current version created 5/23/23

%% Skip if no MarkPoints data

close all; %So that we don't save pre-existing figures

if ~isfield(block.setup, 'hasMarkPoints')
    return
end

% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

% Load previous markpoints experiment: %0 = don't look, 1 = load from file, 2 = use existing block.markpoints (if it exists, e.g. during recompile) to find ensemble names
addParameter(p,'LoadPreviousMarkpointsExperiment', 0); 

% Where to save MarkpointsFigures? If default, they will be in each block folder, otherwise they will be saved in provided path
addParameter(p,'SaveFigPath', 'default'); 

% Maximum distance in microns of a cell from target position to be considered target
addParameter(p,'MaxTargetDistance', 20); 

% Apply Suite2p x/y offset to target positions? [0 0] to skip offset correction, [] to perform, otherwise input desired x/y offset,  ex [-2 1]
% (USE CAREFULLY: only use right now is for making sure x/y spatial control offset is the same)
addParameter(p,'UserOffset', []); 

parse(p, varargin{:});

ops = p.Results; 
% ------- end parse varargin

%% Accommodate multiplane data

%First case: block contains multiple planes
if isfield(block, 'MultiplaneData')
    warning('MarkPoints experiment not defined for multiplane data. It is recommended to compile as a single plane.')
    return
end

%Second case: multiplane data run as single plane (recommended to run this way to extract data more easily later!)
nPlanes = block.setup.XML.nPlanes;
if nPlanes > 1
    suite2p_path = char(block.setup.suite2p_path);
    if strcmp(suite2p_path(end-6:end),'suite2p')
        %Option 1: data converted to single plane data using
        %convert_multiplane_to_singleplane and path directed to suite2p folder
        planeTag = '';
    else
        %Option 2: data processed through suite2p as above OR as multiplane and
        %directed to suite2p/plane# folder
        planeTag = ['_' suite2p_path(end-5:end)]; %this is the plane# part of the path
    end
else
    planeTag = '';
end

%% Look for previously saved MarkpointsExperiment.mat file
%Every subsequent time, if this file is found we can load the details and skip this step

disp('Defining MarkPoints experiment...');

filepath = strcat(block.setup.block_path, '\', 'MarkpointsExperiment', planeTag, '.mat');
figurepath = strcat(block.setup.block_path, '\MarkpointsExperimentFigures', planeTag);

if (ops.LoadPreviousMarkpointsExperiment == 1) && isfile(filepath)
    disp('MarkpointsExperiment found');
    load(filepath);
    block.MarkpointsExperiment = MarkpointsExperiment;
    return
elseif ops.LoadPreviousMarkpointsExperiment == 2
    %Use current block.MarkpointsExperiment for ensemble name only!
    if isfield(block,'MarkpointsExperiment')
        prevStimTable = block.MarkpointsExperiment.StimTable;
    end
end

%% Make new MarkpointsExperiment struct to record experiment parameters

if ~isfield(block.setup, 'markpoints_ensemble')
    block.setup.markpoints_ensemble = nan;
end

MarkpointsExperiment = struct;
MarkpointsExperiment.block_filename = block.setup.block_filename;
MarkpointsExperiment.RecordingType = "";
MarkpointsExperiment.ControlType = "";
MarkpointsExperiment.UncagingShutter = "";
MarkpointsExperiment.Multiplane = block.setup.XML.nPlanes > 1;
MarkpointsExperiment.HasTosca = ~ismissing(block.setup.Tosca_path);
MarkpointsExperiment.HasVoltageRecording = ~any(ismissing(block.setup.VR_filename));
MarkpointsExperiment.StartsWithDummyRow = nan;
MarkpointsExperiment.Ensemble = block.setup.markpoints_ensemble;
MarkpointsExperiment.OffsetXY = [nan, nan];
MarkpointsExperiment.MarkPoints = block.setup.XML.markPoints;

%Determine recording type
if any(ismember(block.setup.XML.sequenceType, 'MarkPoints'))
    MarkpointsExperiment.RecordingType = "Trial-based";
else
    MarkpointsExperiment.RecordingType = "Continuous";
end

%Determine dummy row
if isfield(block.setup.XML.markPoints, 'startsWithDummyRow')
    MarkpointsExperiment.StartsWithDummyRow = block.setup.XML.markPoints.startsWithDummyRow;
end

%Determine uncaging shutter status and control type
uncagingShutter = block.parameters.uncagingShutter';
if all(uncagingShutter == 1)
    MarkpointsExperiment.UncagingShutter = "Open";
    MarkpointsExperiment.ControlType = "TBD";
elseif all(uncagingShutter == 0)
    MarkpointsExperiment.UncagingShutter = "Closed";
    MarkpointsExperiment.ControlType = "ShutterClosed";
else
    MarkpointsExperiment.UncagingShutter = "Variable";
    MarkpointsExperiment.ControlType = "ShutterClosed";
end

%Display parameters
disp(struct2table(MarkpointsExperiment));

%% Get data from block

%only include green cells in online BOT analysis
if isequal(block.setup.analysis_name, block.setup.block_name)
    keep = block.redcell == 0;
else
    keep = true(size(block.redcell));
end

cellNumbers = block.cell_number(keep);
redcell = block.redcell(keep);
nCells = length(cellNumbers);
stat = block.stat(keep);
Trials = block.aligned_stim.zscore(keep,:,:);
nTrials = length(block.Sound_Time);
baseline_length = block.setup.constant.baseline_length; %seconds
after_length = block.setup.constant.after_stim; %seconds
fs = ceil(block.setup.framerate);
nBaselineFrames = round(baseline_length*fs); %frames
xtick = 0:fs:ceil(baseline_length + after_length)*fs;
xticklabel = 0:1:ceil(baseline_length + after_length);
xlin = linspace(0, baseline_length+after_length, size(Trials,3));

%Tosca groups
if isfield(block, 'errorData') %Blocks with Tosca data only
    ToscaGroups = block.errorData.Group';
    ToscaGroups(block.errorData.error_trials) = [];
else
    ToscaGroups = strings(nTrials,1);
end

%Pixels to micron conversion
micronsPerPixelX = str2double(block.setup.XML.micronsPerPixelX);
micronsPerPixelY = str2double(block.setup.XML.micronsPerPixelY);
if ~isequal(micronsPerPixelX, micronsPerPixelY)
    error('X and Y dimensions do not match')
end

%Suite2p ROI coordinates
cellXY = nan(nCells,2);
for c = 1:nCells
    cellXY(c,1) = double(stat{1,c}.med(2))*micronsPerPixelX;
    cellXY(c,2) = double(stat{1,c}.med(1))*micronsPerPixelX;
end

%Average Suite2p x/y offset
xoff = block.ops.xoff;
yoff = block.ops.yoff;
avgxoff = mean(xoff,'omitnan')*micronsPerPixelX;
avgyoff = mean(yoff,'omitnan')*micronsPerPixelX;

%% Set up StimTable

StimTable = table;
StimTable.EnsembleName   =  strings(nTrials,1);
StimTable.TargetedCells  =  cell(nTrials,1);
StimTable.Distance       =  cell(nTrials,1);
StimTable.PointIndex     =  cell(nTrials,1); %Which XY positions do the targeted cells correspond to
StimTable.nCells         =  nan(nTrials,1); %Can be different than nPoints if not all targets had cells < 15um from them
StimTable.TrialType      =  nan(nTrials,1);
StimTable.Stim_V1        =  nan(nTrials,1);
StimTable.Stim_V2        =  nan(nTrials,1);
StimTable.ToscaGroup     =  ToscaGroups;
StimTable.UncShutter     =  uncagingShutter;
StimTable.IgnoreTrial    =  block.parameters.trialsToIgnore';
StimTable.Duration       =  block.setup.XML.markPoints.Duration;
StimTable.Power          =  block.setup.XML.markPoints.Power;
StimTable.Revs           =  block.setup.XML.markPoints.SpiralRevolutions;
StimTable.PointNames     =  block.setup.XML.markPoints.PointNames;
StimTable.nPoints        =  nan(nTrials,1);
StimTable.Size           =  nan(nTrials,1);
StimTable.X              =  cell(nTrials,1);
StimTable.Y              =  cell(nTrials,1);

for i = 1:nTrials
    %Spiral size is the same for all points
    StimTable.Size(i) = block.setup.XML.markPoints.Points{i,1}.SpiralSizeInMicrons;
    
    %nPoints can vary per ensemble
    tempPoints = block.setup.XML.markPoints.Points(i,:);
    tempPoints = tempPoints(~cellfun('isempty',tempPoints)); %Remove empty elements
    StimTable.nPoints(i) = size(tempPoints,2);
    
    %Get X/Y location for each point
    [StimTable.X{i}, StimTable.Y{i}] = deal(nan(1,StimTable.nPoints(i)));
    for j = 1:StimTable.nPoints(i)
        StimTable.X{i}(j) = block.setup.XML.markPoints.Points{i,j}.X*512*micronsPerPixelX;
        StimTable.Y{i}(j) = block.setup.XML.markPoints.Points{i,j}.Y*512*micronsPerPixelX;
    end
end

%% Apply x/y offset correction
% We don't know whether to add/subtract the Suite2p xoff and yoff from each point
% Figure out which combo will minimize the distance from each point to the nearest cell
% This is to accommodate for any systematic shift in cell position due to Suite2p registration

Groups = unique(StimTable.PointNames);

if isempty(ops.UserOffset)
    
    Offset = [0, 0;...
            -avgxoff, -avgyoff;...
            -avgxoff, avgyoff;...
            avgxoff, -avgyoff;...
            avgxoff, avgyoff];

    OffsetDistance = [];
    for O = 1:height(Offset)
        allXY = [];
        for G = 1:length(Groups)
            GroupRows = strcmp(StimTable.PointNames,Groups{G});
            GroupTrials = find(GroupRows);
            pointXY = [StimTable.X{GroupTrials(1)}',StimTable.Y{GroupTrials(1)}'];
            allXY = [allXY; pointXY];
        end
        allXY(:,1) = allXY(:,1) + Offset(O,1); %Apply x offset
        allXY(:,2) = allXY(:,2) + Offset(O,2); %Apply y offset

        %Compute distance from each offset point to nearest cell
        distanceFromTarget = nan(nCells,size(allXY,1));
        for p = 1:size(allXY,1)
            for c = 1:nCells
                distanceFromTarget(c,p) = norm(cellXY(c,:) - allXY(p,:)); %Euclidean distance
            end
        end
        [closestDistance, ~] = min(distanceFromTarget,[],1);
        OffsetDistance = [OffsetDistance; closestDistance];
    end

    %Which offset minimizes point distance UNDER ops.MaxTargetDistance?
    OffsetDistance(OffsetDistance > ops.MaxTargetDistance) = nan;
    allNanColumns = sum(isnan(OffsetDistance),1) == height(OffsetDistance);
    OffsetDistance(:,allNanColumns) = [];
    if isempty(OffsetDistance)
        FinalOffset = [0 0];
    else
        [~,ind] = min(OffsetDistance,[],1, 'omitnan');
        modeInd = mode(ind);
        [~,meanInd] = min(mean(OffsetDistance,2));
        % if modeInd == meanInd
            OffsetInd = modeInd;
        % else
        %     error('Troubleshoot')
        % end
        FinalOffset = Offset(OffsetInd,:);
    end
else
    FinalOffset = ops.UserOffset;
end

%Apply offset to all points
for i = 1:height(StimTable)
    StimTable.X{i} = StimTable.X{i} + FinalOffset(1);
    StimTable.Y{i} = StimTable.Y{i} + FinalOffset(2);
end

MarkpointsExperiment.OffsetXY = FinalOffset;
disp(['Applied ' num2str(round(FinalOffset(1))) ' micron X offset and ' num2str(round(FinalOffset(2))) ' micron Y offset'])

%% Figure out which cell(s) were targeted

for G = 1:length(Groups)
    GroupRows = strcmp(StimTable.PointNames,Groups{G});
    GroupTrials = find(GroupRows);
    ToscaGroup = unique(StimTable.ToscaGroup(GroupRows));
    nPoints = max(StimTable.nPoints(GroupRows));
    pointXY = [StimTable.X{GroupTrials(1)}',StimTable.Y{GroupTrials(1)}'];
    distanceFromTarget = nan(nCells,nPoints);
    
    %Catch issues
    for t = 1:length(GroupTrials)
        temp_pointXY = [StimTable.X{GroupTrials(t)}',StimTable.Y{GroupTrials(t)}'];
        if ~isequal(pointXY, temp_pointXY)
            error('Target XY position varies: new situation')
        end
    end

    %Find distance between targeted area and each cell
    for p = 1:nPoints
        for c = 1:nCells
            distanceFromTarget(c,p) = norm(cellXY(c,:) - pointXY(p,:)); %Euclidean distance
        end
    end
    [temp_closestCellDistance, temp_closestCellIndex] = min(distanceFromTarget,[],1);
    pointIndex = 1:nPoints;
        
    %Check if any of the target points are outside of the cropped Suite2p
    %region (if so, then we won't be able to find the correct target cells there)
    if isfield(block.img, 'xrange') %Online BOTS don't have xrange
        %First find region of cropped rectangle
        cropx = double(block.img.xrange);
        cropy = double(block.img.yrange);
        cropWidth = length(cropx(1):cropx(2));
        cropHeight = length(cropy(1):cropy(2));
        cropSquareX = [cropx(1):cropx(2), zeros(1,cropHeight) + cropx(2), fliplr(cropx(1):cropx(2)), zeros(1,cropHeight) + cropx(1)];
        cropSquareY = [zeros(1,cropWidth) + cropy(1), cropy(1):cropy(2), zeros(1,cropWidth) + cropy(2), fliplr(cropy(1):cropy(2))];

        %Check if all points are inside the rectangle
        insideCroppedRegion = nan(1,nPoints);
        for p = 1:nPoints
            xInPixels = pointXY(p,1)./(micronsPerPixelX);
            yInPixels = pointXY(p,2)./(micronsPerPixelX);
            insideCroppedRegion(p) = inpolygon(xInPixels,yInPixels,cropSquareX,cropSquareY);
        end

        %If not, AND the closest found cell is > ops.MaxTargetDistance away, discard the closest cell as a potential target
        croppedPoints = ((temp_closestCellDistance > ops.MaxTargetDistance) + ~insideCroppedRegion) == 2;
        temp_closestCellDistance(croppedPoints) = [];
        temp_closestCellIndex(croppedPoints) = [];
        pointIndex(croppedPoints) = [];
    end

    %Discard any cells that are > ops.MaxTargetDistance away from the target
    removePoints = temp_closestCellDistance > ops.MaxTargetDistance;
    temp_closestCellDistance(removePoints) = [];
    temp_closestCellIndex(removePoints) = [];
    pointIndex(removePoints) = [];
    
    %Finally, make sure that each cell can only be selected once and sort numerically
    [closestCellIndex, sortind] = unique(temp_closestCellIndex); %Sorts in numerical order at the same time
    closestCellDistance = nan(size(closestCellIndex));
    for c = 1:length(closestCellIndex)
        closestCellDistance(c) = min(temp_closestCellDistance(temp_closestCellIndex == closestCellIndex(c)));
    end
    TargetedCells = cellNumbers(closestCellIndex); %Targeted cells are the closest cells to the target region
    pointIndex = pointIndex(sortind);
    nClosestCells = length(closestCellIndex);
    
    %% FIGURE: Plot traces and FOV to get user to confirm
    figure('units','normalized','outerposition',[0 0 1 1]); hold on
        
    %FOV
    ax(1) = subplot(3, nClosestCells + 4, [1:3; (1:3) + nClosestCells + 4; (1:3) + (nClosestCells + 4)*2]); hold on
    img = imadjust(int16(block.img.meanImg));
    imagesc(img);
    if isfield(block.img, 'xrange')
        plot(cropSquareX, cropSquareY, '--w') %Plot cropped region
    end
    xlim([0 512])
    ylim([0 512])
    colormap(ax(1), 'bone')
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'YDir','reverse')
    title([num2str(nPoints) ' targeted spots and detected target cells'])
    subtitle(strcat('Uncaging shutter: ', {' '}, MarkpointsExperiment.UncagingShutter))

    %Plot X over targeted regions
    for p = 1:nPoints
        xInPixels = pointXY(p,1)./(micronsPerPixelX);
        yInPixels = pointXY(p,2)./(micronsPerPixelX);
        scatter(xInPixels, yInPixels, 500, 'r', 'x', 'LineWidth', 4)
    end
    
    %Plot cell boundaries
    for p = 1:length(closestCellIndex)
        if isfield(stat{1,closestCellIndex(p)}, 'xpix') %Online BOTs don't have xpix
            xpoints = double(stat{1,closestCellIndex(p)}.xpix);
            ypoints = double(stat{1,closestCellIndex(p)}.ypix);
            roi = zeros(size(img));
            roi((xpoints-1)*512+ypoints) = 1; 
            bound = bwboundaries(roi);  
            x = []; y = [];
            for ib = length(bound)
               x = [x; bound{ib}(:, 2)];
               y = [y; bound{ib}(:, 1)];
            end
        else
            x = double(stat{1,closestCellIndex(p)}.xcirc);
            y = double(stat{1,closestCellIndex(p)}.ycirc);
        end
        
        if redcell(closestCellIndex(p))
            plot(x, y, 'Linewidth', 1.5, 'Color', 'r');
        else
            plot(x, y, 'Linewidth', 1.5);
        end
        text(max(x),max(y),num2str(cellNumbers(closestCellIndex(p))), 'Color', 'w');
    end
    
    %TRACES
    ylimits = []; h = {};
    if nClosestCells > 0
        subplot(3, nClosestCells + 4, 5:(4+nClosestCells)); hold on
        for p = 1:nClosestCells
            F7 = block.F7(closestCellIndex(p),:);
            mean_F7 = mean(F7, 'omitnan');
            DFF = (F7 - mean_F7)./mean_F7;
            plot(block.timestamp, DFF*.15 + p,'LineWidth',1);
        end
        timestamp_no_nan = block.timestamp(~isnan(block.timestamp));
        xlim([1 timestamp_no_nan(end)])
        set(gca,'YTick',1:nClosestCells)
        set(gca,'YTickLabel', num2str(cellNumbers(closestCellIndex)))
        set(gca,'YAxisLocation', 'right')
        %vline(block.Sound_Time(GroupTrials),'r')
        ylabel('Cell')
        xlabel('Seconds')
        title(strcat('Markpoints Group: ', {' '}, Groups{G}, {' '}, ' Tosca Group: ', {' '}, ToscaGroup))

        %AVERAGES AND RASTERS
        for p = 1:nClosestCells
            cellTrials = squeeze(Trials(closestCellIndex(p),GroupTrials,:));
            cellTrials(isinf(cellTrials)) = nan;

            h{p} = subplot(3, nClosestCells + 4, 8 + nClosestCells + p); hold on
            y = mean(cellTrials,1,'omitnan');
            if block.setup.framerate == 30
                y = smooth(y,5);
            end
            y_sem = std(cellTrials,1,'omitnan')./sqrt(length(y)-1);
            shadedErrorBar(xlin, y, y_sem)
            vline(baseline_length, 'r')
            xlim([xlin(1) xlin(end)])
            xlabel('Seconds')
            title(['Cell ' num2str(cellNumbers(closestCellIndex(p)))])
            ylimits = [ylimits; h{p}.YLim];
            if p == 1; ylabel('Z'); end

            ax(2) = subplot(3, nClosestCells + 4, 12 + 2*nClosestCells + p); hold on
            smoothCellTrials = cellTrials;
            for c = 1:size(cellTrials,1)
                smoothCellTrials(c,:) = smooth(cellTrials(c,:),5);
            end
            imagesc(cellTrials)
            xlim([0.5 size(cellTrials,2)])
            ylim([0.5 size(cellTrials,1)])
            caxis([-3 6])
            colormap(ax(2), bluewhitered(256))
            set(gca, 'XTick', xtick)
            set(gca, 'XTickLabel', xticklabel)
            xlabel('Seconds')
            if p == 1; ylabel('Trials'); end
        end
        
        %Set ylims on plots
        ymax = max(ylimits(:,2));
        ymin = min(ylimits(:,1));
        for p = 1:nClosestCells
            axes(h{p})
            try
            ylim([ymin ymax])
            catch
            end
            vline(nBaselineFrames, 'k')
        end
    else
        subplot(3, 1 + 4, 5:(4+1)); hold on
        title(strcat('Markpoints Group: ', {' '}, Groups{G}, {' '}, ' Tosca Group: ', {' '}, ToscaGroup))
        subtitle('No ROIs found within cropped region')
    end
    sgtitle(regexprep(block.setup.block_supname,'_',' ','emptymatch'))
        
    %% Ask user if this was a sham group (e.g. part of the FOV with no cells in it)
    % If not, get user to input Ensemble Name: try as hard as possible to make these standardized!
        
    clear 'ensembleName' %Make sure this variable does not exist prior to this section
    sham_group = 99; %Number that is not 0 or 1
    
    %Check if ensemble name is specified in info sheet
    if ~ismissing(MarkpointsExperiment.Ensemble)
        ensembleName = MarkpointsExperiment.Ensemble;
        if strcmp(ensembleName,'Sham')
            disp(strcat('Ensemble name found: ', {' '}, ensembleName))
            disp('Sham stimulation detected from Info sheet')
            sham_group = 1;
        elseif ~strcmp(ensembleName,'Sham') && (length(Groups) == 1)
            sham_group = 0;
        end
    end
        
    %Check if prevStimTable exists
    %This is so user doesn't have to type in ensemble name again
    if exist('prevStimTable', 'var')
        for p = 1:size(prevStimTable,1)
            if isequal(TargetedCells, prevStimTable.TargetedCells{p})
                ensembleName = prevStimTable.EnsembleName(p);
                disp(strcat('Ensemble name found: ', {' '}, ensembleName))
                sham_group = 0;
                break;
            elseif isequal(pointXY(:,1)', prevStimTable.X{p})
                %Try to identify sham ensemble here
                if strcmp(prevStimTable.EnsembleName(p), 'Sham')
                    ensembleName = "Sham";
                    sham_group = 1;
                    disp('Sham stimulation detected from block');
                    break;
                end
            end
        end
    end
    
    if sham_group == 1
        %Go to next section
    elseif sham_group == 0
        %Go to next section
    elseif strcmp(MarkpointsExperiment.UncagingShutter, "Closed")
        %If shutter was closed, automatically make Sham group
        sham_group = 1;
        disp('Sham stimulation detected due to closed shutter')
    elseif strcmp(MarkpointsExperiment.UncagingShutter, "Variable")
        %If shutter was variable, automatically make Not Sham group
        sham_group = 0;
    else
        %If shutter was open and ensemble not defined in info sheet, ask user
        while ~any(ismember(sham_group, [0 1])) %Only allow input of 0 or 1
            sham_group = input('Is this a sham stimulation? 1 for yes, 0 for no: ');
        end
    end

    %Update stim table
    if ~sham_group
        if ~ismissing(MarkpointsExperiment.Ensemble)
            ensembleName = MarkpointsExperiment.Ensemble;
            disp(strcat('Ensemble name found: ', {' '}, ensembleName))
        else
            if ~exist('ensembleName','var') %It will exist if it was found in prevStimtable
                ok = 0;
                while ok ~= 1
                    ensembleName = input('Please input an ensemble name: ', 's');
                    disp(ensembleName);
                    ok = input('Confirm name (1) or re-enter (0): ');
                end
            end
        end
        StimTable.EnsembleName(GroupTrials) = ensembleName;

        %Sort targeted cells numerically
        for t = 1:length(GroupTrials)
            StimTable.TargetedCells{GroupTrials(t)} = TargetedCells;
            StimTable.Distance{GroupTrials(t)} = closestCellDistance;
            StimTable.PointIndex{GroupTrials(t)} = pointIndex;
            StimTable.nCells(GroupTrials(t)) = nClosestCells;
        end
    else
        ensembleName = "Sham";
        StimTable.EnsembleName(GroupTrials) = "Sham";
    end
    
    %Retitle figure before closing/saving with ensemble name
    sgtitle(strcat(regexprep(block.setup.block_supname,'_',' ','emptymatch'), {' - Ensemble Name: '}, ensembleName))

end

%Save figures (DO NOT CLOSE FIGURES)
if strcmp(ops.SaveFigPath, 'default')
    %Make Markpoints directory and save figures in each block
    mkdir(figurepath)
    save_plots(figurepath, [char(block.setup.block_filename) planeTag], 0);
else
    %Save all figures in one folder specified by ops.SaveFigPath
    save_plots(ops.SaveFigPath, [char(block.setup.block_filename) planeTag], 0);
end

%% Determine trial type

var1 = block.parameters.variable1';
var2 = block.parameters.variable2';
[v1, v2, ~, ~, blank_ind, trialsToIgnore, StimInfo] = simple_prepare_variables(block); %standardize stim

if any(isnan(uncagingShutter))
    error('Uncaging shutter should not be NaN. Troubleshoot.')
end

%Update sham trials based on uncaging shutter
StimTable.EnsembleName(uncagingShutter == 0) = 'Sham';
StimTable.TargetedCells(uncagingShutter == 0) = {[]};
StimTable.Distance(uncagingShutter == 0) = {[]};
sham_trials = strcmp(StimTable.EnsembleName, 'Sham');

%Determine trial type
trialType = nan(nTrials,1);

if isequal([var1, var2], [0 0])
    % NO SOUND: ONLY SHAM OR ACTIVATION
    trialType(sham_trials) = 0;
    trialType(~sham_trials) = 1;
else
    StimTable.Stim_V1 = v1;
    StimTable.Stim_V2 = v2;
    stim = (~sham_trials + blank_ind) == 2; %Stim and no sound
    sham = (sham_trials + blank_ind) == 2; %Sham and no sound
    sound = (sham_trials + ~blank_ind) == 2; %Sham and sound
    stim_and_sound = (~sham_trials + ~blank_ind) == 2; %Stim and sound
    
    trialType(sham) = 0;
    trialType(stim) = 1;
    trialType(sound) = 2; %Sound trials are actually sound + sham
    trialType(stim_and_sound) = 3;
    trialType(trialsToIgnore == 1) = 9; %Uncaging shutter opening
end
    
StimTable.TrialType = trialType;

%Update TBD control type
if strcmp(MarkpointsExperiment.ControlType, "TBD")
    if ~any(sham_trials)
        MarkpointsExperiment.ControlType = "NotInThisBlock";
    else
        MarkpointsExperiment.ControlType = "ControlPointsWithShutterOpen";
    end
    %Display updated parameters
    disp(struct2table(MarkpointsExperiment));
end


%% Save to block and save .mat file in block folder

MarkpointsExperiment.StimInfo = StimInfo;
MarkpointsExperiment.StimTable = StimTable;
block.MarkpointsExperiment = MarkpointsExperiment;
save(filepath, 'MarkpointsExperiment')

end %end function