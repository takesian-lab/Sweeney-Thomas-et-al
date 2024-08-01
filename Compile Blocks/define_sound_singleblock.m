function [block] = define_sound_singleblock(block, WFtrigChange) 
% This function obtains the Bruker timestamp from the BOT and VR csv files,
% aligns the timestamp with the Tosca data (if available), and stores all of this in block
% 
% Argument(s):
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes: This function uses the Voltage recording file. Over time, we have
% added addition inputs into our Bruker system. As of 6/25/20, they should
% be the following
% 1) timestamp
% 2) State trigger
% 3) Rep trigger
% 4) trial trigger (we use this to align the timestamps of each trial
% 5) frame trigger from camera (for widefield control only)
% 6) speaker trigger
%   - older data sets will only have columns 1-4
%
% TODO:
% Search 'TODO'
%
% VERSIONS:
% v1 = archived May 2023
% v2 = current version, improves parsing of XML files

%% Skip this function if Bruker data is not available

if ismissing(block.setup.block_path) && ismissing(block.setup.VR_path)
    disp('Skipping Bruker data...');
    return
end

%Set WFtrigChange to 0 by default if not specified by user
if nargin < 2
    WFtrigChange = 0;
end

%Set substitute_ToscaTimes_for_VR to 0 by default if not defined in block
if ~isfield(block.setup, 'substitute_ToscaTimes_for_VR')
    block.setup.substitute_ToscaTimes_for_VR = 0;
end

disp('Pulling out Bruker data...');

%% Go to block and VR folders and pull out XML, BOT and voltage recording files
%The voltage recording is what we use to synchronize BOT data and Tosca data

%LOAD XML
cd(block.setup.block_path) %2p
disp('Reading XML file')
filedir = dir;
filenames = {filedir(:).name};
XML_files = filenames(contains(filenames,'.xml'));
X = readstruct(XML_files{1}, 'AttributeSuffix', '');
    
%WIDEFIELD
if ~ismissing(block.setup.VR_name)
    cd(block.setup.VR_path) %Widefield
    filedir = dir;
    filenames = {filedir(:).name};
    
    if WFtrigChange ==1
        %Load separate VR XML
        vrXML_files = filenames(contains(filenames,'.xml'));
        vrX = readstruct(vrXML_files{1}, 'AttributeSuffix', '');
        botCamTime = string(X.Sequence.time);
        vrCamTime = string(vrX.Sequence.time);
        timeRecorded = [vrCamTime; botCamTime];
        td = duration(timeRecorded,'inputformat','hh:mm:ss.SSSSSSS');
        out = seconds(diff(td));
    end
end

%Find voltage recording files
csvFiles = filenames(endsWith(filenames, '.csv')); %find all .csv files
VR_files = csvFiles(contains(csvFiles,'VoltageRecording')); %find all files that say VoltageRecording
VR_filename = VR_files(~contains(VR_files, 'botData.csv')); %remove any BOT files (which are also .csv)

%Check for markpoints
hasMarkPoints = false;
MarkPoints_files = filenames(contains(filenames,'MarkPoints.xml'));
if ~isempty(MarkPoints_files)
    block.setup.hasMarkPoints = true;
    hasMarkPoints = true;
end

%Record voltage recording name in setup
if isempty(VR_filename) && hasMarkPoints %markpoints recordings don't always have a voltage recording
    block.setup.VR_filename = nan;
elseif ~isempty(VR_filename)
    if length(VR_filename) == 1
        VR_filename = string(VR_filename{:}); %Convert to string if only one VR file
    end
    block.setup.VR_filename = VR_filename; %If more than one VR file, this will record all of the filenames
elseif ~block.setup.substitute_ToscaTimes_for_VR
    %If data was corrupted and voltage recording is missing, you can try to
    %recover it by using Tosca Times instead. Instructions TBD
    error('Voltage Recording missing')
end

%% Align Tosca data with Bruker times
%now lets get the relevant data from Voltage recording. This is when Tosca
%sends a 5V trigger to the Bruker system. The first time this occurs is t=0 on Tosca

%Lucas ephys data - get voltage recording triggers without aligning to Tosca
if block.setup.stim_protocol == 17
    %Load voltage recording
    disp(['Loading ' char(VR_filename)])
    M = csvread(VR_filename,1,0); %see notes at top for description of M
    
    % find the start of each trial (ephys or laser activation) using ACT TRIGGER, and align times to it
    Bruker_trial_triggers = M(M(:,4) >= 1)./1000;
    diffTrigger = diff([0; Bruker_trial_triggers]); %append 0 to the start so that we won't be off by 1 when using diff function
    Bruker_trial_time = Bruker_trial_triggers(diffTrigger > 1);
    
    %Save in block
    block.Sound_Time = Bruker_trial_time;
    block.voltage_recording = M;
    
    %Make block compatible with other scripts that look for stim parameters
    block.parameters.variable1 = 0; 
    block.parameters.variable2 = 0;
    block.parameters.stimLength = nan(1, length(Bruker_trial_time));
end

%Align VR triggers with Tosca data
if ~ismissing(block.setup.Tosca_path) && ~hasMarkPoints %Skip if Tosca info is missing or if MarkPoints data

    %Read VR csv to get Tosca trigger times
    if ~block.setup.substitute_ToscaTimes_for_VR
        disp(['Loading ' char(VR_filename)])
        M = csvread(VR_filename,1,0); %see notes at top for description of M
        if all(isinf(M(:,4))); error('Voltage recording file corrupted'); end %Catch problem where VR file is corrupted

        %Triggers are 5V signals that are active for 19-20ms
        Bruker_trial_triggers = M(M(:,4) > 3)./1000; %Detect all timepoints in the VR that are > 3V and convert to seconds
        diffTrigger = diff([0; Bruker_trial_triggers]); %append 0 to the start so that we won't be off by 1 when using diff function
        Bruker_trial_time = Bruker_trial_triggers(diffTrigger > 1);

        %Check for empty VR file
        if isempty(Bruker_trial_time)
            error('No triggers found in Voltage recording')
        end

    else
        %Option to salvage blocks without voltage recording file by temporarily replacing with Tosca Times
        %WARNING: This is not aligned with with Bruker yet! Must run script salvage_voltagerecording after block is compiled
        Bruker_trial_time = get_Sound_Time_behavior(block)';
    end
    
    %FYI: plot figure to visualize VR
    visualize_VR = 0;
    if visualize_VR
        visualize_VoltageRecording;
    end
        
    %WIDEFIELD ONLY
    if ~ismissing(block.setup.VR_path)
        if WFtrigChange ==1
            Bruker_trial_time = Bruker_trial_time - out;
        end
    end
    
    %Align Tosca variables (Loco, Licks) to Bruker timestamp
    %(This function is now where Bruker_trial_time becomes Sound_Time)
    [block] = align_Tosca_to_Bruker(block, Bruker_trial_time);
end

%% Find the timestamp for each frame using BOTs
% Let's find the time stamp for each frame. This requires to pull out
% the BOT data and correct for the short (>5ms delay) between Voltage recording and BOT

cd(block.setup.block_path)
filedir = dir;
filenames = {filedir(:).name};
BOT_files = filenames(endsWith(filenames, 'botData.csv'));

if length(BOT_files) == 1 %Regular BOT
    BOT_filename = BOT_files{1}; %Convert to string
    block.setup.BOT_filename = BOT_filename;
    disp(['Loading ' BOT_filename])
    try
        frame_data = csvread(BOT_filename, 1,0);
        timestamp = frame_data(:,1)-frame_data(1,1);% this is where we make that small correction
    catch
        %If the BOT cannot be read for some reason, use XML for timestamp
        block.setup.BOT_filename = nan;
    end
    % For z-corrected multiplane data we actually don't want to make this
    % correction so we will fix that below in the XML section
    singleBOT = true;

else %T-series or multiplane BOT

    block.setup.BOT_filename = nan;
    singleBOT = false;
    
    %This is how you would load all of the BOTs for a T-series markpoints or multiplane experiment
    %BUT it is actually faster (and maybe more accurate?) to use the XML
    %than to load all of these BOT files, so don't do it this way:
%     disp(['Loading ' num2str(length(BOT_files)) ' BOT files'])
%     frame_data = [];
%     for b = 1:length(BOT_files)
%         temp_frame_data = csvread(BOT_files{b}, 1,0);
%         frame_data = [frame_data; temp_frame_data];
%     end   
%     timestamp = frame_data(:,1)-frame_data(1,1);
end

%% Save XML data for widefield

if ~ismissing(block.setup.VR_name)
    
    PVStateTable = X.PVStateShard.PVStateValue;
    FrameTable = struct2table(X.Sequence.Frame);

    XML = struct;
    XML.filename            = string(XML_files{1});
    XML.date                = datestr(X.date);
    XML.systemID            = X.SystemIDs.SystemID_1.SystemID;
    XML.version             = X.version;
    XML.activeMode          = getXMLAttribute(PVStateTable, 'activeMode');
    XML.camera_binFactor    = getXMLAttribute(PVStateTable, 'camera_binFactor');
    XML.camera_exposureMode = getXMLAttribute(PVStateTable, 'camera_exposureMode');
    XML.camera_exposureTime = getXMLAttribute(PVStateTable, 'camera_exposureTime');
    XML.framePeriod         = getXMLAttribute(PVStateTable, 'framePeriod');
    XML.opticalZoom         = getXMLAttribute(PVStateTable, 'opticalZoom');
    
    %Determine pixel size in microns
    %Formula: Image pixel size = camera pixel size x binning / (obj. mag x lens mag)
    obj_mag = 4; %4x objective - we can't take this from the XML because sometimes the user does not select 4x from the drop-down menu
    camera_pixel_size = 6.5; %microns https://www.hamamatsu.com/us/en/product/cameras/cmos-cameras/C11440-42U30.html
    image_pixel_size = (camera_pixel_size*str2double(XML.camera_binFactor))/(obj_mag*str2double(XML.opticalZoom));
    XML.micronsPerPixelX = string(image_pixel_size); %Convert to string to keep consistent with other code down the line
    XML.micronsPerPixelY = string(image_pixel_size); 
    
    %Record timestamp
    XML.absoluteTime = FrameTable.absoluteTime;
    XML.relativeTime = FrameTable.relativeTime;

    framerate = ceil(1/str2double(XML.framePeriod));
    if framerate ~= 20
        warning(['Check frame rate. Detected rate is: ' num2str(1/str2double(XML.framePeriod))])
    else
        disp(['Framerate: ' num2str(framerate)]);
    end
end

%% Save XML data for 2P

if ismissing(block.setup.VR_name)
    
    PVStateTable = X.PVStateShard.PVStateValue;
    SequenceTable = struct2table(X.Sequence);
    
    XML = struct;
    XML.filename         = string(XML_files{1});
    XML.date             = datestr(X.date);
    XML.systemID         = X.SystemIDs.SystemID_1.SystemID;
    XML.version          = X.version;
    XML.activeMode       = getXMLAttribute(PVStateTable, 'activeMode');
    XML.framePeriod      = str2double(getXMLAttribute(PVStateTable, 'framePeriod'));
    XML.laserPower       = round(getXMLAttribute(PVStateTable, 'laserPower', 'Imaging'));
    XML.PMT1_Gain        = getXMLAttribute(PVStateTable, 'pmtGain', 'PMT 1 HV');
    XML.PMT2_Gain        = getXMLAttribute(PVStateTable, 'pmtGain', 'PMT 2 HV');
    XML.linesPerFrame    = str2double(getXMLAttribute(PVStateTable, 'linesPerFrame'));
    XML.micronsPerPixelX = string(getXMLAttribute(PVStateTable, 'micronsPerPixel', 'XAxis')); %Use string to round number
    XML.micronsPerPixelY = string(getXMLAttribute(PVStateTable, 'micronsPerPixel', 'YAxis')); %Use string to round number
    XML.opticalZoom      = getXMLAttribute(PVStateTable, 'opticalZoom');
    XML.x                = string(getXMLAttribute(PVStateTable, 'positionCurrent', [], 'XAxis'));
    XML.y                = string(getXMLAttribute(PVStateTable, 'positionCurrent', [], 'YAxis'));
    XML.z                = string(getXMLAttribute(PVStateTable, 'positionCurrent', [], 'ZAxis', 'Z Focus'));

    %Warn user
    if XML.linesPerFrame ~= 512
        warning('Number of pixels per frame is not 512. Check settings in PV before next imaging session.')
    end
        
    %Correct micronsPerPixelX/Y if wrong objective was selected
    objective = getXMLAttribute(PVStateTable, 'objectiveLens');
    if ~strcmp(objective, "Nikon 16x")
        warning('User did not have 16x objective selected, replacing micronsPerPixel with default value')
        if XML.linesPerFrame ~= 512
            error('Number of pixels per frame is not 512. Cannot use the same default value as usual')
        end
        XML.micronsPerPixelX = "1.6154";
        XML.micronsPerPixelY = "1.6154";
    end
    
    %% DETERMINE DATA TYPE
    sequenceTypes = unique(SequenceTable.type);
    
    if length(sequenceTypes) == 1 %Only one sequence type found
        if strcmp(sequenceTypes, "BrightnessOverTime") || strcmp(sequenceTypes, "TSeries Brightness Over Time Element")
            XML.DataType = "BOT"; %Single plane BOT
        elseif contains(sequenceTypes, "TSeries ZSeries")
            XML.DataType = "T-Series"; %Multiplane T-Series
        end
    else %More than one sequence type found, e.g. BOT and MarkPoints
        if any(strcmp(sequenceTypes, "TSeries Brightness Over Time Element"))
            XML.DataType = "T-Series"; 
        end
    end
    if ~isfield(XML, 'DataType'); error('New data type. Contact 2P slack for help'); end

    %% GET TIMESTAMP FROM XML
    
    if strcmp(XML.DataType, "BOT")
        XML.nPlanes = 1;
        FrameTable = struct2table(X.Sequence.Frame);
        XML.absoluteTime = FrameTable.absoluteTime;
        XML.relativeTime = FrameTable.relativeTime;
        framerate = 1/XML.framePeriod;
        
        %XML.relativeTime is virtually equivalent to timestamp:
        %isequal(round(XML.relativeTime, 1), round(timestamp,1)); %FYI
        
        %If BOT.csv could not be loaded for any reason, use XML
        if isnan(block.setup.BOT_filename)
            timestamp = XML.relativeTime;
        end
            
    elseif strcmp(XML.DataType, "T-Series") && ~hasMarkPoints
        
        try
            XML.nPlanes = size(SequenceTable.Frame{1},2);
        catch
            XML.nPlanes = size(SequenceTable.Frame(1,:),2);
        end
        XML.bidirectionalZ = strcmp(unique(SequenceTable.bidirectionalZ), "True");

        %Preallocate
        nFrames = 0;
        for s = 1:height(SequenceTable)
            try
                nFrames = nFrames + size(SequenceTable.Frame{s},2);
            catch
                nFrames = nFrames + size(SequenceTable.Frame(s,:),2);
            end
        end
        [XML.absoluteTime, XML.relativeTime, XML.cycle, XML.index, framePeriod] = deal(nan(nFrames,1));
        
        count = 1;
        for s = 1:height(SequenceTable)
            try
                FrameTable = struct2table(SequenceTable.Frame{s});
            catch
                FrameTable = struct2table(SequenceTable.Frame(s,:));
            end
            for f = 1:height(FrameTable)
                XML.cycle(count) = SequenceTable.cycle(s);
                XML.absoluteTime(count) = FrameTable.absoluteTime(f);
                XML.relativeTime(count) = FrameTable.relativeTime(f);
                XML.index(count) = FrameTable.index(f);
                framePeriod(count) = getXMLAttribute(FrameTable.PVStateShard(f).PVStateValue, 'framePeriod');
                count = count + 1;
            end
        end
        
        XML.multiplaneFramePeriod = mode(framePeriod);
        framerate = (1/XML.multiplaneFramePeriod)/XML.nPlanes;
        %BOT framerate does not account for the retrace frame that PV adds between cycles
        %XML.relativeTime and absoluteTime show that each cycle takes an extra frame to process
        %This is the same for both one-directional and bidirectional data
        %BOT_fs = unique(round(diff(timestamp),5)); %FYI
        %XML_fs = round(diff(XML.relativeTime),5); %FYI
        %XML_fs_betweenCycles = unique(XML_fs(XML.nPlanes:XML.nPlanes:end)); %FYI
        
        %Replace timestamp with XML.relativeTime ONLY if the BOT data hasn't already been corrected using zcorrect_multiplane (in which
        %case there will only be a single BOT in the folder)
        
        if singleBOT %Z-corrected multiplane data (now single plane)
            block.ZCorrected = true;
            timestamp = frame_data; %Use frame_data that hasn't been corrected since the times were originally generated from the XML
        else %Multiplane data or multiplane data that has been converted to single plane
            block.MultiplaneData = true;
            timestamp.combined = XML.relativeTime - XML.relativeTime(1); %Align timestamp to 0 and save
            for n = 1:XML.nPlanes
                plane = n - 1;
                planeName = ['plane' num2str(plane)];
                if  XML.bidirectionalZ
                    nFrames = length(XML.index);
                    nReps = ceil(nFrames/(XML.nPlanes*2));
                    index = repmat([1:XML.nPlanes, fliplr(1:XML.nPlanes)],1,nReps);
                    index = index(1:nFrames)';
                    planeInd = find(index == n);
                else
                    planeInd = find(XML.index == n);  
                end
                timestamp.(planeName) = timestamp.combined(planeInd);
            end
        end
        
    elseif strcmp(XML.DataType, "T-Series") && hasMarkPoints %T-Series MarkPoints
        
        %Determine if multiplane or singleplane data
        nFrames = 0; %should be same as length(timestamp)
        types = SequenceTable.type;
        planesPerSequence = nan(size(types));
        for s = 1:height(SequenceTable)
            if strcmp(types(s), "TSeries Brightness Over Time Element")
                planesPerSequence(s) = 1;
                nFrames = nFrames + size(SequenceTable.Frame{s},2);
            elseif strcmp(types(s), 'TSeries ZSeries Element')
                XML.bidirectionalZ = strcmp(unique(SequenceTable.bidirectionalZ), "True");
                try
                    planesPerSequence(s) = size(SequenceTable.Frame{s},2);
                catch
                    planesPerSequence(s) = size(SequenceTable.Frame(s,:),2);
                end
                nFrames = nFrames + planesPerSequence(s);
            end
        end
        XML.nPlanes = mode(planesPerSequence);
        if XML.nPlanes > 1; block.MultiplaneData = true; end
        
        %Preallocate
        XML.sequenceType = strings(height(SequenceTable),1); %define as either BOT or MarkPoints
        XML.parameterSet = strings(height(SequenceTable),1); %used to open uncaging shutter
        XML.savetime = SequenceTable.time; %time when the sequence was reached in T-series and saved as XML/BOT file
        [XML.sequence, XML.absoluteTime, XML.relativeTime, XML.index, framePeriod] = deal(nan(nFrames,1));
        XML.markPoints.filenames = strings(sum(strcmp(types, 'TSeries Mark Points')),1);
        XML.markPoints.sequence = nan(sum(strcmp(types, 'TSeries Mark Points')),1);
        
        %Loop through all sequences and frames
        count = 1;
        countMP = 1;
        for s = 1:height(SequenceTable)
            if strcmp(types(s), "TSeries Brightness Over Time Element") || strcmp(types(s), "TSeries ZSeries Element")
                XML.sequenceType(s) = 'BOT';
                try
                    FrameTable = struct2table(SequenceTable.Frame{s});
                catch
                    FrameTable = struct2table(SequenceTable.Frame(s,:));
                end
                XML.parameterSet(s) = unique(FrameTable.parameterSet);
                for f = 1:height(FrameTable)
                    XML.sequence(count) = s;
                    XML.absoluteTime(count) = FrameTable.absoluteTime(f);
                    XML.relativeTime(count) = FrameTable.relativeTime(f);
                    XML.index(count) = FrameTable.index(f);
                    framePeriod(count) = getXMLAttribute(FrameTable.PVStateShard(f).PVStateValue, 'framePeriod');
                    count = count + 1;
                end
                if isstruct(SequenceTable.MarkPoints{s})
                    XML.sequenceType(s) = 'Combined';
                    nMarkPoints = size(SequenceTable.MarkPoints{s},2);
                    for ss = 1:nMarkPoints
                        XML.markPoints.filenames(countMP) = SequenceTable.MarkPoints{s,ss}.filename;
                        XML.markPoints.sequence(countMP) = s;
                        countMP = countMP + 1;
                    end
                end
            elseif strcmp(types(s), 'TSeries Mark Points')
                XML.sequenceType(s) = 'MarkPoints';
                nMarkPoints = size(SequenceTable.MarkPoints,2);
                for ss = 1:nMarkPoints
                    XML.markPoints.filenames(countMP) = SequenceTable.MarkPoints{s,ss}.filename;
                    XML.markPoints.sequence(countMP) = s;
                    countMP = countMP + 1;
                end
            else
                error('New sequence type')
            end
        end
        
        %Determine final timestamp and framerate
        if ~isfield(block, 'MultiplaneData') %Single plane
            %XML.relativeTime is relative to the start of each sequence
            %XML.absoluteTime is virtually equivalent to timestamp:
            absoluteTime = XML.absoluteTime - XML.absoluteTime(1); %FYI
            %isequal(round(absoluteTime, 3), round(timestamp, 3)); %FYI
            timestamp = absoluteTime; %Use XML timestamp instead of BOT timestamp
            XML.markpointsFramePeriod = mode(framePeriod);
            framerate = (1/XML.markpointsFramePeriod)/XML.nPlanes;
            
        else %Multiplane
            %XML.multiplaneFramePeriod = mode(framePeriod);
            %framerate = (1/XML.multiplaneFramePeriod)/XML.nPlanes;
            %^I don't think this method is accurate for markpoints
            frameDiffs = diff(XML.relativeTime);
            frameDiffs = frameDiffs(abs(frameDiffs) < 1); %Remove wait time between sequences
            XML.multiplaneFramePeriod = mean(frameDiffs);
            framerate = (1/XML.multiplaneFramePeriod)/XML.nPlanes;

            %Replace timestamp with XML.absoluteTime (BOTs are not accurate for multiplane data)
            timestamp = struct;
            timestamp.combined = XML.absoluteTime - XML.absoluteTime(1);
            for n = 1:XML.nPlanes
                plane = n - 1;
                planeName = ['plane' num2str(plane)];
                if  XML.bidirectionalZ
                    nReps = ceil(nFrames/(XML.nPlanes*2));
                    index = repmat([1:XML.nPlanes, fliplr(1:XML.nPlanes)],1,nReps);
                    index = index(1:nFrames)';
                    planeInd = find(index == n);
                else
                    planeInd = find(XML.index == n);  
                end
                timestamp.(planeName) = timestamp.combined(planeInd);
            end
        end
    end

    %% DISPLAY FRAMERATE
    %For single plane data, make framerate conform to 30 
    % - The actual framerate varies slightly from 30, here we are rounding to 30 for simplicity
    % - If at any point we need to have the exact framerate we can remove this
    if round(framerate) == 30 || ceil(framerate) == 30 || floor(framerate) == 30
        framerate = 30;
    elseif round(framerate) == 15 || ceil(framerate) == 15 
        framerate = 15; %(*Some data was collected with framerate = 15 in the past)
    elseif ~isfield(block, 'MultiplaneData') && ~isfield(block, 'ZCorrected')
        warning(['Check frame rate. Detected rate is: ' num2str(1/XML.framePeriod)]);
    end
    disp(['Framerate: ' num2str(framerate)]);
    
    %% MARKPOINTS
    
    if hasMarkPoints
        
        %Read environment file to get info about uncaging shutter
        %UncagingShutter = true/1 means it was open during block
        %UncagingShutter = false/0 means it was closed during block
        try
            ENV = readstruct(strcat(XML_files{1}(1:end-4),'.env'), 'Filetype', 'xml', 'AttributeSuffix', '');
            utilityButtonInfo = getXMLAttribute(ENV.PVStateShard.PVStateValue, 'utilityButton', 'Unc Shutter');
            XML.uncagingShutter = strcmp(utilityButtonInfo, "True");
        catch
            %Some older .env files are not read by MATLAB ("invalid XML file")
            %If you still get an error with readtable and you know what
            %state the uncaging shutter was in, you can assign
            %XML.uncagingShutter = true or false and skip this try/catch loop
            ENV = readtable(strcat(XML_files{1}(1:end-4),'.env'), "Filetype", "text", "VariableNamingRule", "preserve");
            rowind = find(strcmp(ENV.(ENV.Properties.VariableNames{1}), '<PVStateValue key="utilityButton">'), 1, 'last')+1;
            if ~contains(ENV.(ENV.Properties.VariableNames{1}){rowind}, "True") && ~contains(ENV.(ENV.Properties.VariableNames{1}){rowind}, "False")
                error('Troubleshoot: This may not have grabbed the right line')
            end
            XML.uncagingShutter = contains(ENV.(ENV.Properties.VariableNames{1}){rowind}, "True");
        end
        
        %If markpoints was recorded within a BOT, initialize MarkPoints struct
        if strcmp(XML.DataType, "BOT") 
            %Save info about sequence to be consistent with T-Series MarkPoints
            XML.sequenceType{1,1} = 'Combined';
            XML.savetime{1,1} = X.Sequence.time;
            XML.markPoints = struct;
            XML.markPoints.filenames{1,1} = X.Sequence.MarkPoints.filename;
            XML.sequence = ones(height(FrameTable),1);
        end
        
        %Loop through markpoints files to save information about each activation
        if isequal(size(XML.markPoints.filenames), [1 1]) %One Markpoints file with row elements triggered by Tosca
            
            XX = readstruct(XML.markPoints.filenames{1,1}, 'AttributeSuffix', '');
            PVMarkPointTable = struct2table(XX.PVMarkPointElement);
            XML.markPoints.Iterations = XX.Iterations;
            XML.markPoints.IterationDelay = XX.IterationDelay;
            
            %Find dummy rows and remove
            DummyRows = (strcmp(PVMarkPointTable.UncagingLaser, 'None') + strcmp(PVMarkPointTable.TriggerFrequency, 'None')) > 0;
            XML.markPoints.startsWithDummyRow = 0;
            if ~isempty(DummyRows)
                if DummyRows(1) == 1
                    XML.markPoints.startsWithDummyRow = 1;
                end
                PVMarkPointTable(DummyRows == 1,:) = [];
            end
            
            %If the first element in the table was NOT a dummy row, then we need to rotate the order of all the iterations by 1
            %This is a 'peculiarity' of PV where the first row of Markpoints gets triggered by the BOT starting, meaning that Tosca triggers start on row 2        
            %TODO: Does this happen in T-Series multiplane as well???
            nElements = height(PVMarkPointTable);
            elementOrder = 1:nElements;
            if ~XML.markPoints.startsWithDummyRow
                elementOrder = [elementOrder(2:end), elementOrder(1)];
            end
            
            %Include enough Iterations to account for the number of Tosca trials
            %The number of iterations should be set arbitrarily higher than
            %the number of Tosca trials to make sure we don't run out.
            %If the trial ended up being an error in Tosca, Markpoints was still activated 
            nToscaTrials = length(block.errorData.result); %Include error trials
            errorTrials = block.errorData.error_trials;
            nToscaTrialsWithoutErrors = nToscaTrials - length(errorTrials);
            if XML.markPoints.Iterations*nElements < nToscaTrials
                error('User did not input enough MarkPoints iterations for the number of Tosca Trials. Add code to account for this situation')
            end
            XML.markPoints.Iterations = ceil(nToscaTrials/nElements);
            
            %Loop through all markpoints elements and repeat n Iterations
            count = 1;
            toscaTrial = 1;
            for i = 1:XML.markPoints.Iterations
                for mm = 1:length(elementOrder)
                    m = elementOrder(mm);
                    
                    %Stop when we have reached the number of Tosca Trials including errors
                    if count > nToscaTrialsWithoutErrors
                        continue;
                    end
                    
                    %Skip error trials
                    if any(ismember(errorTrials, toscaTrial)) 
                        toscaTrial = toscaTrial + 1;
                        continue;
                    end
                    
                    try
                        PVGalvoPointElement = struct2table(PVMarkPointTable.PVGalvoPointElement{m});
                    catch
                        PVGalvoPointElement = struct2table(PVMarkPointTable.PVGalvoPointElement(m));
                    end
                    XML.markPoints.Iteration(count,1)          = i;
                    XML.markPoints.Repetitions(count,1)        = PVMarkPointTable.Repetitions(m);
                    XML.markPoints.Power(count,1)              = PVMarkPointTable.UncagingLaserPower(m);
                    XML.markPoints.TriggerCount(count,1)       = PVMarkPointTable.Repetitions(m);
                    XML.markPoints.TriggerFrequency(count,1)   = PVMarkPointTable.TriggerCount(m);
                    XML.markPoints.TriggerSelection(count,1)   = PVMarkPointTable.TriggerFrequency(m);
                    XML.markPoints.Duration(count,1)           = PVGalvoPointElement.Duration;
                    XML.markPoints.SpiralRevolutions(count,1)  = PVGalvoPointElement.SpiralRevolutions;
                    XML.markPoints.InitialDelay(count,1)       = PVGalvoPointElement.InitialDelay;
                    XML.markPoints.InterPointDelay(count,1)    = PVGalvoPointElement.InterPointDelay;
                    XML.markPoints.AllPointsAtOnce(count,1)    = PVGalvoPointElement.AllPointsAtOnce;
                    XML.markPoints.PointNames(count,1)         = PVGalvoPointElement.Points;
                    XML.markPoints.Indices(count,1)            = PVGalvoPointElement.Indices;

                    %Loop through all points
                    PointStruct = PVGalvoPointElement.Point;
                    for p = 1:size(PointStruct,2)
                        XML.markPoints.Points{count,p} = PointStruct(p);
                    end
                    count = count + 1;
                    toscaTrial = toscaTrial + 1;
                end
            end
            
        else %multiple markpoints files, use each one as a single activation trigger
            
            %Loop through all markpoints.xml files
            XML.markPoints.Points = cell(1,1);
            for m = 1:length(XML.markPoints.filenames)            
                XX = readstruct(XML.markPoints.filenames{m}, 'AttributeSuffix', '');
                PVMarkPointTable = struct2table(XX.PVMarkPointElement);
                
                XML.markPoints.IterationDelay(m,1)            = XX.IterationDelay;
                XML.markPoints.Iteration(m,1)                 = XX.Iterations;
                XML.markPoints.Repetitions(m,1)               = PVMarkPointTable.Repetitions;
                XML.markPoints.Power(m,1)                     = PVMarkPointTable.UncagingLaserPower;
                XML.markPoints.Duration(m,1)                  = PVMarkPointTable.PVGalvoPointElement.Duration;
                XML.markPoints.SpiralRevolutions(m,1)         = PVMarkPointTable.PVGalvoPointElement.SpiralRevolutions;
                XML.markPoints.InitialDelay(m,1)              = PVMarkPointTable.PVGalvoPointElement.InitialDelay;
                XML.markPoints.InterPointDelay(m,1)           = PVMarkPointTable.PVGalvoPointElement.InterPointDelay;
                XML.markPoints.AllPointsAtOnce(m,1)           = PVMarkPointTable.PVGalvoPointElement.AllPointsAtOnce;
                XML.markPoints.PointNames(m,1)                = PVMarkPointTable.PVGalvoPointElement.Points;
                XML.markPoints.Indices(m,1)                   = PVMarkPointTable.PVGalvoPointElement.Indices;

                %Loop through all points
                nPoints = size(PVMarkPointTable.PVGalvoPointElement.Point,2);
                for p = 1:nPoints
                    XML.markPoints.Points{m,p} = PVMarkPointTable.PVGalvoPointElement.Point(p);
                end
            end
        end
        
        %Estimate markpoints timestamps based on XML
        %If markpoints was triggered externally, this will be corrected in define_markpoints_timestamp
        nMarkpoints = size(XML.markPoints.Points,1);
        time = zeros(nMarkpoints,1); %(ms)
        for t = 1:nMarkpoints
            %For now, this code is only written for 1 point or 1 group of points targeted with SLM 
            %Here we would need to accommodate InterPointDelay
            if length(XML.markPoints.Points{t,p}) > 1 && ~strcmp(XML.markPoints.AllPointsAtOnce{t},'True')
                error('Code does not accommodate multiple points yet')
            end

            %Determine timing of each activation
            if t == 1
                time(t) = XML.markPoints.InitialDelay(t);
                currentIteration = XML.markPoints.Iteration(t);
            else
                %Find transition between iterations
                if XML.markPoints.Iteration(t) ~= currentIteration
                    IterationDelay = XML.markPoints.IterationDelay;
                    currentIteration = XML.markPoints.Iteration(t);
                else
                    IterationDelay = 0;
                end

                time(t) = time(t-1) + XML.markPoints.Duration(t-1) + IterationDelay + XML.markPoints.InitialDelay(t);
            end
        end
        
        %Convert to seconds and align with start of timestamp
        if XML.nPlanes == 1
            XML.markPoints.time = time./1000 + timestamp(1); 
        else
            XML.markPoints.time = time./1000 + timestamp.combined(1); 
        end
    end
end

%% Confirm that Tosca data and Bruker data were recorded on the same day
%  This is to catch data entry mistakes in the Info spreadsheet

if ~ismissing(block.setup.Tosca_path) %Skip if Tosca info is missing
    tosca_day = datestr(block.setup.Tosca_date,2); %Get day from Tosca date
    tosca_time = datestr(block.setup.Tosca_date,13); %Get time from Tosca date
    bruker_day = datestr(XML.date,2); %Get day from Bruker date
    bruker_time = datestr(XML.date,13); %Get time from Bruker date
 
    %Check day
    if ~isequal(tosca_day, bruker_day)
        error(['Tosca (' tosca_day ') and Bruker (' bruker_day ') data come from different days. Check Info spreadsheet.'])
    end
    
    %Check time to within 5 minutes (Tosca computer has a different clock than Bruker computer,
    %so it is possible that time could be very off even if blocks are correct)
    if abs(etime(datevec(datenum(bruker_time)), datevec(datenum(tosca_time)))) > 5*60
        warning(['Possible issue with Tosca (' tosca_time ') and Bruker (' bruker_time ') start times. Check run numbers in Info spreadsheet.'])
    end
end

%% Record XML data, block data, and filenames

block.setup.XML = XML;
block.setup.framerate = framerate;
block.timestamp = timestamp;

end
