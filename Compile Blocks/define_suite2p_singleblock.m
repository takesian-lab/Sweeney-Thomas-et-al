function [block] = define_suite2p_singleblock(block,varargin)
% This function accesses the Fall.mat suite2p data and stores it in block
%
% Argument(s):
%   block (struct)
% Varargin:
%   UpsampleTo30 (Optional 0 or 1 to upsample data < 30Hz)
%   CheckOps (Optional struct) = user-specified ops to check against Fall.ops
%   DisplayTable (Optional 0 or 1 to display table of blocks & frame numbers from the function get_frames_from_Fall)
%
% Returns:
%   block (struct)
%
% Notes:
%
% Uses the function:
% -get_frames_from_Fall
%
% TODO:
% Search 'TODO'

%% Setup

% ------- parse varargin
p = inputParser; 

%USAGE: addOptional(p,'parametername',defaultvalue);

%Default is to display block information in table format in get_frames_from_Fall
addOptional(p, 'DisplayTable', 1); % Option (0 or 1) to print table of blocks and frame numbers to command line

%Option to upsample data that is < 30Hz
%User can overwrite uption by providing new value, otherwise take value stored in block
if ~isfield(block.setup, 'UpsampleTo30')
    %Old blocks might not have this field
    block.setup.UpsampleTo30 = 0;
end
addOptional(p, 'UpsampleTo30', block.setup.UpsampleTo30);

%Option to check fall.ops against user-specified ops.mat file
user_ops.checkOps = 0;
addOptional(p, 'UserOps', user_ops);

parse(p, varargin{:});

functionOps = p.Results; 
% ------- end parse varargin

checkOps = user_ops.checkOps;

%% Skip this function if Suite2p data is not available

if ismissing(block.setup.suite2p_path) || isequal(block.setup.suite2p_path, block.setup.block_path)
    disp('Skipping Suite2p data...');
    return
end

disp('Pulling out Suite2p data...');

%% Go to Suite2p folder and load plane0/Fall.mat
%  For multi-plane recordings, there will be one 'plane#' folder per plane

cd(block.setup.suite2p_path)

[~,currentFolder,~] = fileparts(pwd);
folderRoot = currentFolder(1:5); %Get 'plane' or 'suite' only

if strcmp(folderRoot, 'plane')
    %If the user enters a filepath ending in 'plane#', treat as single plane data
    nPlanes = 1;
    Fall = load('Fall.mat');  
else
    %If filepath is not a suite2p or plane# folder, set to \suite2p to attempt to find the folder
    if ~strcmp(folderRoot, 'suite')
        newpath = strcat(block.setup.suite2p_path, '\suite2p');
        if ~isfolder(newpath)
            error('Could not find suite2p folder. Check suite2p path')
        end
        block.setup.suite2p_path = string(newpath);
        cd(newpath)
        [~,currentFolder,~] = fileparts(pwd);
    end
    
    %If the user enters a filepath ending in 'suite2p', check how many planes there are
    currentDirectory = dir;
    names = {};
    for n = 1:size(currentDirectory,1)
        names{n,1} = currentDirectory(n).name;
    end
    planeFolders = names(contains(names, 'plane'));
    nPlanes = length(planeFolders);

    %Regardless of nPlanes, cd to plane0 and use as the default folder for getting ops info
    cd(strcat(block.setup.suite2p_path, '\plane0'))
    Fall = load('Fall.mat');
end

%% Determine if data is singleplane or multiplane and find Frame_set
%If the user enters a filepath ending in 'plane0', treat the data as if there is only one plane
%For multiplane data, this will allow for planes to be processed separately by specifying plane
%in the analysis_name column of info sheet
    
if isfield(block, 'MultiplaneData')
    %Multiplane timestamp found

    if nPlanes == 1
        %Treat as singleplane data
        disp('Processing multiplane data as single plane')
        
        if strcmp(currentFolder, 'suite2p')
            %Data was converted to single plane before compiling (because there is only one plane folder)
            %To find the original plane, use block name. Converted blocks should have the prefix 'Plane#-'
            split_block_name = split(block.setup.block_name, '-'); %split string to obtain prefix
            currentFolder = lower(split_block_name{1}); %make prefix all lowercase
            [Frame_set,~] = get_frames_from_Fall(Fall.ops, block.setup.block_name, functionOps.DisplayTable); %Pull out Frame_set for single plane
        else
            planeNumber = str2double(currentFolder(6:end));
            [~,Frame_set] = get_frames_from_Fall(Fall.ops, block.setup.block_name, functionOps.DisplayTable, planeNumber); %Use Plane_set for Frame_set
        end
            
        block.timestamp = block.timestamp.(currentFolder); %Replace with timestamp for current plane
        block = rmfield(block,'MultiplaneData'); %Remove multiplane field from block
        
        %Correct timestamp for Suite2p extra plane problem
        %(If the last cycle has nFrames < nPlanes, suite2p does not assign the extra frames correctly)
        %Plane_set only includes correct planes, so adjust timestamp to match the length of Plane_set
        %https://docs.google.com/presentation/d/1D32WRyEAAMaBo9NYC7Q9kB_p9UuGU262AKHHhZLscQs/edit?usp=sharing
        block.timestamp = block.timestamp(1:length(Frame_set));
    else
        %Keep multiplane data
        disp(['Found ' num2str(nPlanes) ' suite2p planes'])
        [Frame_set,~] = get_frames_from_Fall(Fall.ops, block.setup.block_name, functionOps.DisplayTable, 0); %Pull out Frame_set for multiplane data
    end
    
else
    %Single plane timestamp found
    [Frame_set,~] = get_frames_from_Fall(Fall.ops, block.setup.block_name, functionOps.DisplayTable); %Pull out Frame_set for single plane
end

%% Check for dropped frames when using multipage tiffs
%It seems like Suite2p sometimes drops frames for multipage tiff data
%We will try to detect and correct that (conservatively) here

%Determine if dataset has multipage tiffs
multipageTIFFs = false;
if isfield(Fall.ops, 'frames_per_file') %Old suite2p data do not have this variable (but therefore are also not going to be multipage)
    if any(Fall.ops.frames_per_file > 1)
        multipageTIFFs = true;
    end
end

%If it does, check for errors and correct
if multipageTIFFs
    if length(Frame_set) < length(block.timestamp)
        nFramesDiff = length(block.timestamp) - length(Frame_set);
        if nFramesDiff > 2
            error('Large difference between Suite2p frame number and timestamp. Consider whether to keep data.')
        end
        block.timestamp = block.timestamp(1:end-nFramesDiff); %Remove last frame of timestamp
        warning(['Suite2p error for multipage tiffs detected and corrected. Suite2p was ' num2str(nFramesDiff), ' frame(s) off.'])
    end
end

%% Check that Frame_set matches timestamp from Bruker function
% There should be an identical number of frames in the Bruker data as the Suite2p data

if ismissing(block.setup.block_path) && ismissing(block.setup.VR_path)
    %If Bruker data is not available (not included in info sheet)
    warning('Frame_set could not be checked against timestamp')
elseif isfield(block, 'MultiplaneData') && nPlanes > 1
    if length(Frame_set) ~= length(block.timestamp.combined)
        error('Frame_set does not match timestamp')
    end
else
    if length(Frame_set) ~= length(block.timestamp)
        %NaN-correction for some markpoints data
        if isfield(block, 'timestamp_idx')
            temp_timestamp = block.timestamp(block.timestamp_idx);
            
            if length(Frame_set)~= length(temp_timestamp)
                error('Frame_set does not match timestamp')
            end
        else
            error('Frame_set does not match timestamp')
        end
    end
end

%% Check ops

%Save ops in block
block.ops = get_abridged_ops(Fall.ops);

%Check that framerate matches block.ops.fs
if round(block.ops.fs) ~= round(block.setup.framerate)
    warning(['Suite2p framerate (' num2str(block.ops.fs) ') does not match block framerate (' num2str(block.setup.framerate) ')'])
end

%Check that neucoeff matches block.ops.neucoeff
if ~isequal(block.setup.constant.neucoeff, block.ops.neucoeff)
    warning(['User neucoeff (' num2str(block.setup.constant.neucoeff) ') does not match Suite2p neucoeff (' num2str(block.ops.neucoeff) ')'])
end

%Optional check ops against user-specified values
if checkOps
    [unmatchingOps, userOps, blockOps] = deal([]);
    fields = fieldnames(user_ops);
    for f = 1:numel(fields)
        currentField = fields{f};
        if strcmp(currentField, 'checkOps')
            continue
        elseif isfield(block.ops, currentField)
            if user_ops.(currentField) ~= block.ops.(currentField)
                unmatchingOps = [unmatchingOps; string(currentField)];
                userOps = [userOps; user_ops.(currentField)];
                blockOps = [blockOps; block.ops.(currentField)];
            end
        end
    end
    
    if ~isempty(unmatchingOps)
        warning('Some ops do not match the user file.')
        format long
        disp(table(unmatchingOps, userOps, blockOps))
    else
        disp('Ops match the user file.')
    end
end


%% Pull out data from Fall: SINGLE-PLANE DATA
% Fall is too big to keep in its entirety (a couple GB), so just keep the data we'll need

if nPlanes == 1
    
    %Images for visualization
    block.img.meanImg = Fall.ops.meanImg;
    block.img.refImg = Fall.ops.refImg;
    block.img.max_proj = Fall.ops.max_proj;
    block.img.meanImgE = Fall.ops.meanImgE;
    block.img.Vcorr = Fall.ops.Vcorr;
    block.img.xrange = Fall.ops.xrange; %Refers to cropped x/y of vcorr and max_proj
    block.img.yrange = Fall.ops.yrange;
    
    %Cell and neuropil data
    %Only keep data from 'is cells' within the block's frame set
    keep_ind = find(Fall.iscell(:,1) == 1);
    if isempty(keep_ind)
        error('No ROIs found in Suite2p data.')
    end
    block.cell_number = keep_ind-1; %Subtract 1 for the python to matlab correction of cell label
    block.stat = Fall.stat(1,keep_ind);
    block.ops.badframes = Fall.ops.badframes(Frame_set);
    block.ops.xoff = Fall.ops.xoff(:,Frame_set);
    block.ops.yoff = Fall.ops.yoff(:,Frame_set);
    block.F = Fall.F(keep_ind,Frame_set);
    block.Fneu = Fall.Fneu(keep_ind,Frame_set);
    block.spks = Fall.spks(keep_ind,Frame_set);
    block.F7 = block.F - block.setup.constant.neucoeff*block.Fneu; %Neuropil corrected trace
    block.df_f = (block.F7 - mean(block.F7, 2, 'omitnan'))./mean(block.F7, 2, 'omitnan'); %DF/F = (total-mean)/mean
    
    %Catch manual labeling problem
    %If allow_overlap is set to 0 in your ops and you attempt manual labeling,
    %if the new ROI is fully overlapping with another ROI then the F trace will be all zeros
    if any(sum(block.F,2) == 0)
        error('Detected Manual Labeling problem. Change Ops.npy and re-exract ROIs')
    end
    
    %Save channel 2 data
    if isfield(Fall, 'F_chan2') && ~isempty(Fall.F_chan2)
        block.F_chan2 = Fall.F_chan2(keep_ind,Frame_set);
        block.Fneu_chan2 = Fall.Fneu_chan2(keep_ind,Frame_set);
        block.F7_chan2 = block.F_chan2 - block.setup.constant.neucoeff*block.Fneu_chan2; %Neuropil corrected trace
        block.df_f_chan2 = (block.F7_chan2 - mean(block.F7_chan2))./mean(block.F7_chan2); %DF/F = (total-mean)/mean
    end
    
    %Save red channel data
    try
        redcell = Fall.redcell;
    catch
        error('At least one of the expected Fall.mat fields/variables not found; go back to suite2p and save to .mat (File > save to .mat)');
    end
    
    %Look for manual redcells.mat variable
    %This variable is designed to let you define your own red cells (e.g. if you only recorded one channel but you'd like to label red cells)
    %redcells.mat must be a N x 1 vector of suite2p ROI labels where N is the number of red cells
    if isfile('redcells.mat')
        disp('Found manual red cell labels')
        load('redcells.mat', 'redcells');
        manual_redcells = redcells;
        redcell(:) = 0;
        redcell(manual_redcells + 1) = 1; %Add 1 for python to MATLAB transformation
    end

    if ~isempty(redcell)% Added by VA (04/05/22)
        block.redcell = redcell(keep_ind);
        if any(block.redcell)
            disp('Found red cells')
        end
    end    
       
    %Update zcorr frame set
    if isfield (Fall.ops, 'zcorr')
        disp('Found zcorr')
        block.zcorr = Fall.ops.zcorr(:,Frame_set); %Dimensions are z-stack position vs. frame
    end
    
    %Determine distances between neurons
    block.distance_in_microns = find_cell_distance(block);
end
    
%% Pull out data from Fall: MULTIPLANE DATA
    
if nPlanes > 1
    for n = 0:nPlanes-1
        currentPlane = strcat('plane', num2str(n));
        disp(['Processing ' currentPlane])
        
        %CD to current plane folder and load Fall.mat to get Plane_set
        cd(strcat(block.setup.suite2p_path, '/', currentPlane))
        Fall = load('Fall.mat');
        [~,Plane_set] = get_frames_from_Fall(Fall.ops, block.setup.block_name, 0, n);

        %Images for visualization
        block.img.(currentPlane).meanImg = Fall.ops.meanImg;
        block.img.(currentPlane).refImg = Fall.ops.refImg;
        block.img.(currentPlane).max_proj = Fall.ops.max_proj;
        block.img.(currentPlane).meanImgE = Fall.ops.meanImgE;
        block.img.(currentPlane).Vcorr = Fall.ops.Vcorr;
        block.img.(currentPlane).xrange = Fall.ops.xrange; %Refers to cropped x/y of vcorr and max_proj
        block.img.(currentPlane).yrange = Fall.ops.yrange;
    
        %Prepare to save xoff and yoff per plane
        if n == 0
            block.ops.badframes = struct;
            block.ops.xoff = struct;
            block.ops.yoff = struct;
        end
        
        %Cell and neuropil data
        %Only keep data from 'is cells'
        keep_ind = find(Fall.iscell(:,1) == 1);
        if isempty(keep_ind)
            warning('No ROIs found in Suite2p data.')
        end
        block.cell_number.(currentPlane) = keep_ind-1; %Subtract 1 for the python to matlab correction of cell label
        block.stat.(currentPlane) = Fall.stat(1,keep_ind);
        block.ops.badframes.(currentPlane) = Fall.ops.badframes(Plane_set);
        block.ops.xoff.(currentPlane) = Fall.ops.xoff(:,Plane_set);
        block.ops.yoff.(currentPlane) = Fall.ops.yoff(:,Plane_set);
        block.F.(currentPlane) = Fall.F(keep_ind,Plane_set);
        block.Fneu.(currentPlane) = Fall.Fneu(keep_ind,Plane_set);
        block.spks.(currentPlane) = Fall.spks(keep_ind,Plane_set);
        block.F7.(currentPlane) = block.F.(currentPlane) - block.setup.constant.neucoeff*block.Fneu.(currentPlane); %Neuropil corrected trace
        block.df_f.(currentPlane) = (block.F7.(currentPlane) - mean(block.F7.(currentPlane), 2, 'omitnan'))./mean(block.F7.(currentPlane), 2, 'omitnan'); %DF/F = (total-mean)/mean

        %Correct timestamp for Suite2p extra plane problem
        %(If the last cycle has nFrames < nPlanes, suite2p does not assign the extra frames correctly)
        %Plane_set only includes correct planes, so adjust timestamp to match the length of Plane_set
        %https://docs.google.com/presentation/d/1D32WRyEAAMaBo9NYC7Q9kB_p9UuGU262AKHHhZLscQs/edit?usp=sharing
        block.timestamp.(currentPlane) = block.timestamp.(currentPlane)(1:length(Plane_set));

        %Catch manual labeling problem
        %If allow_overlap is set to 0 in your ops and you attempt manual labeling,
        %if the new ROI is fully overlapping with another ROI then the F trace will be all zeros
        if any(sum(block.F.(currentPlane),2) == 0)
            error('Detected Manual Labeling problem. Change Ops.npy and re-exract ROIs')
        end

        %Save channel 2 data
        if isfield(Fall, 'F_chan2')
            block.F_chan2.(currentPlane) = Fall.F_chan2(keep_ind,Plane_set);
            block.Fneu_chan2.(currentPlane) = Fall.Fneu_chan2(keep_ind,Plane_set);
            block.F7_chan2.(currentPlane) = block.F_chan2.(currentPlane) - block.setup.constant.neucoeff*block.Fneu_chan2.(currentPlane); %Neuropil corrected trace
            block.df_f_chan2.(currentPlane) = (block.F7_chan2.(currentPlane) - mean(block.F7_chan2.(currentPlane)))./mean(block.F7_chan2.(currentPlane)); %DF/F = (total-mean)/mean
        end

        %Save red channel data
        try
            redcell = Fall.redcell;
        catch
            error('At least one of the expected Fall.mat fields/variables not found; go back to suite2p and save to .mat (File > save to .mat)');
        end

        %Look for manual redcells.mat variable
        if isfile('redcells.mat')
            disp('Found manual red cell labels')
            load('redcells.mat', 'redcells');
            manual_redcells = redcells;
            redcell(:) = 0;
            redcell(manual_redcells + 1) = 1; %Add 1 for python to MATLAB transformation
        end

        block.redcell.(currentPlane) = redcell(keep_ind);
        if any(block.redcell.(currentPlane))
            disp('Found red cells')
        end

        %Update zcorr frame set
        if isfield (Fall.ops, 'zcorr')
            if n == 0
                block.zcorr = [];
            end
            block.zcorr.(currentPlane) = Fall.ops.zcorr(:,Plane_set); %Dimensions are z-stack position vs. frame
            disp('Found zcorr')
        end

        %Determine distances between neurons
        block.distance_in_microns.(currentPlane) = find_cell_distance(block, 'Plane', n);
    end
end

%% Bad frames correction
%TODO: Add a correction for badframes (e.g. interpolate)? If so, do this before any resampling

%% Upsample data to 30Hz if needed

if functionOps.UpsampleTo30 && block.setup.framerate < 30
    %If data has markpoints, upsample to 30 will happen in the define_markpoints_timestamp script
    if ~isfield(block.setup, 'hasMarkPoints')
        block = upsample_suite2p_data(block,30);
        disp('Data upsampled to 30Hz')
    end
end

end