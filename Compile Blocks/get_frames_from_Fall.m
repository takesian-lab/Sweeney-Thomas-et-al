function [Frame_set, Plane_set] = get_frames_from_Fall(ops, block_name, displayTable, plane) 
% This function gets frame numbers for each recording block from the Suite2p Fall.mat file
% Returns Frame_set and/or Plane_set for the block specified by block_name.
% Plane_set applies to multiplane data only
% Optionally displays a complete table of block names and frame numbers for the current Fall.mat
%
% Argument(s): 
%   ops (struct from Fall.mat)
%   block_name (optional string) - if no block_name is included, function will display table  without returning anything
%   displayTable (optional 0 or 1) - option to display table
%   plane (optional integer corresponding to suite2p plane - e.g. for plane0  enter 0)
%   *If no plane is specified for multiplane data, Frame_set will still be returned
% 
% Returns:
%   Frame_set - full list of frames in current block
%   Plane_set - full list of frames in current cycle (multiplane data only)
% 
% Notes:
%
%
% TODO: 
% Search 'TODO'

%% Define function options

%If no block_name is specified, simply display table
returnFrameSet = 1;
if nargin < 2
    returnFrameSet = 0;
    displayTable = 1;
end

%Default is to display block information in table format
if nargin < 3
    displayTable = 1;
end

%Accommodate multiplane data
if nargin < 4
    plane = nan;
    multiplaneData = false;
else
    plane = plane + 1; %python to matlab conversion
    multiplaneData = true;
end

%% Determine if multipage tiffs

multipageTIFFs = false;
if isfield(ops, 'frames_per_file') %Old suite2p data do not have this variable (but therefore are also not going to be multipage)
    if any(ops.frames_per_file > 1)
        multipageTIFFs = true;
    end
end

if multipageTIFFs && multiplaneData
    error('This has not been coded for yet')
end

%% Get cycle, frame, block, and channel numbers from the filenames in ops

%The filelist format should be as follows:
%D:/FILEPATH/BOT_FILENAME-###/BOT_FILENAME-###_Cycle#####_Ch#_######.ome.tif
% # is the channel number (1 or 2)
% ### is the block number
% ##### is the cycle number (used for t-series)
% ###### is the tif/frame number

frameA = 13; %number of chars to start of frame number
frameB = 8;  %number of chars to end of frame number
blockA = 32; %number of chars to start of block number
blockB = 30; %number of chars to end of block number
chanA  = 15; %number of chars to channel number
cycleA = 23; %number of chars to cycle number
cycleB = 19; %number of chars to end of cycle number

redChannel = 1; %Number of the red channel, i.e. Ch1 or Ch2

% Get filelist from Fall.mat:
filelist_char = ops.filelist;
filelist_str = string(ops.filelist);

[frames, blocks, chans, cycles] = deal(nan(size(filelist_str))); %Prep variables
paths = strings(size(filelist_str)); %Prep variables

%Process filenames with the same length at the same time (this will be faster than looping through all files)
filelist_lengths = strlength(deblank(filelist_str)); %Get length of each filename
unique_lengths = unique(filelist_lengths); %How many unique lengths are there

for i = 1:length(unique_lengths)
    L = unique_lengths(i);
    currentFiles = filelist_lengths == L;
    frames(currentFiles,:) = double(string(filelist_char(currentFiles,L-frameA:L-frameB)));
    blocks(currentFiles,:) = double(string(filelist_char(currentFiles,L-blockA:L-blockB)));
    paths(currentFiles,:) = string(filelist_char(currentFiles,1:L-blockB));
    chans(currentFiles,:) = double(string(filelist_char(currentFiles,L-chanA)));
    cycles(currentFiles,:) = double(string(filelist_char(currentFiles,L-cycleA:L-cycleB)));
end

%If data has both red and green channels, delete channel 1 (F only includes green frames)
if ops.nchannels > 1 % keep red channel if using RCaMP or other red indicator. 
    if sum(chans == redChannel) > 0
        frames(chans == redChannel,:) = [];
        blocks(chans == redChannel,:) = [];
        paths(chans == redChannel,:) = [];
        cycles(chans == redChannel,:) = [];
    end
end

%% For multiplane data, correct frame numbers that are > nPlanes
%This happens when the last cycle has extra frames
%e.g. Last cycle has 2 extra frames: 1 2 3 / 1 2 3 / ... / 1 2 3 4 5
if multiplaneData
    %In single plane data, frames count up. In multiplane data, frames represent planes and cycles count up
    T = tabulate(frames); 
    nFrames = mode(T(:,2)); %Frames per plane
    nPlanes = sum(T(:,2) == nFrames); %Number of planes with an equal number of frames

    %SUITE2P EXTRA PLANE PROBLEM
    %The following section is commented out because Suite2p does not assign
    %extra frames as we would expect. Since we cannot be sure about the identify of the frames
    %suite2p assigned instead, we can't correct the frame numbers
    %Instead, there is a catch in define_suite2p_singleblock that will edit
    %timestamp to only include the frames we know the identity of
    %https://docs.google.com/presentation/d/1D32WRyEAAMaBo9NYC7Q9kB_p9UuGU262AKHHhZLscQs/edit?usp=sharing
            
%     %Identify and correct extra frames
%     for i = 1:size(T,1)
%         if T(i,1) > nPlanes
%             %Replace extra frames with Frame - nPlanes
%             %CAUTION: This will not work for bidirectional data, where
%             %frames are 1 2 3 / 3 2 1 / 1 2 3 ...... but bidirectional data
%             %should not be run through suite2p. Use the code convert_bidirectional_tiffs first
%             currentPlane = T(i,1);
%             ind = find(frames == currentPlane);
%             frames(ind) = frames(ind) - nPlanes;
%         end
%     end
end

%% Get first and last frame of each block by finding unique filepaths

uniquePaths = unique(paths); %Unique paths corresponds to # blocks

%Order paths by block number
uniquePathBlocks = nan(size(uniquePaths));
for i = 1:length(uniquePaths)
    pathBlocks = blocks(strcmp(uniquePaths{i},paths));
    uniquePathBlocks(i) = unique(pathBlocks);
end
[blockNumbers, sortInd] = sort(uniquePathBlocks);
uniquePaths = uniquePaths(sortInd);

%Find first and last frame associated with each block
%Use ind as opposed to frame # because frames are counted differently for single plane and multiplane data
blockFrames = nan(length(uniquePaths),2); %Prep variable
for i = 1:length(uniquePaths)
    currentPath = uniquePaths(i);
    matchingPaths = strcmp(paths, currentPath);
    blockFrames(i,1) = find(matchingPaths, 1, 'first');
    blockFrames(i,2) = find(matchingPaths, 1, 'last');
end

%If multipage tiffs found, adjust frame numbers (Currently this only works with single plane data)
if multipageTIFFs
   newBlockFrames = nan(length(uniquePaths),2); %Prep variable
   for i = 1:length(uniquePaths)
        currentPath = uniquePaths(i);
        matchingPaths = strcmp(paths, currentPath);
        nBlockFrames = sum(ops.frames_per_file(matchingPaths));
        
        firstMatchingPath = find(matchingPaths, 1, 'first');
        if firstMatchingPath > 1
            nPreviousFrames = sum(ops.frames_per_file(1:firstMatchingPath-1));
        else
            nPreviousFrames = 0;
        end
        
        newBlockFrames(i,1) = nPreviousFrames + 1;
        newBlockFrames(i,2) = nPreviousFrames + nBlockFrames;
    end 
    blockFrames = newBlockFrames;
end

%Remove blocks that are only one frame long (this could be from a SingleImage file in the same folder as the other data)
%I don't think this should happen unlesss you used the 'Look one down' option in suite2p
%or you accidentally including a SingleImage folder instead of a block, but I'll leave here for now
isOneFrame = (blockFrames(:,2) - blockFrames(:,1)) == 0;
blockNumbers(isOneFrame,:) = [];
blockFrames(isOneFrame,:) = [];
uniquePaths(isOneFrame,:) = [];

%% Get first and last frame of each block for a specific plane (plane_set)

if multiplaneData
    %Only look at current plane
    planeFrames = frames(frames == plane);
    planeCycles = cycles(frames == plane);
    planeBlocks = blocks(frames == plane);
    planePaths = paths(frames == plane);
    
    blockPlaneFrames = nan(length(uniquePaths),2); %Prep variable
    for i = 1:length(uniquePaths)
        currentPath = uniquePaths(i);
        matchingPaths = strcmp(planePaths, currentPath);
        blockPlaneFrames(i,1) = find(matchingPaths, 1, 'first');
        blockPlaneFrames(i,2) = find(matchingPaths, 1, 'last');
    end
    
else
    blockPlaneFrames = nan;
end
    

%% Get block names by splitting filepaths based on both / and \ and return frame set
if size(uniquePaths,1) == 1
    blockNames_temp1 = split(uniquePaths, '/');
    blockNames_temp2 = split(blockNames_temp1(end), '\');
    blockNames = blockNames_temp2(end);
else
    blockNames_temp1 = split(uniquePaths, '/');
    blockNames_temp2 = split(blockNames_temp1(:,end), '\');
    blockNames = blockNames_temp2(:,end);
end

%Display table of block numbers and frames to user
if displayTable
    format long
    if multiplaneData
        disp(table(blockNames, blockNumbers, blockFrames, blockPlaneFrames))
    else
        disp(table(blockNames, blockNumbers, blockFrames))
    end
end

%Convert blockFrames and blockPlaneFrames into Frame_set and Plane_set based on block_name
if returnFrameSet

    %Accommodate Z-corrected and converted bidirectional (1D) data
    %and multiplane data that has been converted to single plane
    %The reason we have to do this is because we added a prefix to the
    %block folder name during conversion, but this doesn't match the Tiffs inside the folder
    %If length(char_block_name) is less than the number of characters we're
    %comparing it to, there will be an error. Pad the end of the comparison string to avoid this
    char_block_name = char(block_name);
    padded_char_block_name = [char_block_name,'STRINGPADDING'];
    if isequal(padded_char_block_name(1:3),'1D-')
        block_name = char_block_name(4:end);
    elseif isequal(padded_char_block_name(1:14),'Zcorrected-1D-')
        block_name = char_block_name(15:end);
    elseif isequal(padded_char_block_name(1:11),'Zcorrected-')
        block_name = char_block_name(12:end);
    elseif isequal(padded_char_block_name(1:5),'Plane')
        block_name = char_block_name(8:end); %If plane# is ever >9 this will break
    end
        
    matching_blocks = strcmp(blockNames, block_name);
    
    if sum(matching_blocks) == 0
        error('block_name is not contained in dataset')
    elseif sum(matching_blocks) > 1
        error('Multiple blocks with the same block_name');
    end
    currentFrames = blockFrames(matching_blocks,:);
    Frame_set = currentFrames(1,1):currentFrames(1,2);
    
    if multiplaneData
        currentPlaneFrames = blockPlaneFrames(matching_blocks,:);
        Plane_set = currentPlaneFrames(1,1):currentPlaneFrames(1,2);
        disp(['Current plane set is:   ' num2str(currentPlaneFrames)]) 
    else
        Plane_set = nan;
        disp(['Current frame set is:   ' num2str(currentFrames)])
    end
end


