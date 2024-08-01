function block = fr_add(block,i,nFrames)
% tosca occasionally failes to record a timestamp for a single frame. this
% function will add a time to a single frame and update the pupil time
% stamps. 
% input:    block
%           error_i ....trial index in block.errorData.PupilFrameData
% output: block with corrected timestamps and trial numbers

disp(['......adding frame time to trial ', num2str(i), '......']);

if nargin < 3
    nFrames = 1;
end

ftime = block.errorData.PupilFrameData{1,i}.frame_times;
frames = block.errorData.PupilFrameData{1,i}.frames;
fs = mean(diff(ftime));

%Ind where new data will be added (if nFrames is 1, A and B will be equal)
A = length(ftime) + 1;
B = length(ftime) + nFrames;

%Pad end of array with the number of frames needed
nan_add = nan(nFrames,1);
ftime = [ftime; nan_add];
frames = [frames; nan_add];

%Fill array with new timestamps and trial numbers
ftime(A:B) = ftime(A-1)+fs:fs:ftime(A-1)+(nFrames*fs);
frames(A:B) = frames(A-1)+1:frames(A-1)+nFrames;

%replace error data
block.errorData.PupilFrameData{1,i}.frame_times = ftime;
block.errorData.PupilFrameData{1,i}.frames = frames;

%replace trial data
%Figure out the trial number
trials = 1:length(block.errorData.PupilFrameData);
trials(block.errorData.error_trials) = [];
ii = find(trials == i);

if ~isempty(ii)
    block.PupilFrameData{1,ii}.frame_times = ftime;
    block.PupilFrameData{1,ii}.frames = frames;
end

%Remove previous pupil_timestamp and regenerate
if isfield(block, 'pupil_timestamp')
    block = rmfield(block, 'pupil_timestamp');
end
block = generate_pupil_timestamp(block);

end