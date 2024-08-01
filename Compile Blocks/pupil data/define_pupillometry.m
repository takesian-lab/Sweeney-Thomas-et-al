function [block] = define_pupillometry(block)
% This function reads DeepLabCut .csv files and extracts pupil radius 
%
% Version history:
% - V1 = current version
%
%
% TODO: Confirm how to store concatenated pupil radius and corresponding
% time to block, add pupil radius to align_to_stim, figure out trial start
% time with respect to pupillometry timescale
% Search 'TODO'

%% Skip this function if pupil data is not available

if ~isfield(block, 'PupilFrameData')
    return
end
    
disp('Pulling out Pupillometry data...');

%% Go to Tosca folder and look for DLC .csv files

tosca_path = block.setup.Tosca_path;

%If enhanced videos exist, cd to enhanced videos folder
if isfolder(strcat(tosca_path, '\enhanced videos'))
    dlc_path = strcat(tosca_path, '\enhanced videos');
    cd(dlc_path);
else
    cd(tosca_path);
end

%Look for DLC .csv files
csvfiles = dir('*.csv'); %all of the csv files in directory
csvfiles = {csvfiles.name}'; %transform to cell array

%Only keep .csv files for the current Tosca session and run
tosca_csv_name = strcat(block.setup.mousename, '-Session', num2str(block.setup.Tosca_session), '-Run', num2str(block.setup.Tosca_run));
csvfiles = csvfiles(contains(csvfiles, tosca_csv_name));

%Only keep .csv files for the selected DLC training algorithm
%We only have one for now, but we may add more in the future
DLC_training_name = 'DLC_resnet101_Pupillometry_01_15_23Jan15shuffle1_500000';
csvfiles = csvfiles(contains(csvfiles, DLC_training_name));

% filter to only include block of interest:
tosca_data = strcat('Session',string(block.setup.Tosca_session),'-Run',string(block.setup.Tosca_run),'.');
csvfiles = csvfiles(contains(csvfiles, tosca_data));

%If there are no files left after all this, skip processing for now
if isempty(csvfiles)
    disp('No DeepLabCut files found. Skipping...')
    return
end

%% Pull out pupil radius and eyelid area

disp('DeepLabCut files found. Computing pupil size...')
csvfiles(1) = []; %000 trial, do not include in data
[Radius, ~] = get_pupil_radius(csvfiles, DLC_training_name, 0);

%% Save in block

%Before removing errors, there should be the same number of trials as csv files
all_data = block.errorData.PupilFrameData;
error_trials = block.errorData.error_trials;

if length(Radius) ~= length(all_data)
    error('Number of trials do not match')
end

%Check that number of frames match for all trials and correct differences

%Temporarily record frame numbers for troubleshooting
FrameRecord = table;
FrameRecord.TrialNumber = (1:length(all_data))';
FrameRecord.Radius = nan(length(all_data),1);
FrameRecord.Tosca = nan(length(all_data),1);
for i = 1:length(all_data)
    RadiusLength = length(Radius{i});
    ToscaFrameLength = length(all_data{1,i}.frames);
    
    FrameRecord.Radius(i) = RadiusLength;
    FrameRecord.Tosca(i) = ToscaFrameLength; 
end
FrameRecord.FrameDiff = abs(FrameRecord.Radius - FrameRecord.Tosca); %TIP: Sort FrameRecord by FrameDiff to find all non-matching trials

%Perform fixes if necessary
for i = 1:length(all_data)
    RadiusLength = length(Radius{i});
    ToscaFrameLength = length(all_data{1,i}.frames);
    
    if RadiusLength ~= ToscaFrameLength
        if RadiusLength - ToscaFrameLength == 1
            %Tosca data is missing a frame, add back with this function:
            block = fr_add(block,i);
        elseif RadiusLength < ToscaFrameLength
            error('Pupil data is missing frames. Check that avi is not corrupt. If not, manually correct DLC csv')
        elseif ismember(i, error_trials)
            %Error trials are more likely to be off in frame numbers
            %Here, allow to correct for these larger differences
            if RadiusLength > ToscaFrameLength
                block = fr_add(block,i,RadiusLength - ToscaFrameLength);
            end
        else
            %If the trial was not an error, look into it:
            %Check video and Tosca files for any discrepancies
            error('mismatch of DLC and tosca frames')
        end
    end
end

%Now add corrected pupil data to PupilFrameData for non-error trials only
all_data = block.errorData.PupilFrameData; %Refresh in case errorData got corrected in functions above
trials = 1:length(all_data);
trials(error_trials) = [];

for i = 1:length(trials)
    T = trials(i);
    currentTrace = Radius{T};
    RadiusLength = length(currentTrace);
    ToscaFrameLength = length(block.PupilFrameData{1,i}.frames);
    %Catch any remaining problems
    if RadiusLength ~= ToscaFrameLength
        error('Frames lengths are still mismatched')
    end
    block.PupilFrameData{1,i}.radius = single(currentTrace);
end

%Store concatenated pupil radius in block
concat_pupil = [];
for i = 1:length(Radius)
    concat_pupil = [concat_pupil, Radius{i}'];
end

%Make sure concatenated DLC trace matches pupil timestamp
%WARNING: pupil_timestamp will change again during define_sound_fibPhot,
%but this allows us to catch any errors with DLC data earlier
[block] = generate_pupil_timestamp(block);
if length(block.pupil_timestamp) ~= length(concat_pupil)
    error('Timestamps do not match')
end

%Store concatenated pupil data in block
block.concat_pupil = single(concat_pupil);
