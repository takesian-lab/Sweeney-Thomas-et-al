function [block] = define_whiskerpad(block)
% Set whiskerpad ROI on video and extract average pixel difference for each frame
% This will prompt a user to set the ROI the first time a block is compiled and save
% ROI.mat in the Tosca folder, which will be automatically loaded on subsequent runs
%
% Version history:
% - V1 = current version
%
%
% TODO: 
% Search 'TODO'

%% Skip this function if video data is not available

if ~isfield(block, 'PupilFrameData')
    return
end

if block.setup.skip_avi
    warning('whisker pad data not computed because user chose to skip avi analysis in tosca')
    return
end
    
disp('Pulling out whisker pad data...');

look_for_previous_ROIdata = 1;

%% Go to Tosca folder and load ROIdata.mat file if it exists
% Otherwise, create new ROIdata.mat to define ROI and pull out trace

tosca_path = block.setup.Tosca_path;
cd(tosca_path);

tosca_name = strcat(block.setup.mousename, '-Session', num2str(block.setup.Tosca_session), '-Run', num2str(block.setup.Tosca_run),'.');
filename = strcat('ROIdata_', tosca_name);
filename = erase(filename,".");

if look_for_previous_ROIdata && isfile(strcat(filename, '.mat'))
    disp('ROI data found.');
    load(strcat(filename, '.mat')); %loads ROIdata variable
else
    %Look for .avi files
    %**Enhanced videos may exist, for now use the raw videos for whisker pad movement
    avifiles = dir('*.avi'); %all of the avi files in directory
    avifiles = {avifiles.name}'; %transform to cell array

    %Only keep .avi files for the current Tosca session and run
    avifiles = avifiles(contains(avifiles, tosca_name));

    %If there are no files left after all this, skip processing for now
    if isempty(avifiles)
        disp('No videos found. Skipping...')
        return
    end

    % Extract ROIdata and save in Tosca folder
    avifiles(1) = []; %000 trial, do not include in data
    
    save_ROI_figure = 1; %Save figure with ROI in Tosca folder for reference
    [ROIdata] = get_orofacial_ROIs_takesian(avifiles, filename, save_ROI_figure, 0);
    save(strcat(filename, '.mat'), 'ROIdata');
end

%% Save in block

%Before removing errors, there should be the same number of trials as ROIData
all_data = block.errorData.PupilFrameData;
error_trials = block.errorData.error_trials;

if length(ROIdata.TraceByTrial) ~= length(all_data)
    error('Number of trials do not match')
end

%Check that number of frames match for all trials and correct differences

%Temporarily record frame numbers for troubleshooting
FrameRecord = table;
FrameRecord.TrialNumber = (1:length(all_data))';
FrameRecord.ROI = nan(length(all_data),1);
FrameRecord.Tosca = nan(length(all_data),1);
for i = 1:length(all_data)
    ROILength = length(ROIdata.TraceByTrial{i});
    ToscaFrameLength = length(all_data{1,i}.frames);
    
    FrameRecord.ROI(i) = ROILength;
    FrameRecord.Tosca(i) = ToscaFrameLength; 
end
FrameRecord.FrameDiff = abs(FrameRecord.ROI - FrameRecord.Tosca); %TIP: Sort FrameRecord by FrameDiff to find all non-matching trials

%Perform fixes if necessary

%If pupil data already exists in block, compare whisker data to pupil data
%define_pupillometry should have already fixed any mismatches between Tosca and avi timestamps
%If any discrepancies remain, it means there was a mismatch between DLC and ROIdata
%Review that here and choose to favor DLC data as long as the discrepancy isn't too big
if isfield(block, 'concat_pupil')
    
    %Compare all ROI data to DLC data
    for i = 1:length(all_data)
        ROILength = length(ROIdata.TraceByTrial{i});
        RadiusLength = length(block.errorData.PupilFrameData{1,i}.frames); %Do not use block.PupilFrameData.radius (missing error trials)

        if ROILength ~= RadiusLength
            if ROILength - RadiusLength == 1
                %ROI data is longer by 1 frame, correct to favor pupillometry data:
                ROIdata = fr_add_ROIdata(block, ROIdata, filename, i);
            elseif ROILength < RadiusLength
                error('This has not happened yet, code for it if needed')
                %If you KNOW that pupil and ROIdata do not match due to a
                %corrupted video file and you already manually fixed the pupil data,
                %you can use fr_add_ROIdata to make ROIdata match the pupil data: 
                %ROIdata = fr_add_ROIdata(block, ROIdata, filename, i);
            elseif ismember(i, error_trials)
                error('This has not happened yet, code for it if needed')
            else
                %If the trial was not an error, look into it:
                %Check video and Tosca files for any discrepancies
                error('mismatch of whisker and pupil data')
            end
        end
    end
else
    %If pupil data does not exist, proceed with fixing Tosca data to match ROI data
    for i = 1:length(all_data)
        ROILength = length(ROIdata.TraceByTrial{i});
        ToscaFrameLength = length(all_data{1,i}.frames);

        if ROILength ~= ToscaFrameLength
            if ROILength - ToscaFrameLength <= 3 %Allow for very small errors in frame length
                %Tosca data is missing a frame, add back with this function:
                block = fr_add(block,i,ROILength - ToscaFrameLength);
            elseif ROILength < ToscaFrameLength
                %Whisker data is missing frames, replace with NaNs:
                ROIdata = fr_add_ROIdata(block, ROIdata, filename, i);
            elseif ismember(i, error_trials)
                %Error trials are more likely to be off in frame numbers
                %Here, allow to correct for these larger differences
                if ROILength > ToscaFrameLength
                    block = fr_add(block,i,ROILength - ToscaFrameLength);
                end
            else
                %If the trial was not an error, look into it:
                %Check video and Tosca files for any discrepancies
                error('mismatch of AVI and tosca frames')
            end
        end
    end
end

%Now add corrected whiskerpad data to PupilFrameData for non-error trials only
all_data = block.errorData.PupilFrameData; %Refresh in case errorData got corrected in functions above
trials = 1:length(all_data);
trials(error_trials) = [];

for i = 1:length(trials)
    T = trials(i);
    currentTrace = ROIdata.TraceByTrial{1,T};
    ROILength = length(currentTrace);
    ToscaFrameLength = length(block.PupilFrameData{1,i}.frames);
    %Catch any remaining problems
    if ROILength ~= ToscaFrameLength
        error('Frames lengths are still mismatched')
    end
    block.PupilFrameData{1,i}.whiskerpad = single(currentTrace');
end

%Make sure concated whisker trace matches pupil trace
%WARNING: pupil_timestamp will change again during define_sound_fibPhot,
%but this allows us to catch any errors with whisker data earlier
concat_whisker = ROIdata.Trace; %This should be corrected along with the other data during fr_add_ROIdata
[block] = generate_pupil_timestamp(block);
if length(block.pupil_timestamp) ~= length(concat_whisker)
    error('Timestamps do not match')
end

%Store concatenated ROI data in block
block.concat_whisker = single(concat_whisker);