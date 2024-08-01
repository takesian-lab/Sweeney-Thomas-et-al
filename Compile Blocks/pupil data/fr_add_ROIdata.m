function ROIdata = fr_add_ROIdata(block, ROIdata, filename, N)
%Sometimes the u8 avi and Tosca data are fine, but the VideoConverter
%software runs into an issue with creating the video
%e.g. Write Frame error "Incompatible image size"
%If this happens, replace the whisker data with NaNs

disp(['......adding frame time to trial ' num2str(N) '......']);
originalROIdata = ROIdata;

trial = ROIdata.TraceByTrial{N};

    
%Figure out how many frames are missing:
nFrames = length(block.errorData.PupilFrameData{1,N}.frames);
n = nFrames - length(trial);

if n > 1 %ADD FRAMES
    %Replace with nans
    nandata = single(nan(1,nFrames-length(trial)));
    new_trial = [trial, nandata];

    %Replace
    ROIdata.TraceByTrial{N} = new_trial;
    
elseif n < 1 %REMOVE FRAMES
    new_trial = trial(1:end + n);

    %Replace
    ROIdata.TraceByTrial{N} = new_trial;
end

%Remake concat_trace
trace = [];
for i = 1:length(ROIdata.TraceByTrial)
    trace = [trace, ROIdata.TraceByTrial{i}];
end
ROIdata.Trace = trace;
    
%Save over old ROIdata so you don't have to do this again
save(strcat(filename, '.mat'), 'ROIdata')

%Save original ROIdata if the frame difference was large
if abs(n) > 3
    save(strcat(filename, '_original.mat'), 'originalROIdata')
end