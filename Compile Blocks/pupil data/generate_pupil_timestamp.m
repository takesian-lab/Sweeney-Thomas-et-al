function [block] = generate_pupil_timestamp(block)

% we compiled a lot of data before we had pupil timestamps. to avoid
% recompiling blocks, we created this script to generate a pupil timestamp
% and put it in a block. Eventually, you will not need this; however, it is
% a quick work around for the already compiled blocks. 

if isfield(block,'pupil_timestamp')
    return
end

if ~isfield(block,'PupilFrameData')
    disp('no pupil data for this block')
    return
end

%%

pupil_timestamp = [];
data = block.errorData.PupilFrameData;

if isfield(block, 'Sound_Time') % If you dont have FP or S2P data, you will not have a Sound_time variable

    error_trials = block.errorData.error_trials;
    if isempty(error_trials)
        for i = 1:length(block.Sound_Time)
            t = data{1,i}.frame_times;
            st = block.New_sound_times(i);
            Sound_Time = block.Sound_Time(i);
            M = Sound_Time - st;
            Pupil_dataTime{i} = t(:) + M;
            pupil_timestamp = [pupil_timestamp; Pupil_dataTime{i}(:)];
        end
    else        
        count = 0; %Count for non-error trials
        for i = 1:length(data)
            t = data{1,i}.frame_times;
            if ismember(i, error_trials)
                if i == 1
                    error('This method does not work if error trial is #1')
                end
                Mdiff = block.errorData.start_time(i) - block.errorData.start_time(i - 1); %Estimate diff between start time of error trial and previous trial
                M = M + Mdiff; %Add to previous trial's M
            else
                %Non-error trial sound times have been adjusted to match Bruker or FibPhot times
                count = count + 1;
                st = block.errorData.New_sound_times(i);
                Sound_Time = block.Sound_Time(count);
                M = Sound_Time - st;
            end
            Pupil_dataTime{i} = t(:) + M;
            pupil_timestamp = [pupil_timestamp; Pupil_dataTime{i}(:)];
        end
    end
    
else % generate a relative timestamp in absence of neural data
       
    for i = 1:length(data)
        t = data{1,i}.frame_times;
        
        %Correct nans (assuming there will just be 1 or 2)
        if any(isnan(t))
            nanind = find(isnan(t));
            fs = mean(diff(t),'omitnan');
            for n = 1:length(nanind)
                if nanind(n) == 1
                    t(nanind(n)) = 0;
                else
                    t(nanind(n)) = t(nanind(n)-1) + fs;
                end
            end
        end
                    
        if i > 1
            M = pupil_timestamp(end);
            if isnan(M)
                error('Timestamp has NaNs');
            end
        else
            M = 0;
        end

        Pupil_dataTime{i} = t + M;
        pupil_timestamp = [pupil_timestamp; Pupil_dataTime{i}];
    end
end

block.pupil_timestamp = pupil_timestamp;

end
