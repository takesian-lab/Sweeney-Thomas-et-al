function visualize_extracted_ephys_block(ExtractedData)
%FUNCTION IN PROGRESS
%it shows individual traces for ephys data after extracting, and mean +-
%SEM of those traces


%%
%Obtain 'time' variable to plot in x axis from frame rate and number of frames:
trace_duration = ExtractedData.StimInfo.nFrames/ExtractedData.StimInfo.fs;
step = (ExtractedData.StimInfo.nFrames/ExtractedData.StimInfo.fs)/ExtractedData.StimInfo.nFrames;
time = 0:step:trace_duration-step;
%%
%Find out which row in ExtractedData.CellDataByStim.All(i).StimTraces
%contains a stim (which could be identified in
%ExtractedData.StimInfo.Combined_Label)
numRows = length(ExtractedData.CellDataByStim.All);
nonEmptyRow_index = zeros(numRows,1);
for i = 1:numRows
    numRows2 = size(ExtractedData.CellDataByStim.All(i).StimTraces, 1);   
    for j = 1:numRows2
        if ~isempty(ExtractedData.CellDataByStim.All(i).StimTraces{j,1})            
            nonEmptyRow_index(i) = j;
            stimTitle(i,1) = ExtractedData.StimInfo.Combined_Label(j);
        end
    end
end


%%
%Plot all individual traces for each stim:
figure;
for i = 1:numRows    
    subplot(round(numRows/2),2,i)           
         numRows3 = size(ExtractedData.CellDataByStim.All(i).StimTraces{nonEmptyRow_index(i),1},1);
         for j = 1:numRows3       
            hold on
            plot(time,ExtractedData.CellDataByStim.All(i).StimTraces{nonEmptyRow_index(i),1}(j,:))
            vline(0.5,'r')
            xlabel('time (sec)')
            ylabel('Z-score')
            str = strcat('Row', num2str(i), '-', stimTitle(i));
            title(str);
         end       
end
%%
%Plot traces average for each stim with SEM:
figure;
for i = 1:numRows
    subplot(round(numRows/2),2,i)
    temp_data = ExtractedData.CellDataByStim.All(i).StimTraces{nonEmptyRow_index(i)};
    means = mean(temp_data,1); % calculate means of each column
    temp_data_with_mean = vertcat(temp_data,means);
    sem = std(temp_data)/sqrt(size(temp_data,1)); % calculate SEM
    hold on
    plot(time,temp_data_with_mean(end,:),'b','LineWidth',2);
    shadedErrorBar(time,temp_data_with_mean(end,:),sem);
    
    vline(0.5,'r')
    xlabel('time (sec)')
    ylabel('Z-score')
    str = strcat('Row', num2str(i), '-', stimTitle(i));
    title(str);
end

%%
%clearvars i j nonEmptyRow_index numRows numRows2 numRows3 step str time trace_duration
