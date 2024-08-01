function BlockInfo = add_loco_to_BlockInfo(BlockInfo, TrialData)
%Add percent of time spent running and average speed per block to BlockInfo
%This info could be used to control for differences in running between groups

%**Loco from TrialData (instead of block) will already be cropped to
%match the 2P or FibPhot trace so we won't be including the parts of the loco
%trace that happen before/after functional data**

%% Loop through each block

for b = 1:height(BlockInfo)

    %Determine percentage of time mouse is running
    [loco_percent, loco_mean] = get_running_percent_and_speed(TrialData(b));

    %BlockInfo
    BlockInfo.LocoPercent(b) = loco_percent;
    BlockInfo.LocoSpeed(b) = loco_mean;
end