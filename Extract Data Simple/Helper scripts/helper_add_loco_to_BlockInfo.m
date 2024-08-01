% Helper function:
% Add loco % and average speed to BlockInfo

ExtractedData_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice';

ExtractedDataFiles = {'ExtractedData_NoiseITI_20230412-001517.mat'...,
    'ExtractedData_FM_20230412-010119.mat'...,
    'ExtractedData_RF_20230412-031712.mat'...,
    'ExtractedData_SAMfreq_20230412-041442.mat'...,
    'ExtractedData_SAMnoise_20230412-031314.mat'};

for f = 1:length(ExtractedDataFiles)
    cd(ExtractedData_path);
    disp(ExtractedDataFiles{f})
    load(ExtractedDataFiles{f})
    
    BlockInfo = ExtractedData.BlockInfo;
    TrialData = ExtractedData.TrialData;
    
    BlockInfo = add_loco_to_BlockInfo(BlockInfo, TrialData);
    
    if f == 1
        StoreBlockInfo = BlockInfo;
    else
        StoreBlockInfo = [StoreBlockInfo; BlockInfo];
    end
end

%Remove refImg and meanImg from BlockInfo before saving excel spreadsheet
StoreBlockInfoCSV = StoreBlockInfo;
StoreBlockInfoCSV = removevars(StoreBlockInfoCSV,{'refImg', 'meanImgE'});

return %Maker user use next section intentionally

%% Save BlockInfo 

save('BlockInfoWithLoco.mat', 'StoreBlockInfo');
writetable(StoreBlockInfoCSV,'BlockInfoWithLoco.csv','Delimiter',',');
