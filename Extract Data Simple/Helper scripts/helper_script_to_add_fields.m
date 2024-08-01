% Helper function: Add Field1 and Field2 columns to ExtractedData

ExtractedData_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR z 1 reliability blanks';

ExtractedDataFiles = {'ExtractedData_FM_20230119-135425.mat',...
    'ExtractedData_NoiseITI_20230116-171548.mat'...,
    'ExtractedData_RF_20230116-212251.mat'...,
    'ExtractedData_SAM_20230117-135226.mat'...,
    'ExtractedData_SAMfreq_20230117-164541.mat'};

for f = 1:length(ExtractedDataFiles)
    cd(ExtractedData_path);
    disp(ExtractedDataFiles{f})
    load(ExtractedDataFiles{f})
    
    Ops = ExtractedData.Ops;
    
    %Load info sheet and select files corresponding to stim_protocol
    Info = tak_read_info_table(Ops.info_file_path, 'StimProtocol', Ops.stim_protocol, 'StripFields', false);

    %Recreate BlockInfo
    [BlockInfo, ~, ~] = fillSimpleDataTableFromInfo(Info, Ops.blocks_path, 'StimProtocol', Ops.stim_protocol, 'ReturnBlocks', 0, 'ReturnCellList', 0);

    %Add Field1 and Field2 variable to all of ExtractedData
    for b = 1:height(BlockInfo)
        block = BlockInfo.Block(b);
        field1 = BlockInfo.Field1(b);
        field2 = BlockInfo.Field2(b);
        
        %BlockInfo
        ind = strcmp(block, ExtractedData.BlockInfo.Block);
        ExtractedData.BlockInfo.Field1(ind) = field1;
        ExtractedData.BlockInfo.Field2(ind) = field2;
    
        %CellList
        ind = strcmp(block, ExtractedData.CellList.Block);
        ExtractedData.CellList.Field1(ind) = field1;
        ExtractedData.CellList.Field2(ind) = field2;

        %Summary
        loco_fields = fields(ExtractedData.Summary);
        for L = 1:length(loco_fields)
            ind = strcmp(block, ExtractedData.Summary.(loco_fields{L}).Block);
            ExtractedData.Summary.(loco_fields{L}).Field1(ind) = field1;
            ExtractedData.Summary.(loco_fields{L}).Field2(ind) = field2;
        end
    end
        
    %Save new ExtractedData file
    disp('Saving ExtractedData...')
    cd(ExtractedData_path);
    filename = ExtractedData.Filename;
    save([filename(1:end-4) '_edited.mat'], 'ExtractedData', '-v7.3');
    clear('ExtractedData')
end
