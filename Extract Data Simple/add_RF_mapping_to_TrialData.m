function TrialData = add_RF_mapping_to_TrialData(StimInfo, TrialData)
% Add StimID field to TrialData which corresponds to stim order in StimInfo.Combined
% StimInfo.Combined is a vertical list of all stim combinations [V1,V2]
% Blanks will be coded as StimID 0

for b = 1:size(TrialData,2)
    
    %Get trial information from TrialData
    isBlank = TrialData(b).IsBlank;
    stim_v1 = TrialData(b).Stim_V1;
    stim_v2 = TrialData(b).Stim_V2;

    stim_id = nan(size(isBlank)); %store stimID (number corresponding to StimInfo.Combined)
    
    for v = 1:size(StimInfo.Combined,1)
        v1 = StimInfo.Combined(v,1);
        v2 = StimInfo.Combined(v,2);

        %Find rows corresponding to current stim
        if any(StimInfo.StimProtocol == [21 211])
            stim_rows = intersect(find(strcmp(stim_v1,v1)),find(strcmp(stim_v2,v2)));
        else
        stim_rows = intersect(find(stim_v1 == v1), find(stim_v2 == v2));
        end
            
        %Record stimID
        stim_id(stim_rows) = v;
    end

    stim_id(isBlank) = 0;

    if any(isnan(stim_id))
        if StimInfo.StimProtocol ~=7
        error('Troubleshoot')
        end
    end
    
    TrialData(b).StimID = stim_id;
end