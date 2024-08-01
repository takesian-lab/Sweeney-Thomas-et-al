function StimInfo = add_RF_mapping_to_StimInfo(StimInfo, combined_stim)
% Combined_stim is a vertical list of all stim combinations [V1,V2]

%Find unique stim combinations
unique_rows = unique(combined_stim, 'rows');

%Remove NaNs
if ~any(StimInfo.StimProtocol == [21 211]) % V1/V2 are strings
if any(any(isnan(unique_rows)))
    %Since the unique function doesn't treat NaNs like numbers, we will
    %temporarily recode NaNs as 99 to find unique stim combos and then change back after
    if any(any(unique_rows == 99))
        %If 99 is one of the stim identities, this won't work
        error('Use different number to recode NaNs')
    end
    
    %Temporarily recode for NaNs using 99
    unique_rows(isnan(unique_rows)) = 99;
    unique_rows = unique(unique_rows, 'rows');
    
    %Change back to NaNs
    unique_rows(unique_rows == 99) = nan;
end


%Record unique stim
StimInfo.V1 = single(unique(unique_rows(:,1)));
StimInfo.V2 = single(unique(unique_rows(:,2)));
else
    StimInfo.V1 = unique(unique_rows(:,1)); % Stimprotocol 21/211 cant be stored as a single
    StimInfo.V2 = unique(unique_rows(:,2)); % Stimprotocol 21/211 cant be stored as a single
end 




%Order dB from loudest to softest
if any(strcmp(StimInfo.Units, 'dB'))
    if find(strcmp(StimInfo.Units, 'dB')) == 1
        StimInfo.V1 = flipud(StimInfo.V1);
    elseif  find(strcmp(StimInfo.Units, 'dB')) == 2
        StimInfo.V2 = flipud(StimInfo.V2);
    end
end

%Record combination of unique stim as final Combined list
%unique_rows might not contain all stim combinations, so we will create them here

if any(StimInfo.StimProtocol == [21 211]) % V1/V2 are strings and linked. 
    Combined = unique_rows;
else
Combined = nan(length(StimInfo.V1)*length(StimInfo.V2),2);
count = 1;
for vv = 1:length(StimInfo.V2)
    for v = 1:length(StimInfo.V1)
        Combined(count,1) = StimInfo.V1(v);
        Combined(count,2) = StimInfo.V2(vv);
        count = count + 1;
    end
end
end
StimInfo.Combined = Combined;

%Record Labels for each stim (strings that can be used in figure axes)
StimInfo.Combined_Label = strings(size(StimInfo.Combined,1),1);
StimInfo.V1_Label = strings(size(StimInfo.Combined,1),1);
StimInfo.V2_Label = strings(size(StimInfo.Combined,1),1);
for i = 1:length(StimInfo.Combined_Label)
    StimInfo.V1_Label(i) = strcat(num2str(StimInfo.Combined(i,1)), StimInfo.Units{1});
    StimInfo.V2_Label(i) = strcat(num2str(StimInfo.Combined(i,2)), StimInfo.Units{2});
    StimInfo.Combined_Label(i) = strcat(StimInfo.V1_Label(i), {' '}, StimInfo.V2_Label(i));
end

%Determine receptive field mapping based on V1 and V2
if any(StimInfo.StimProtocol == [21 211]) % V1/V2 are and always matched
    StimInfo.RF_V1 = nan(length(StimInfo.V1),1);
    StimInfo.RF_V2 = nan(length(StimInfo.V2),1);
    StimInfo.RF_Map = nan(size(unique_rows,1),2);
else
StimInfo.RF_V1 = single(nan(length(StimInfo.V2),length(StimInfo.V1)));
StimInfo.RF_V2 = single(nan(length(StimInfo.V2),length(StimInfo.V1)));
StimInfo.RF_Map = single(nan(size(unique_rows,1),2));

count = 1;
for vv = 1:length(StimInfo.V2)
    for v = 1:length(StimInfo.V1)
        StimInfo.RF_V1(vv,v) = StimInfo.V1(v);
        StimInfo.RF_V2(vv,v) = StimInfo.V2(vv);
        StimInfo.RF_Map(count,1) = vv; %Y
        StimInfo.RF_Map(count,2) = v; %X
        count = count + 1;
    end
end
end

