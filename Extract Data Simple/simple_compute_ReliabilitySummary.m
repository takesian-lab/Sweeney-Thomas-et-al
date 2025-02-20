function ReliabilitySummary = simple_compute_ReliabilitySummary(CellDataByStim)
%CellDataByStim is for a SINGLE CELL
%Example of how to run: ReliabilitySummary = simple_compute_ReliabilitySummary(CellDataByStim.(Loco{L})(c));
%
%Version control:
%V1 = archived 2/27/24 by Maryse
%   -R1 = average R for stim with sig. peaks
%   -R2 = average R for reliable stim with sig. peaks
%   -R3 = average R for non-reliable stim with sig. peaks
%V2 = current version (updated to save more meaningful R values)
%   -R1 = average R for all stim
%   -R2 = average R for isRF stim
%% 

ReliabilitySummary = table; 

%Check for reliability data

if isfield(CellDataByStim, 'ReliabilityData')
    ReliabilityType = "Pearsons";
elseif isfield(CellDataByStim, 'XcorrData')
    ReliabilityType = "XCorr";
else
    %No reliability data, return empty table
    return;
end

%% Compute average R values for ReliabilityData

isRF = CellDataByStim.RF == 1; %Convert nans to zeros

if ReliabilityType == "Pearsons"
    
    R = [CellDataByStim.ReliabilityData(:).R]';
    ReliabilitySummary.R1 = mean(R,'omitnan');
    ReliabilitySummary.R2 = mean(R(isRF),'omitnan');
        
elseif ReliabilityType == "XCorr"
    
    cc_zero = CellDataByStim.XcorrData.cc_zero;
    ReliabilitySummary.R1 = mean(cc_zero, 'omitnan');
    ReliabilitySummary.R2 = mean(cc_zero(isRF), 'omitnan');
end
