function type = determine_RF_response_type(types, IsRF)
%% Determine overall RF response type (excitatory, inhibitory, mixed)
%
%  Inputs:
%  - types (string array)
%
%  Outputs:
%  - type (string)
%
% MET August 2022 for use with simple_extract_data

%%

%Optional allow user to filter by IsRF
if nargin < 2
    IsRF = ones(size(types));
else
    %Check that IsRF dimensions match types
    if ~isequal(size(IsRF), size(types))
        error('IsRF dimensions do not match')
    end
end

types_orig = types;
types = types(IsRF == 1);

%Return immediately if no types to check
if isempty(types)
    if any(strcmp(types_orig,'none'))
        type = 'none';
    else
        type = 'undetermined';
    end
    return;
end

allTypes = {'undetermined', 'none', 'excitatory', 'inhibitory'};

%Simplify types to excitatory or inhibitory
types(strcmp(types, 'activated')) = 'excitatory';
types(strcmp(types, 'prolonged')) = 'excitatory';
types(strcmp(types, 'suppressed')) = 'inhibitory';

type = [];
            
%If all of the values in the RF are the same, classify as that type
for t = 1:length(allTypes)
    if all(strcmp(types, allTypes{t}))
        type = allTypes{t};
    end
end

if isempty(type)
    %If that is not the case, remove undetermined and try again
    types(strcmp(types, 'undetermined')) = [];

    for t = 1:length(allTypes)
        if all(strcmp(types, allTypes{t}))
            type = allTypes{t};
        end
    end
end

if isempty(type)
    %If that is still not the case, remove 'none' and try again
    types(strcmp(types, 'none')) = [];

    for t = 1:length(allTypes)
        if all(strcmp(types, allTypes{t}))
            type = allTypes{t};
        end
    end
end

%Now we know that we will have a 'mixed' response type
if isempty(type)
    type = 'mixed';
end

type = convertCharsToStrings(type);