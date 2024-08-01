function [block] = add_custom_parameters_to_block(block)
% This function will overwrite parameters v1, v2, v3 and stimLength based on
% manually-entered columns in info sheet (e.g. for ephys data)
% 
% Argument(s): 
%   block (struct)
% 
% Returns: 
%   block (struct)
% 
% Notes:
%
%
% TODO:
% Search 'TODO'

%%  Skip this function if v1, v2, v3, or v4 columns are empty or not present in Info
%   These variables are saved in setup when the info sheet is loaded
%   *You must have 'stripfields' set to false when calling tak_read_info_table*

headers = {'v1', 'v2', 'v3', 'stimLength'};

%First check to see if header exists and is not empty
headerExists = zeros(size(headers));
for h = 1:length(headers)
    headerExists(h) = isfield(block.setup, headers{h});
    if headerExists(h)
        if ismissing(block.setup.(headers{h}))
            headerExists(h) = 0;
        end
    end
end

%Quit function early if all not found or empty
if ~any(headerExists)
    return;
end

%% Replace parameters with variables

disp('Replacing parameters with Info variables...');

parameters = struct;
nTrials = length(block.Sound_Time);
vMat = nan(1, nTrials); %Mimic variable format for other data types

%Variable 1
if headerExists(1)
    parameters.variable1 = vMat;
    parameters.variable1(:) = block.setup.v1;
end
    
%Variable 2
if headerExists(2)
    parameters.variable2 = vMat;
    parameters.variable2(:) = block.setup.v2;
end

%Variable 3
if headerExists(3)
    parameters.variable3 = vMat;
    parameters.variable3(:) = block.setup.v3;
end

%stimLength [in milliseconds!]
if headerExists(4)
    parameters.stimLength = vMat;
    parameters.stimLength(:) = block.setup.stimLength;
end

%% Confirm that parameters format is valid
%if not supplied, use default

%We must have V1 and V2
if ~isfield(parameters, 'variable1') && ~isfield(parameters, 'variable2')
    parameters.variable1 = 0;
    parameters.variable2 = 0;
end

%If V1 exists and not V2...
if isfield(parameters, 'variable1') && ~isfield(parameters, 'variable2')
    parameters.variable2 = parameters.variable1;
    parameters.variable2(:) = nan;
end

%If V2 exists and not V1...
if isfield(parameters, 'variable2') && ~isfield(parameters, 'variable1')
    parameters.variable1 = parameters.variable2;
    parameters.variable1(:) = nan;
end

%If no stimLength, make nans
if ~isfield(parameters, 'stimLength')
    parameters.stimLength = vMat;
end

%It's okay for variable3 to not exist

%% Save parameters in block

block.parameters = parameters;

end %end function

