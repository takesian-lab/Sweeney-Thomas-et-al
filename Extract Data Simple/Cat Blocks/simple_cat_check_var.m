% when concatenating blocks for ExtractData, sometimes a variable is
% missing from one block but not others. Here, we cat the variables,
% check for missing data, and pad the matrix appropriately.


% this function works for missing pupil,whisker, or velocity. These
% parameters might be missing and it will pad NaN in the missing data.

% TODO: There might be a need to add missing functional data in the future.
%
function [dataout] = simple_cat_check_var(data1, data2, datafield)


% data1 : aligned to stim data
% data2: data to add to data1
if isfield(data2,datafield) & ~isempty(data1.([datafield])) % add the aligned data to the existing concat matrix
    dataout    = [data1.([datafield]);data2.([datafield])];
elseif ~isempty(data1.([datafield]))  %if datafield is missing from the block but has been added to the cat block already
    if isfield(data2,'velocity') % most likely missing whisker/pupil because of camera crashes. Velocity is likely to be there
        dataout = [data1.([datafield]);nan(size(data2.velocity))];
        warning('using size of aligned velocity to pad matrix')
    elseif isfield(data2,'whisker')
        dataout = [data1.([datafield]);nan(size(data2.whisker))];
        warning('using size of aligned whisker to pad matrix')
    elseif isfield(data2,'licks')
        dataout = [data1.([datafield]);nan(size(data2.licks))];
        warning('using size of aligned licks to pad matrix')
    else
        error('Missing variable for concat')
    end
elseif isempty(data1.([datafield])) & isfield(data2,([datafield])) % data are in this block but not previous ones
    if ~isempty(data1.F_stim) %Fibphot data should have this
        %set size of matrix to pad with nan
        size1 = size(data1.F_stim,2);
        size2 = size(data1.F_stim,3);
    elseif ~isempty(data1.F7_stim) %S2P data will have this
        size1 = size(data1.F7_stim,1);
        size2 = size(data1.F7_stim,2);
    elseif ~isempty(data1.licks) % no functional data might point to a purely behavioral experiment which might have licks
        size1 = size(data1.licks,1);
        size2 = size(data1.licks,2);
    else
        % TODO: add new data to blocks if it was missing, but there are no
        % functional data either.
        error('data exist in current block, but not previous. Error in concat. also, no functional data')
    end
    
    size3 = size(data2.([datafield]),1);
    % the functional data is concatenated before this step, so it is bigger
    % by 1 block in dim2. Fix that here:
    newsize1 = size1-size3;
   
    nanmat = nan(newsize1,size2);
    dataout = [nanmat;data2.([datafield])];
else
    dataout = data1.([datafield]);
end
end

