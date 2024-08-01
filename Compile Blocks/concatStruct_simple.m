function [concatStruct] = concatStruct_simple(inputStruct) 

blockFields = fieldnames(inputStruct); %field names
concatStruct = inputStruct(1); %use the first block as template

for f = 1:length(blockFields)
    thisBlockField = concatStruct.(blockFields{f}); %only one block from the template - just want to check what type of field this is
    
    %% determine if a measure is monotonically increasing (i.e. time measure)
    
    
    %% concatenate blocks based on dimentions
    %if the each block for this field only has one row
    if size(thisBlockField,1)==1 %this includes both arrays that are 1xN and structs 
        %horizontal concat
        concatStruct.(blockFields{f}) = [inputStruct.(blockFields{f})];
    elseif size(size(thisBlockField),2) ==3 ... %if there are 3 dimentions 
            && size(thisBlockField,1) == 3 % and the first dimention is 3 (3 channels)
        concatStruct.(blockFields{f}) = cat(2,inputStruct.(blockFields{f}));
    else % vertical cat
        concatStruct.(blockFields{f}) = cat(1,inputStruct.(blockFields{f}));
    end
    
    clear thisBlockField

end

% now if a given field has all the same value; write it as the first value
for f = 1:length(blockFields)
    thisValue = concatStruct.(blockFields{f});
    if ~isempty(thisValue)
        % check if every value is the same
        if isa(thisValue,'string')
            if all(thisValue(:)==thisValue(1))
                concatStruct.(blockFields{f}) = thisValue(1);
            end
        end
        
        % check for doubles that should be monotonically increasing (i.e.
        % time)
        if isa(thisValue,'double') 
            if length(thisValue(1)) == 1 && all(thisValue(:)==thisValue(1))
                concatStruct.(blockFields{f}) = thisValue(1);
            end 
            
            increaseLogic = diff(inputStruct(1).(blockFields{f}))>0;
            if sum(increaseLogic)> length(inputStruct(1).(blockFields{f}))*0.9
                %if 90% of values in the first block are increasing (i.e. time)
                modifiedValue = [];
                previousEndValue = 0;
                for n = 1:length(inputStruct)
                    currentValues = inputStruct(n).(blockFields{f});
                    dim=find(size(currentValues)~=1);
                    modifiedValue = cat(dim,modifiedValue,(currentValues+previousEndValue));
                    previousEndValue = previousEndValue + currentValues(end);
                    clear currentValues
                end
                concatStruct.(blockFields{f}) = modifiedValue;
                clear modifiedValue
            end
            
        end

    end
end 

end 