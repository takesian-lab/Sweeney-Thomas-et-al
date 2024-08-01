
function result = check_complete_triplet(V1,result)
for ii = 1:length(V1)

    % start with error trials and fix the triplet if one of them is broken:
    if result(ii) =='Error'
        if V1(ii) =='A'
            if ii ~= length(V1)
                if V1(ii+1) == 'B'; result(ii+1) = 'Error'; end
                if V1(ii+2) == 'C1' | V1(ii+2) == 'C2'; result(ii+2) = 'Error'; end
            end
        elseif V1(ii) =='B'
            if V1(ii-1) == 'A'; result(ii-1) = 'Error'; end
             if ii ~= length(V1)
            if V1(ii+1) == 'C1' | V1(ii+1) == 'C2'; result(ii+1) = 'Error'; end
             end
        elseif V1(ii) == 'C1' | V1(ii) == 'C2'
            if V1(ii-2) == 'A'; result(ii-2) = 'Error'; end
            if V1(ii-1) == 'B'; result(ii-1) = 'Error'; end
        elseif V1(ii) == 'X'
                 if ii ~= length(V1)
            if V1(ii+1) == 'Y'; result(ii+1) = 'Error'; end
            if V1(ii+2) == 'Z1' | V1(ii+2) == 'Z2'; result(ii+2) = 'Error'; end
                 end
        elseif V1(ii) == 'Y'
            if V1(ii-1) == 'X'; result(ii-1) = 'Error'; end
                 if ii ~= length(V1)
            if V1(ii+1) == 'Z1' | V1(ii+1) == 'Z2'; result(ii+1) = 'Error'; end
                 end
        elseif V1(ii) == 'Z1' | V1(ii) == 'Z2'
            if V1(ii-2) == 'X'; result(ii-2) = 'Error'; end
            if V1(ii-1) == 'Y'; result(ii-1) = 'Error'; end
        else
            error('error in triplet')
        end
    end
    
    
    
    % double check that the triplets went in order (no missing trials):
    if V1(ii) == 'A'
        if ii == length(V1) %doesnt end in complete triplet
            result(ii) = 'Error';
        elseif V1(ii+1) ~= 'B' % not followed by correct value
            result(ii) = 'Error';
        end
    end
    
    if V1(ii) == 'B'
        if ii == length(V1) %doesnt end in complete triplet
            result(ii) = 'Error';
            result(ii-1) = 'Error'
        elseif V1(ii +1) ~='C1' & V1(ii +1) ~='C2'
            result(ii) = 'Error'
            if V1(ii-1) == 'A'
                result(ii-1) = 'Error'
            end
        elseif V1(ii -1) ~='A'
            result(ii) = 'Error';
        end
    end
    
    if V1(ii) == 'C1' | V1(ii) == 'C2'
        if V1(ii-1) ~= 'B'
            result(ii) = 'Error';
        end
    end
    
    
    if V1(ii) == 'X'
        if ii == length(V1)
            result(ii) = 'Error';
        elseif V1(ii+1) ~= 'Y'
            result(ii) = 'Error';
        end
    end
    
    if V1(ii) == 'Y'
        if ii == length(V1)
            result(ii) = 'Error';
            result(ii-1) = 'Error';
        elseif V1(ii +1) ~='Z1' & V1(ii +1) ~='Z2'
            result(ii) = 'Error'
        elseif V1(ii -1) ~='X'
            result(ii) = 'Error';
        end
    end
    
    if V1(ii) == 'Z1' | V1(ii) == 'Z2'
        if V1(ii-1) ~= 'Y'
            result(ii) = 'Error';
        end
    end
end
end
