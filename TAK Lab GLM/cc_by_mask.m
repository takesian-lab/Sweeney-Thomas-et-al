function [LL,CC] = cc_by_mask (varargin)
% function [LL,CC] = cc_by_mask(maskidx,const,MaskTraces,T)

% GLM_ExDat will calculate the trace associated with a given mask. Whether
% or not those were saved will determine how you run this function:

% 1) if you have already calculated the masktraces:
% 'hasMaskTraces' = 1
% inputs are maskidx - index of traces to calculate
% const = constant paramter
% MaskTraces = model associated with each mask
% T = % Ca data to test

% 2) if you need to calculated the mask traces from scratch:
% T = % Ca data to test
% Xtest + other ops - calculated in XXX
% Ws: filter (const = Ws(1))
% mask
% mask idx?

S = inputParser; 
addParameter(S,'loco_filter', 0); 
addParameter(S,'ch', 1:3); 
addParameter(S,'behaviorMeasures', true); 
addParameter(S,'traceType', 'zscore'); 
addParameter(S,'shadeblips', true); 
addParameter(S,'raw_data', true); 
parse(S, varargin{:});

if hasMaskTraces

LL = 0;
for b = 1:length(maskidx)
    beta = MaskTraces(b,:);
    LL = LL +beta;
end
LL = LL +const;
R = corrcoef(LL,T');
CC = R(1,2);

elseif MakeMaskTraces
    
    
    
end

end

