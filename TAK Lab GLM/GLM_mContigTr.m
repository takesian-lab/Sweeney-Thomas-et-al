
% Set up test and training trials as one, long contiguous segrment rather
% than many individual trials

function [tIDX1, tIDX2] = GLM_mContigTr(t,Per1)

C = size(t,1);
endbound = C-(C*Per1-1); % end boundary for selecting point
pt = randsample(1:endbound,1);

X1 = pt:pt+(Per1*C);
 tIDX1=false(length(t),1);
 tIDX1(X1) = true;
 
X2 = [1:pt-1,(tIDX1(end)+1:C)];
 tIDX2=false(length(t),1);
 tIDX2(X2) = true;
end




    
   