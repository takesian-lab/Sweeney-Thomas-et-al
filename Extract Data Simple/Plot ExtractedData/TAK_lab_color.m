function [outputcolor] = TAK_lab_color(ColorsToUse) 
%TAK Lab color schemes for projects, compatible with
% plot_simple_ExtractedData

ColorsToUse = convertCharsToStrings(ColorsToUse);
% Red (NDNF) / Blue (VIP), 3 activity types - Carolyn's posters
if ColorsToUse == 'Reds'
    outputcolor = [ 234 34 39; 242 118 124; 250 210 211];
    
elseif ColorsToUse == 'Blues'
    outputcolor =  [ 47 54 143; 129 130 188; 170 172 210];
elseif ColorsToUse == 'Purples'
    outputcolor =  [ 128 0 128; 171 0 171; 205 102 205];
elseif ColorsToUse == 'Greens'
    outputcolor =  [43 86 0; 85 171 0; 153 205 102];
else
    error ('color not specified, try using colorbrewer instead')
end



 outputcolor = outputcolor./256; % used illustrator colors - convert to make compatible with matlab
end




