% create cell list from ridge regression data
% input RD = ridge data from GLM_ExDat
%       stimtype: which type of stim was modeled ('FM','RF',etc)
function [celldata] = cellList_from_ridge(RD,stimType)

celldata = table;
for rr = 1:height(RD);
    tHeight = length(RD.params(rr).cellnum);
    cc = table;
    cc.BlockName = repelem(RD.BlockName(rr),tHeight)';
    cc.Group = repelem(RD.Group(rr),tHeight)';
    cc.MouseID = repelem(RD.MouseID(rr),tHeight)';
    cc.StimType = repelem(string(stimType),tHeight)';
    cc.cellnum = RD.params(rr).cellnum;
    cc.cellMatch = RD.params(rr).cellMatch;
    cc.cellIDX = (1:tHeight)';
    celldata = [celldata; cc];
end
end