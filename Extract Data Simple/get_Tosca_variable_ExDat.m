% find tosca session in Extracted Data:
function [Sessions] = get_Tosca_variable_ExDat(Blocks)
Sessions = nan(length(Blocks),1);
blocknames = fieldnames(Blocks);
for i = 1:length(blocknames)
    Sessions(i) = Blocks.([blocknames{i}]).setup.Tosca_session;
end
end