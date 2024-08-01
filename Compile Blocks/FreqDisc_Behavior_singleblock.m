function [block] = FreqDisc_Behavior_singleblock(block)
% DOCUMENTATION IN PROGRESS
% 
% What does this function do?
% This function finds "prep trials" (blocks used to get the animal situated
% with the system each day), and labels those blocks. Additionally, it will
% find the hit and FA rate for a single day. 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes:
%
% TODO: Remove magic numbers 
% Search 'TODO'

%% Only run this function if this is a behavior block (i.e. stim_protocol == 7)

% find prep trials
if ~ismissing(block.setup.stim_protocol)
if block.setup.stim_protocol==7
    findPrep = block.trialType;
    r = ismember(0, findPrep);
    if r ==0
        block.prepTrial = 1;
    else block.prepTrial = 0;
    end
end


%% find hitrate; determine if it is above threshold
if block.setup.stim_protocol==7
    Hits = find(block.Outcome==1);
    misses =find(block.Outcome==0);
    FA = find(block.Outcome==4);
    withold = find(block.Outcome==3);
    block.HitRate = length(Hits)./(length(Hits)+length(misses));
    block.FARate = length(FA)./(length(FA)+length(withold));
end

end
end