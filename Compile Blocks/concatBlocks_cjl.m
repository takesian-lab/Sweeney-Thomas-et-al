%%
% function variables % % % TODO
clear all
clc

compiledBlocksPath = '\\apollo\research\ENT\Takesian Lab\Christine\Data\FiberPhotometry\Data\1_VIPsensor\SxG3K052222M2\CompiledBlocks';
concatBlockPath = '\\apollo\research\ENT\Takesian Lab\Christine\Data\FiberPhotometry\Data\1_VIPsensor\SxG3K052222M2\ConcatBlocksPerDay';
sessionID = '05'; %,'06','07','08','09','10','11','12','13','14','15'};
checkStr = '_TaskTraining';
% function [concatBlocks] = concatBlocks_cjl(compiledBlocksPath,concatBlockPath,sessionID,checkStr)
% to concatenate behavior (and fiber functions) fir Carolyn and Christine
% such that each day of training (for trasktraining) will have only one
% block

%% defile the sessions to be loaded
cd(compiledBlocksPath);
allDirs = dir(['*Session' sessionID '*' checkStr '*']); % % % TODO: convert all the

count = 0;
% loop through all the directories and find the ones that have high hit
% rate, and store all the structs into "filteredBlocks" together
for n=1:length(allDirs)
    thisBlock = load(allDirs(n).name);
    thisBlock = thisBlock.block;
    if thisBlock.HitRate > 0.5 && ~isnan(thisBlock.FARate)
        count = count+1;
        try 
            filteredBlocks(count) = thisBlock;
        catch
            warning('the structure you"re adding is does not have the same fields as the previous ones');
        end 
    end
    
end

%% loop through each field and concatenate based on the dimention sizes 

[concatBlocks_draft] = concatStruct_simple(filteredBlocks);

%% if a field in the concatenated block is a structure itself, go through the same function
blockFields = fieldnames(concatBlocks_draft);

for f = 1:length(blockFields)
    % if there's one layer of nested structure
    if isstruct(concatBlocks_draft.(blockFields{f}))
        
        thisStruct = concatBlocks_draft.(blockFields{f});
        thisConcat = concatStruct_simple(thisStruct);
        
        concatBlocks.(blockFields{f}) = thisConcat;
               
    else %if no changes 
        concatBlocks.(blockFields{f}) = concatBlocks_draft.(blockFields{f});
    end 
end 

%% fix special nested structs
% fibphot
ops_temp = concatStruct_simple(concatBlocks.FibPhot.Ops);
params_temp = concatStruct_simple(ops_temp.params);
concatBlocks.FibPhot.Ops = [] ;
concatBlocks.FibPhot.Ops.params = params_temp;

% setup.constant
constant_temp = concatStruct_simple(concatBlocks.setup.constant);
concatBlocks.setup.constant = constant_temp ;

% block_supname name
runRange = "";
blockRange = "";
for i = 1:length(concatBlocks.setup.block_supname)
    thisName = concatBlocks.setup.block_supname(i);
    if i ==1
        runRange = extractBefore(extractAfter(thisName,'Run'),'Block');
        blockRange = extractAfter(thisName,'Block');
    else    
        runRange = runRange  + ',' +  extractBefore(extractAfter(thisName,'Run'),'Block');
        blockRange = blockRange  + ',' +   extractAfter(thisName,'Block');
    end 
end 

block_supname = extractBefore(concatBlocks.setup.block_supname(1),'-Runs')+...
'-Runs' + runRange + '-Blocks' + blockRange;

concatBlocks.setup.block_supname = block_supname;
concatBlocks.setup.block_filename = 'Concat' + ...
    extractAfter(extractBefore(concatBlocks.setup.block_filename(1),'Run'),'Compiled') + ...
    'Runs' + strrep(runRange,',','-') + '_Blocks' + strrep(blockRange,',','-');

%% save 
if ~exist(concatBlockPath)
    mkdir(concatBlockPath)
end 
save(fullfile(concatBlockPath,[concatBlocks.setup.block_filename + ".mat"]),'concatBlocks');
 