function [loco_data, ignoreLoco] = check_loco_data(block, loco_data, redo)
%Check to make sure wheel was functioning during run. If there is no loco
%speed greater than the noise floor of the wheel, the run will get flagged
%and show a figure & pop-up to the user. The user should confirm if this
%was an issue with the wheel or whether the mouse was not running, and 
%verify whether the problem affects just the run or the whole session

%In both cases, a .mat file called IgnoreLoco will be saved in the session folder
%     IgnoreLoco.Session  = 1   ->   User confirms wheel was broken for the entire session
%     IgnoreLoco.Run#     = 1   ->   User confirms wheel was broken for just this run
%     IgnoreLoco.Run#     = 0   ->   User confirms wheel was not broken for just this run

%Upon recompiling, if an IgnoreLoco file already exists, it will be checked
%so that the user does not have to do this process again. Alternatively, the
%user can make and the save the .mat file manually in advance

%If IgnoreLoco == 1, loco_data.speed will be replaced with NaNs

%% Check Loco data

%If redo == 1, show pop-up to user even if a previous IgnoreLoco.mat was found
if nargin < 3
    redo = 0;
elseif nargin == 3
    if ~(redo == 0 || redo == 1)
        error('Unexpected value for redo')
    end
end

ignoreLoco = 0;

% Variables needed from block
mousename = block.setup.mousename;
session = num2str(block.setup.Tosca_session);
run = num2str(block.setup.Tosca_run);

% Look for IgnoreLoco file and load from current directory unless redo == 1
if isfile('IgnoreLoco.mat') && ~redo
    load('IgnoreLoco.mat')
else
    IgnoreLoco = struct;
end

%If IgnoreLoco exists, check fields for this session or run
%Check run first so that any run info will supersede session info
if isfield(IgnoreLoco, ['Run' run])
    if IgnoreLoco.(['Run' run]) == 1
        ignoreLoco = 1;
        disp('Found IgnoreLoco file. Replacing loco with NaNs for this block')
    else
        ignoreLoco = 0;
        disp('Found IgnoreLoco file. Keeping loco data for this block')
    end

elseif isfield(IgnoreLoco, 'Session')
    if IgnoreLoco.Session == 1
        ignoreLoco = 1;
        disp('Found IgnoreLoco file. Replacing loco with NaNs for this block')
    else
        ignoreLoco = 0;
        disp('Found IgnoreLoco file. Keeping loco data for this block')
    end
    
%If IgnoreLoco does not exist, check loco_data trace
elseif ~any(loco_data.speed > block.setup.constant.locoThresh)

    warning('No running detected - possible issue with loco data. Check pop-up dialog and choose whether to replace loco data with NaNs if you think the wheel was not working, or to keep as is if you agree the mouse was not running. DO NOT SKIP THIS STEP. Ask 2P code slack if you need help.')

    figure;
    plot(loco_data.t - loco_data.t(1), loco_data.speed)
    hline(block.setup.constant.locoThresh)
    ylim([0 2])
    ylabel('Running speed')
    xlabel('Seconds')
    title(strcat(mousename, ' Tosca Session', {' '}, session, ' Run', {' '}, run))

    %Structure needed for questdlg function
    question_ops = struct;
    question_ops.Default = 'Loco is fine';
    question_ops.Interpreter = 'none';

    answer = questdlg(strcat('Possible issue with loco data for', {' '}, mousename, ' Tosca Session', {' '}, session, ' Run', {' '}, run, '. Replace with NaNs?'), 'Possible problem with wheel',...
    'Replace loco with NaNs for this run only', 'Replace loco with NaNs for all runs in this session', 'Loco is fine', question_ops);

    %Only 3 boxes allowed for questdlg, so I can't add a 4th option unfortunately

    switch answer
        case 'Replace loco with NaNs for this run only'
            ignoreLoco = 1;
            IgnoreLoco.(['Run' run]) = 1;
            save('IgnoreLoco.mat', 'IgnoreLoco')
        case 'Replace loco with NaNs for all runs in this session'
            ignoreLoco = 1;
            IgnoreLoco.Session = 1; %Don't include session number to make it easy to copy-paste the file
            save('IgnoreLoco.mat', 'IgnoreLoco')
        case 'Loco is fine'
            ignoreLoco = 0;
            IgnoreLoco.(['Run' run]) = 0;
            save('IgnoreLoco.mat', 'IgnoreLoco')
    end
end

%Replace loco data with NaNs
if ignoreLoco
    loco_data.speed(:) = NaN;
end