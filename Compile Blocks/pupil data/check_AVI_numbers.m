
function newAVI = check_AVI_numbers(TL, AVI)
% correct when AVI indicies are greater than actual AVI files. 
% written by Carolyn Sweeney Aug 2022 to work with Ken Hancock's tosca_merge_avi_log

correctAVI = 0;

%Find the number of mismatched trials
extra = length(AVI) - length(TL.trials);

%Get all AVI times
AVItime = [AVI.toscaTime]';

%Get all Tosca times
starttime = nan(length(TL.trials),1);
for i = 1:length(TL.trials)
    starttime(i) = TL.trials{1,i}.start;
end

%Find duration of each trial (goal is that these should match for both AVI and Tosca times)
dAVI = diff(AVItime);
dstart = diff(starttime);

%Attempt two ways to fix the mismatch

%FIX ATTEMPT #1: Remove the first N extra trials
%-------------------
AVItime_fix = AVItime;
AVItime_fix(1:extra) = [];
dAVI_fix = diff(AVItime_fix);

%FIX ATTEMPT #2: Remove N trials with the biggest trial differences
%-------------------
[~, remove_ind] = maxk(dAVI,extra);
AVItime_fix2 = AVItime;
AVItime_fix2(remove_ind) = [];
dAVI_fix2 = diff(AVItime_fix2);

%Check if fix has already been recorded
%If it has, a figure will be found in the Tosca folder
figtitle = TL.filename;
figtitle = figtitle(1:end-4); %Remove '.txt'

AVIfix_Filename = ['AVIfix_' figtitle];
AVIfix2_Filename =  ['AVIfix2_' figtitle];

if exist([AVIfix_Filename '.fig'], 'file')
    disp('Found AVI fix')
    correctAVI = 1;
elseif exist([AVIfix2_Filename '.fig'], 'file')
    correctAVI = 2;
else
    warning('Possible issue with avi labels. Check pop ups')
    %If fix hasn't already been recorded, plot figure and ask user for input
    
    %FIX 1
    %----------------
    fig1 = figure;
    tiledlayout(2,1)
    nexttile
    stem(dAVI); hold on
    stem(dstart)
    title('original AVI times with tosca times')
    legend({'AVI', 'Tosca'})
    
    nexttile
    stem(dAVI_fix); hold on
    stem(dstart)
    title('corrected AVI times with first trial(s) removed')
    legend({'AVI', 'Tosca'})

    tak_suptitle(figtitle)
    %----------------

    answer = questdlg('Does this fix seem correct?', ...
    'Options', ...
    'Yes','No','Yes');

    switch answer
        case 'Yes'
            correctAVI = 1;
            
            %Save figure in folder as record of fix
            saveas(fig1,AVIfix_Filename,'fig');
            saveas(fig1,AVIfix_Filename,'jpg');
            disp('AVI fix saved')
            close(fig1)

        case 'No'
            correctAVI = 0;
    end
    
    if correctAVI == 0
        %FIX 2
        %----------------
        fig2 = figure;
        tiledlayout(2,1)
        nexttile
        stem(dAVI); hold on
        stem(dstart)
        title('original AVI times with tosca times')
        legend({'AVI', 'Tosca'})

        nexttile
        stem(dAVI_fix2); hold on
        stem(dstart)
        title('corrected AVI times with greatest diff trial(s) removed')
        legend({'AVI', 'Tosca'})

        tak_suptitle(figtitle)
        %----------------

        answer = questdlg('Does this fix seem correct?', ...
        'Options', ...
        'Yes','No','Yes');

        switch answer
            case 'Yes'
                correctAVI = 2;

                %Save figure in folder as record of fix
                saveas(fig2,AVIfix2_Filename,'fig');
                saveas(fig2,AVIfix2_Filename,'jpg');
                disp('AVI fix saved')
                close(fig2)

            case 'No'
                correctAVI = 0;
        end
    end
end

%---------------
if correctAVI == 1
    count = extra;
    for i = 1:length(AVI)-extra
        count = count + 1;
        newAVI(i).aviFile = AVI(count).aviFile;
        newAVI(i).toscaTime = AVI(count).toscaTime;
        newAVI(i).states = AVI(count).states;
        newAVI(i).frameRate = AVI(count).frameRate;
    end
elseif correctAVI == 2
    count = 0;
    for i = 1:length(AVI)
        if any(ismember(remove_ind, i))
            continue
        end
        count = count + 1;
        newAVI(count).aviFile = AVI(i).aviFile;
        newAVI(count).toscaTime = AVI(i).toscaTime;
        newAVI(count).states = AVI(i).states;
        newAVI(count).frameRate = AVI(i).frameRate;
    end
else 
    error('Issue with AVI labels')
end
    
end

