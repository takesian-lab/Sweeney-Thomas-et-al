function [block, Bruker_trial_time] = check_Bruker_triggers(Bruker_trial_time, New_sound_times, block)
%Check that Bruker trials match number of Tosca trials and fix if needed
%FYI: use visualize_VoltageRecording to plot raw voltage recording

%% Set up variables and return if trials already match

Sound_Time = get_Sound_Time_behavior(block, 1);
brukerTrialIndex = 1:length(Bruker_trial_time);
toscaTrialIndex = 1:length(Sound_Time);
errorTrials = block.errorData.error_trials;

%If trial numbers already match, return
if length(Bruker_trial_time) == length(Sound_Time)
    if length(Sound_Time) == length(New_sound_times)
        return
    else
        error('Sound_Time and New_sound_times do not match')
    end
end
    
%% Determine how to fix data
%Check if fix has already been recorded for this run
%If it has, a figure will be found in the Tosca folder
%If fix hasn't already been recorded, plot figures and ask user for input

previous_dir = pwd; %Record previous directory so we can send user back there at the end of the script
cd(block.setup.Tosca_path); %Go to Tosca folder

figtitle = char(block.setup.block_supname);
Brukerfix_Filename = ['Brukerfix_' figtitle];
Brukerfix2_Filename = ['Brukerfix2_' figtitle];
Brukerfix4_Filename = ['Brukerfix4_' figtitle];
Brukerfix5_Filename = ['Brukerfix5_' figtitle];

if exist([Brukerfix_Filename '.fig'], 'file')
    disp('Found Bruker trigger fix')
    correctBruker = 1;
elseif exist([Brukerfix2_Filename '.fig'], 'file')
    disp('Found Bruker trigger fix')
    correctBruker = 2;
elseif exist([Brukerfix4_Filename '.fig'], 'file')
    disp('Found Bruker trigger fix')
    correctBruker = 4;
elseif exist([Brukerfix5_Filename '.fig'], 'file')
    disp('Found Bruker trigger fix')
    correctBruker = 5;
else
    warning('Possible issue with Bruker triggers. Check pop ups')
    
    if length(Sound_Time) > length(Bruker_trial_time)
        %%%%%%%%%%% MORE TOSCA TRIALS THAN TRIGGERS FOUND IN BRUKER VR %%%%%%%%%%% 
        %For some reason, we have more Tosca trials than Bruker received triggers for
        %Even if sounds were actually played here, we cannot align them to any time in
        %the 2P data, so we need to remove those trials. Figure out which trials to remove here
    
        %FIX ATTEMPT #1 = REMOVE FIRST N TOSCA TRIGGERS
        %----------------
        N = length(Sound_Time) - length(Bruker_trial_time);

        fig1 = figure;
        subplot(2,1,1); hold on
        stem(diff(Bruker_trial_time))
        stem(diff(Sound_Time))
        if ~isempty(errorTrials)
            for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
            title('Original trigger ISIs - OK for data under green lines to differ')
        else
            title('Original trigger ISIs')
        end
        ylabel('ISI')
        legend({['Bruker N=' num2str(length(Bruker_trial_time))], ['Tosca N=' num2str(length(Sound_Time))]})

        subplot(2,1,2); hold on
        stem(diff(Bruker_trial_time))
        stem(diff(Sound_Time(1+N:end)))
        if ~isempty(errorTrials)
            for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
            title(['First ' num2str(N) ' Tosca trial(s) removed - OK for data under green lines to differ'])
        else
            title(['First ' num2str(N) ' Tosca trial(s) removed'])
        end
        ylabel('ISI')
        legend({['Bruker N=' num2str(length(Bruker_trial_time))], ['Tosca N=' num2str(length(Sound_Time)- N)]})
        tak_suptitle(regexprep(figtitle,'_',' ','emptymatch')) %Replace underscores with space in fig title
        
        answer = questdlg('Does this fix seem correct?', ...
        'Options', ...
        'Yes','No','Yes');

        switch answer
            case 'Yes'
                correctBruker = 1;
                %Save figure in folder as record of fix
                saveas(fig1,Brukerfix_Filename,'fig');
                saveas(fig1,Brukerfix_Filename,'jpg');
                disp('Bruker fix saved')
                close(fig1)

            case 'No'
                correctBruker = 0;
        end

        %FIX ATTEMPT #2 = FIND THE TRIAL WITH THE BIGGEST DIFF TO REMOVE
        %----------------
        if correctBruker == 0
            
            [~, remove_ind] = max(diff(Sound_Time));

            fig2 = figure;
            subplot(2,1,1); hold on
            stem(diff(Bruker_trial_time))
            stem(diff(Sound_Time))
            if ~isempty(errorTrials)
                for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
                title('Original trigger ISIs - OK for data under green lines to differ')
            else
                title('Original trigger ISIs')
            end
            ylabel('ISI')
            legend({['Bruker N=' num2str(length(Bruker_trial_time))], ['Tosca N=' num2str(length(Sound_Time))]})

            subplot(2,1,2); hold on
            stem(diff(Bruker_trial_time))
            temp = Sound_Time;
            temp(remove_ind) = [];
            stem(diff(temp))
            if ~isempty(errorTrials)
                for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
                title(['Tosca trial ' num2str(remove_ind) ' removed - OK for data under green lines to differ'])
            else
                title(['Tosca trial ' num2str(remove_ind) ' removed'])
            end
            ylabel('ISI')
            legend({['Bruker N=' num2str(length(Bruker_trial_time))], ['Tosca N=' num2str(length(Sound_Time)-length(remove_ind))]})
            tak_suptitle(regexprep(figtitle,'_',' ','emptymatch')) %Replace underscores with space in fig title

            answer = questdlg('Does this fix seem correct?', ...
            'Options', ...
            'Yes','No','Yes');

            switch answer
                case 'Yes'
                    correctBruker = 2;

                    %Save figure in folder as record of fix
                    saveas(fig2,Brukerfix2_Filename,'fig');
                    saveas(fig2,Brukerfix2_Filename,'jpg');
                    disp('Bruker fix saved')
                    close(fig2)

                case 'No'
                    correctBruker = 0;
            end
        end

        %FIX ATTEMPT #3 = REMOVE LAST N TOSCA TRIGGERS
        %----------------
        if correctBruker == 0
            
            N = length(Sound_Time) - length(Bruker_trial_time);

            fig3 = figure;
            subplot(2,1,1); hold on
            stem(diff(Bruker_trial_time))
            stem(diff(Sound_Time))
            if ~isempty(errorTrials)
                for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
                title('Original trigger ISIs - OK for data under green lines to differ')
            else
                title('Original trigger ISIs')
            end
            ylabel('ISI')
            legend({['Bruker N=' num2str(length(Bruker_trial_time))], ['Tosca N=' num2str(length(Sound_Time))]})

            subplot(2,1,2); hold on
            stem(diff(Bruker_trial_time))
            stem(diff(Sound_Time(1:end-N)))
            if ~isempty(errorTrials)
                for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
                title(['Last ' num2str(N) ' Tosca trial(s) removed - OK for data under green lines to differ'])
            else
                title(['Last ' num2str(N) ' Tosca trial(s) removed'])
            end
            ylabel('ISI')
            legend({['Bruker N=' num2str(length(Bruker_trial_time))], ['Tosca N=' num2str(length(Sound_Time)-N)]})
            tak_suptitle(regexprep(figtitle,'_',' ','emptymatch')) %Replace underscores with space in fig title

            answer = questdlg('Does this fix seem correct?', ...
            'Options', ...
            'Yes','No','Yes');

            switch answer
                case 'Yes'
                    %In this scenario, it appears that the Bruker data stopped running before the Tosca data
                    %(maybe Prairie View crashed or there was not enough time in the Voltage Recording)
                    %Here, the solution is to manually delete the excess Tosca trials from the Tosca Session folder
                    % 1. Backup session folder
                    % 2. Open Mousename-Session#-Run#.txt and scroll to the bottom. Delete the information for the last N trials
                    % 3. Delete the Mousename-Session#-Run#-Trial##.di file
                    %Next time you compile the block, Ken's Tosca scripts
                    %will automatically correct the loco data(and hopefully pupil data too)
                    %You MUST do this manually (no automatic fix) because other places in the code will be trying to remove error trials > length misc. variables
                    error(['Block has ' num2str(N) ' too many trials at the end. Fix by deleting last ' num2str(N) ' trial text files and info from Tosca Run file (backup first).']) 

                case 'No'
                    correctBruker = 0;
            end
        end

    elseif length(Sound_Time) < length(Bruker_trial_time)
        %%%%%%%%%%% FEWER TOSCA TRIALS THAN TRIGGERS FOUND IN BRUKER VR %%%%%%%%%%% 
        %This time, we have more triggers in the Bruker voltage recording than number
        %of trials recorded by Tosca.
        
        %FIX ATTEMPT #4 = REMOVE FIRST N BRUKER TRIGGERS
        %----------------
        N = length(Bruker_trial_time) - length(Sound_Time);

        fig4 = figure;
        subplot(2,1,1); hold on
        stem(diff(Bruker_trial_time))
        stem(diff(Sound_Time))
        if ~isempty(errorTrials)
            for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
            title('Original triggers ISIs - OK for data under green lines to differ')
        else
            title('Original trigger ISIs')
        end
        ylabel('ISI')
        legend({['Bruker N=' num2str(length(Bruker_trial_time))], ['Tosca N=' num2str(length(Sound_Time))]})

        subplot(2,1,2); hold on
        stem(diff(Bruker_trial_time(1+N:end)))
        stem(diff(Sound_Time))
        if ~isempty(errorTrials)
            for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
            title(['First ' num2str(N) ' Bruker trigger(s) removed - OK for data under green lines to differ'])
        else
            title(['First ' num2str(N) ' Bruker trigger(s) removed'])
        end
        ylabel('ISI')
        legend({['Bruker N=' num2str(length(Bruker_trial_time)-N)], ['Tosca N=' num2str(length(Sound_Time))]})
        tak_suptitle(regexprep(figtitle,'_',' ','emptymatch')) %Replace underscores with space in fig title
        
        answer = questdlg('Does this fix seem correct?', ...
        'Options', ...
        'Yes','No','Yes');

        switch answer
            case 'Yes'
                correctBruker = 4;
                %Save figure in folder as record of fix
                saveas(fig4,Brukerfix4_Filename,'fig');
                saveas(fig4,Brukerfix4_Filename,'jpg');
                disp('Bruker fix saved')
                close(fig4)

            case 'No'
                correctBruker = 0;
        end
        
        %FIX ATTEMPT #5 = REMOVE LAST N BRUKER TRIGGERS
        %----------------
        if correctBruker == 0
            
            N = length(Bruker_trial_time) - length(Sound_Time);

            fig5 = figure;
            subplot(2,1,1); hold on
            stem(diff(Bruker_trial_time))
            stem(diff(Sound_Time))
            if ~isempty(errorTrials)
                for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
                title('Original trigger ISIs - OK for data under green lines to differ')
            else
                title('Original trigger ISIs')
            end
            ylabel('ISI')
            legend({['Bruker N=' num2str(length(Bruker_trial_time))], ['Tosca N=' num2str(length(Sound_Time))]})

            subplot(2,1,2); hold on
            stem(diff(Bruker_trial_time(1:end-N)))
            stem(diff(Sound_Time))
            if ~isempty(errorTrials)
                for i = 1:length(errorTrials); vline(errorTrials(i)-1,'g'); end
                title(['Last ' num2str(N) ' Bruker trigger(s) removed - OK for data under green lines to differ'])
            else
                title(['Last ' num2str(N) ' Bruker trigger(s) removed'])
            end
            ylabel('ISI')
            legend({['Bruker N=' num2str(length(Bruker_trial_time)-N)], ['Tosca N=' num2str(length(Sound_Time))]})
            tak_suptitle(regexprep(figtitle,'_',' ','emptymatch')) %Replace underscores with space in fig title

            answer = questdlg('Does this fix seem correct?', ...
            'Options', ...
            'Yes','No','Yes');

            switch answer
                case 'Yes'
                    correctBruker = 5;
                    %Save figure in folder as record of fix
                    saveas(fig5,Brukerfix5_Filename,'fig');
                    saveas(fig5,Brukerfix5_Filename,'jpg');
                    disp('Bruker fix saved')
                    close(fig5)
                case 'No'
                    correctBruker = 0;
            end
        end
        
        %IF NOTHING WORKED
        if correctBruker == 0
            %The first thing to do is doublecheck that everything is correct in your info sheet.
            %If this isn't the problem, plotting both traces can give you insight about what went wrong.
            error('Bruker voltage recording has more trials than the Tosca data. Check Info sheet to make sure block number matches Tosca run.')

            %Plot all state and trial triggers
            visualize_VoltageRecording(block);
            
            %Plot ISIs to troubleshoot
            figure; hold on
            subplot(3,1,1)
            stem(diff(Bruker_trial_time))
            ylabel('ISI')
            title(['Bruker trigger ISIs - N trials = ' num2str(length(Bruker_trial_time))])

            subplot(3,1,2)
            stem(diff(Sound_Time),'r')
            ylabel('ISI')
            title(['Tosca trial ISIs - N trials = ' num2str(length(Sound_Time))])

            subplot(3,1,3); hold on
            stem(diff(Bruker_trial_time))
            stem(diff(Sound_Time),'r')
            ylabel('ISI')
            title('Both')
            legend({'Bruker', 'Tosca'})
        end
    end
end
        
%% Fix data based on determined correction

if correctBruker == 1
    N = length(Sound_Time) - length(Bruker_trial_time);
    remove_ind = 1:N;
elseif correctBruker == 2
    [~, remove_ind] = max(diff(Sound_Time));
elseif correctBruker == 4
    N = length(Bruker_trial_time) - length(Sound_Time);
    remove_ind = 1:N;
elseif correctBruker == 5
    N = length(Bruker_trial_time) - length(Sound_Time);
    remove_ind = length(Bruker_trial_time)-(N-1):length(Bruker_trial_time);
else
    %No automatic correction was found
    visualize_VoltageRecording(block)
    error('Number of Bruker triggers does not match Tosca trials')
end

%If a correction was found, edit the block by removing the extra trials
if correctBruker > 0
    
    if length(Sound_Time) > length(Bruker_trial_time)
        %More tosca trials than Bruker trials
        
        %Record which Bruker trial indices correspond to which Tosca
        %trial indices after block has been fixed
        toscaTrialIndex(remove_ind) = [];
        block.errorData.corrected_BrukerTrialIndex = brukerTrialIndex;
        block.errorData.corrected_ToscaTrialIndex = toscaTrialIndex;
        
        %Find out if the extra trials were already recorded as error trials
        alreadyError = zeros(size(remove_ind));
        for i = 1:length(remove_ind)
            if ismember(errorTrials,remove_ind(i))
                alreadyError(i) = 1;
            end
        end
        
        %If so, we don't need to remove them from the block
        remove_ind(alreadyError == 1) = [];
        
        %Catch situation we haven't run into yet
        if ~isempty(errorTrials)
            if ~all(remove_ind < min(errorTrials)) 
                %If any of the trials to be removed come after the error
                %trials, then the indexing for which trials to remove from
                %the main block could be messed up (since error trials were
                %already removed from there) -> Code for this when we encounter it!
                error('Code needs to be added here')
            end
        end
        
        %Add removed trials to list of error trials
        block.errorData.error_trials = unique([remove_ind, errorTrials]); 
        
        %Remove extra trial(s) from block
        block.start_time(remove_ind) = [];
        block.Tosca_times(remove_ind) = [];
        block.New_sound_times(remove_ind) = [];
        block.New_sound_idx(remove_ind) = [];
        block.lick_time(remove_ind) = [];
        block.Outcome(remove_ind) = [];
        block.trialType(remove_ind) = [];
        block.rxn_time(remove_ind) = [];
        block.water_delivery(remove_ind) = [];
        if isfield(block,'loc_Trial_times')
            block.loc_Trial_times(remove_ind) = [];
            block.loc_Trial_activity(remove_ind) = [];
        end
        paramfields = fields(block.parameters);
        for f = 1:length(paramfields)
            if ~isempty(block.parameters.(paramfields{f}))
                block.parameters.(paramfields{f})(remove_ind) = [];
            end
        end

    elseif length(Sound_Time) < length(Bruker_trial_time)
        %More Bruker trials than Tosca trials
        %Edit Bruker_trial_time itself
        Bruker_trial_time(remove_ind) = [];
        
        %Record which Bruker trial indices correspond to which Tosca
        %trial indices after block has been fixed
        brukerTrialIndex(remove_ind) = [];
        block.errorData.corrected_BrukerTrialIndex = brukerTrialIndex;
        block.errorData.corrected_ToscaTrialIndex = toscaTrialIndex;
    end
end

%Go back to previous directory
cd(previous_dir)

end