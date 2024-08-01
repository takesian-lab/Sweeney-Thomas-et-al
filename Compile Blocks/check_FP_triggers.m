function new_FibPhot_trial_time = check_FP_triggers(FibPhot_trial_time, block)
%Check that FibPhot triggers match number of block trials and fix if needed

if length(FibPhot_trial_time) == length(block.errorData.Tosca_times)
    %If trial numbers already match, return FibPhot_trial_time
    new_FibPhot_trial_time = FibPhot_trial_time;
    
else
    %Check if fix has already been recorded for this run
    %If it has, a figure will be found in the FibPhot folder

    figtitle = char(block.setup.block_name);
    FPfix_Filename = ['FPfix_' figtitle];
    FPfix2_Filename = ['FPfix2_' figtitle];

    if exist([FPfix_Filename '.fig'], 'file')
        disp('Found FP trigger fix')
        correctFP = 1;
    elseif exist([FPfix2_Filename '.fig'], 'file')
        disp('Found FP trigger fix')
        correctFP = 2;
    else
        warning('Possible issue with FibPhot triggers. Check pop ups')

        if length(FibPhot_trial_time) > length(block.errorData.Tosca_times)
            %If fix hasn't already been recorded, plot figure and ask user for input
            
            %FIX ATTEMPT #1 = REMOVE FIRST TRIGGER
            %----------------
            fig1 = figure;
            subplot(2,1,1); hold on
            stem(diff(FibPhot_trial_time))
            stem(diff(block.errorData.start_time))
            title('Original triggers')
            legend({'FibPhot', 'Tosca'})
    
            subplot(2,1,2); hold on
            stem(diff(FibPhot_trial_time(2:end)))
            stem(diff(block.errorData.start_time))
            title('First FibPhot trigger removed')
            legend({'FibPhot', 'Tosca'})

            tak_suptitle(regexprep(figtitle,'_',' ','emptymatch')) %Replace underscores with space in fig title
            %----------------

            answer = questdlg('Does this fix seem correct?', ...
            'Options', ...
            'Yes','No','Yes');

            switch answer
                case 'Yes'
                    correctFP = 1;

                    %Save figure in folder as record of fix
                    saveas(fig1,FPfix_Filename,'fig');
                    saveas(fig1,FPfix_Filename,'jpg');
                    disp('FP fix saved')
                    close(fig1)

                case 'No'
                    correctFP = 0;
            end
            
            if correctFP == 0
                %FIX ATTEMPT #2 = FIND THE TRIAL WITH THE BIGGEST DIFF TO REMOVE
                %----------------
                
                [~, remove_ind] = max(diff(FibPhot_trial_time));
                
                fig2 = figure;
                subplot(2,1,1); hold on
                stem(diff(FibPhot_trial_time))
                stem(diff(block.errorData.start_time))
                title('Original triggers')
                legend({'FibPhot', 'Tosca'})

                subplot(2,1,2); hold on
                temp = FibPhot_trial_time;
                temp(remove_ind) = [];
                stem(diff(temp))
                stem(diff(block.errorData.start_time))
                title(['FibPhot trigger ' num2str(remove_ind) ' removed'])
                legend({'FibPhot', 'Tosca'})

                tak_suptitle(regexprep(figtitle,'_',' ','emptymatch')) %Replace underscores with space in fig title
                %----------------

                answer = questdlg('Does this fix seem correct?', ...
                'Options', ...
                'Yes','No','Yes');

                switch answer
                    case 'Yes'
                        correctFP = 2;

                        %Save figure in folder as record of fix
                        saveas(fig2,FPfix2_Filename,'fig');
                        saveas(fig2,FPfix2_Filename,'jpg');
                        disp('FP fix saved')
                        close(fig2)

                    case 'No'
                        correctFP = 0;
                end
            end
            
        elseif length(FibPhot_trial_time) < length(block.errorData.Tosca_times)
            figure;
            warning('There is an issue with the FibPHot triggers. Possible reasons are that the block and tosca number are wrong on the info sheet. Another reason could be that a trigger was not sent from Tosca to FibPhot analysis. This will most likely require a manual fix. Use the graph to help figure out what is wrong ')
            stem(diff(FibPhot_trial_time)); hold on
            stem(diff(block.errorData.start_time));
            title(['FibPhot trigger missing. Investigate further...'])
            legend({'FibPhot', 'Tosca'})
            ylabel('difference in time(sec) between trial n and n-1')
            error('It seems like FibPhot ended early. Check Info sheet to make sure block number matches Tosca run. Alternative explanation may be that a TTL was not sent to FP for error trial')
            %If the problem isn't with info sheet, FibPhot may have crashed early
            %Choose whether you prefer to delete the mismatched trials from Tosca data or not use the FibPhot data
            %Example fix: If you have 211 triggers and 215 trials, first plot the figure above to confirm
            %the first 211 triggers look normal, back everything up, then delete the .di, .u8, and .avi files
            %for trials 212:215, delete their record in the Run.txt file, and delete their record in the Run.avi file
        end
    end
    
    if correctFP == 1
        new_FibPhot_trial_time = FibPhot_trial_time(2:end);
    elseif correctFP == 2
        [~, remove_ind] = max(diff(FibPhot_trial_time));
        new_FibPhot_trial_time = FibPhot_trial_time;
        new_FibPhot_trial_time(remove_ind) = [];
    else
        error('Number of FibPhot triggers does not match Tosca trials')
    end
end