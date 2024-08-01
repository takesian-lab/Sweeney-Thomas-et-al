% updating blocks for Pupillometry and DLC data

block_path = '\\apollo\research\ENT\Takesian Lab\Maryse\FibPhot analysis\Troubleshooting\Blip blocks';
save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\FibPhot analysis\Troubleshooting\Fixed blocks';
recompile = 0;
compute_whiskerpad = 1;
redo_fibphot = 0; %To redo fibphot analysis with V3
EFR = 30; %for behavior experiments where FP or S2p failed/werent used. Choose 20 if matching FP data, 30 for matching S2P data ('estimated frame rate')

% setup for Fiber Photometry analysis
Ops.FibPhot.plot_spectrogram = 0; % plot spectrograms? no = 0, yes = 1
Ops.FibPhot.plot_graphs = 0; % plot graphs? no = 0, yes = 1
Ops.FibPhot.params.remove_blips = 1; 
% if 0, will error on data interuptions, 
% if 1, will interpolate data interuptions and give error under block.FibPhot.error
% if 2, will inject NaNs during data interuptions and give error under block.FibPhot.error
Ops.FibPhot.params.end_analysis = 1; % for Fiber Photometry analyis, 0 = continue with error, 1 = abort analysis if spectrogram fails on all channels (eg. no LED on) 

Ops.FibPhot.params.freqStep   = 1; % Step size in Hz for spectrogram
Ops.FibPhot.params.freqStepWidth = 7; % Frequency band around peak for spectrogram
Ops.FibPhot.params.inclFreqWin = 4; % Number of frequency bins to average for signal (on either side of peak freq) 
Ops.FibPhot.params.detrendWindowTime = 60; % in seconds
Ops.FibPhot.params.lowPassCorner = 100; % low pass filter - keep at 100 for no filter
Ops.FibPhot.params.notes = []; % establish notes

% final downsampled freqeuncy, if = [], then program will define it
% automatically based on freqeuncies of LEDs 
Ops.FibPhot.params.finalSampleFreq = 30; % []; % this value can't be above 37Hz!
               
%these values are based on acquisition rates and AIN number
Ops.FibPhot.params.numChannels = 10; % number of AIN channels, defined by hardware 
Ops.FibPhot.params.scanRate = 2000; % acquisition scan rate, determined by program

%%
cd(block_path)
allblocks=dir('*Compiled*');


for i = 1:length(allblocks)
     filename = allblocks(i).name;
     disp(filename);
     if ~recompile
        cd(save_path)
        if isfile(filename)
            disp('Skipping (already compiled)');
            continue
        end
    end
     
     
     
    cd(block_path)
    load(filename);

    %IF YOU WANT TO REDO FIBPHOT ONLY, UNCOMMENT THIS SECTION AND COMMENT THE REST
%     [block] = FiberPhotometryAnalysis_blocks(block,Ops.FibPhot);
%     if ~any(contains(block.FibPhot.Ops.params.notes,'No fiber photometry data to analyze!!'))
%         [block] = define_sound_fibPhot(block);
%         [block] = align_to_stim_FibPhot(block);
%     end
%         
%     cd(save_path)
%    save(strcat(save_path, '\', filename), 'block');

    if isfield(block,'PupilFrameData') % we are only updating blocks with video data....
        block = define_behavior_singleblock(block);
        %block = generate_pupil_timestamp(block); %MET: I don't think we need this
        %because I think it runs in define_pupillometry and/or define_whiskerpad
        block = define_pupillometry(block);
        if compute_whiskerpad
            [block] = define_whiskerpad(block);
        end

        %Redo align_to_stim
        if strcmp(block.setup.analysis_name, 'FibPhot')
            if redo_fibphot
                end_analysis = 1;
                [block] = FiberPhotometryAnalysis_blocks(block,Ops.FibPhot);
            end

            if ~any(contains(block.FibPhot.params.notes,'No fiber photometry data to analyze!!'))
                [block] = define_sound_fibPhot(block);
                [block] = align_to_stim_FibPhot(block);
            else % generate stim-aligned variables if no FP
                [block] = align_to_stim_behavior(block, block.setup.framerate);
            end

        elseif ismissing(block.setup.block_path) && ismissing(block.setup.VR_path) % not FP and not S2P
            [block] = align_to_stim_behavior(block,EFR);
        else
            %find the stim-aligned traces for S2P data
            [block] = align_to_stim(block);
        end

        cd(save_path)
        save(strcat(save_path, '\', filename), 'block');

    end

end
disp('done correcting blocks');