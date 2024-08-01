function [V1, V2, V3, stimLength] = define_stim_parameters(Data, Params, setup)
% Get stim parameters from Tosca
% Old version of this script is commented at the very bottom

% these are arbitrary variable placeholders that vary protocol-to-protocol 
% in theory, there could be more than 3... but in practice 3 is sufficient
V1 = [];
V2 = []; 
V3 = [];
stimLength = []; % in ms

%% 1 = Noiseburst (retired)

if setup.stim_protocol == 1
   
    print('We have phased out stim_protocol = 1 in favor of stim_protocol = 10 for Noiseburst. Compiling will pursue normally')
    setup.stim_protocol = 10;
%     [V1, V2, stimLength] = deal(nan(1,length(Data)));
    
    %Stim Length
%     for m = 1:length(Data)
%         if isfield(Params, 'Output_States')
%             stimLength(1,m) = Params.Output_States(2).StimChans.Burst.Duration_ms;
%             V1(1,m)= 1;
%             V2(1,m) =  Params.Output_States(2).StimChans.Level.Level;% modified on 03/04/22
%         elseif isfield(Params, 'Tosca')
%             stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels.Gate.Duration_ms;
%             V1(1,m)= 1; % modified on 03/04/22
%             V2(1,m)= 1;
%         else
%             error('StimLength not found for noiseburst. Ask 2P Slack for help.')
%         end
%     end
end

%% 2 = Receptive Field

if setup.stim_protocol == 2
    
    [V1, V2, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        
        if isfield(Data{m},'Sound')
            
            %Frequency
            if isfield(Data{m}.Sound.Signal, 'Tone')
                V1(1,m)  = Data{m}.Sound.Signal.Tone.Frequency_kHz;
            elseif isfield(Data{m}.Sound.Signal, 'Waveform')
                V1(1,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
            elseif isfield(Params, 'Tosca')
                if size(Params.Tosca.Schedule.Families,2) > 1
                    %If true, user has more than one stim Group set up in schedules
                    %and the group missing Frequency is the blank/sham group
                    V1(1,m)  = 0; %Sham stim
                else
                    error('Frequency not found for Receptive Field. Ask 2P Slack for help.')
                end
            else
                error('Frequency not found for Receptive Field. Ask 2P Slack for help.')
            end
            
            %Intensity
            V2(1,m) = Data{m}.Sound.Signal.Level.dB_SPL;
            
            %StimLength
            if V2(1,m) <= 0
                stimLength(1,m) = 0;
            else
                if isfield(Params, 'Output_States')
                    stimLength(1,m) = Params.Output_States(2).StimChans.Burst.Duration_ms;
                elseif isfield(Params, 'Tosca')
                    stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels.Gate.Duration_ms;
                else
                    error('StimLength not found for Receptive Field. Ask 2P Slack for help.')
                end
            end
            
        elseif isfield(Data{m},'cue')
            
            %Frequency
            if isfield(Data{m}.cue.Signal,'Tone')
                V1(1,m)  = Data{m}.cue.Signal.Tone.Frequency_kHz;
            elseif isfield(Data{m}.cue.Signal,'Waveform')
                V1(1,m)  = Data{m}.cue.Signal.Waveform.Frequency_kHz;
            else
                V1(1,m) = 0; %Blank trials?
            end

            %Intensity
            V2(1,m) = Data{m}.cue.Signal.Level.dB_SPL;
            
            %StimLength
            if V2(1,m) <= 0
                stimLength(1,m) = 0;
            else
                if isfield(Params, 'Tosca')
                    stimLength(1,m) =  Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
                else
                    error('StimLength not found for Receptive Field. Ask 2P Slack for help.')
                end
            end
            
        else
            error('Receptive Field parameters not found. Ask 2P Slack for help.')
        end
    end
end

%% 3 = FM Sweep

if setup.stim_protocol == 3
    
    [V1, V2, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        
        if isfield(Data{m},'Sound')
            
            %Rate
            if isfield(Data{m}.Sound.Signal,'FMSweep')
                V1(1,m) = Data{m}.Sound.Signal.FMSweep.Rate_oct_s;
            elseif isfield(Data{m}.Sound.Signal, 'Sweep')
                V1(1,m) = Data{m}.Sound.Signal.Sweep.Rate;
            else
                V1(1,m) = 0; %Blank trial
            end

            %Intensity
            V2(1,m)  = Data{m}.Sound.Signal.Level.dB_SPL;

            %Stim Length
            if  V2(1,m) == 0
                stimLength(1,m) = 0;
            else
                speed = abs(V1(1,m));
                
                if isfield(Params, 'Output_States')
                    Fmin = Params.Output_States(2).StimChans.Waveform.FMSweep.Fmin_Hz;
                    Fmax = Params.Output_States(2).StimChans.Waveform.FMSweep.Fmax_Hz;
                elseif isfield(Params, 'Tosca')
                    Fmin = Params.Tosca.Flowchart(1).State.SigMan.Channels.Waveform.Fmin_Hz;
                    Fmax = Params.Tosca.Flowchart(1).State.SigMan.Channels.Waveform.Fmax_Hz;
                else
                    error('Fmax and Fmin not found for FM sweep. Ask 2P Slack for help.')
                end
                stimLength(1,m) = ((log2(Fmax) -log2(Fmin))./speed)*1000;
            end

        elseif isfield(Data{m},'cue')
            
            %Rate
            if isfield(Data{m}.cue.Signal,'Sweep')
                V1(1,m) = Data{m}.cue.Signal.Sweep.Rate;
            else
                V1(1,m) = 0; %Blank trials
            end
            
            %Intensity
            V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;
            
            %StimLength
            if  V2(1,m) == 0
                stimLength(1,m) = 0;
            else
                speed = abs(V1(1,m));
                % original FM sweep files -- signal is in the 3rd flowchart item
                if ~isempty(Params.Tosca.Flowchart(3).State.SigMan.Channels)
                    Fmin = Params.Tosca.Flowchart(3).State.SigMan.Channels.Waveform.Fmin_Hz;
                    Fmax = Params.Tosca.Flowchart(3).State.SigMan.Channels.Waveform.Fmax_Hz;
                else % Christine's new FM sweep files (with opto) signal is the second flowchart
                    Fmin = Params.Tosca.Flowchart(2).State.SigMan.Channels(1).Waveform.Fmin_Hz;
                    Fmax = Params.Tosca.Flowchart(2).State.SigMan.Channels(1).Waveform.Fmax_Hz;
                end 
                stimLength(1,m) = ((log2(Fmax) -log2(Fmin))./speed)*1000;
            end
        else
            error('FM sweep parameters not found. Ask 2P Slack for help.')
        end
    end
end
        
%% 4 = Widefield (NO LONGER USED)

if setup.stim_protocol == 4
    error('We have phased out stim_protocol = 4 for widefield. Use Receptive Field (stim_protocol = 2) instead.')
end

%% 5 = SAM

if setup.stim_protocol == 5
    
    [V1, V2, V3, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)

        if isfield(Data{m},'Sound')

            if isfield(Data{m}.Sound.Signal,'SAM')
                %stim trials...
                if isfield(Data{m}.Sound.Signal.SAM, 'Frequency_Hz')
                    V1(1,m)  = Data{m}.Sound.Signal.SAM.Frequency_Hz;
                elseif isfield(Data{m}.Sound.Signal.SAM, 'Rate_Hz')
                    V1(1,m)  = Data{m}.Sound.Signal.SAM.Rate_Hz;
                else
                    error('Rate not found for SAM. Ask 2P Slack for help.')
                end
                
                V2(1,m) = Data{m}.Sound.Signal.SAM.Depth_0_minus1;
                V3(1,m) = Data{m}.Sound.Signal.Level.dB_SPL;
                
                if isfield(Params, 'Output_States')
                    stimLength(1,m) = Params.Output_States(2).StimChans.Burst.Duration_ms;
                elseif isfield(Params, 'Tosca')
                    stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels.Gate.Duration_ms;
                else
                    error('StimLength not found for SAM. Ask 2P Slack for help.')
                end
            else
                %blank trials...
                V1(1,m) = NaN;
                V2(1,m) = NaN;
                V3(1,m) = 0;
                stimLength(1,m) = 0;
            end

        elseif isfield(Data{m},'cue')
                
            if isfield(Data{m}.cue.Signal, 'SAM')
                
                %Frequency/Rate
                if isfield(Data{m}.cue.Signal.SAM, 'Frequency_Hz') || isfield(Data{m}.cue.Signal.SAM, 'Rate_Hz')
                    %stim trials...
                    if isfield(Data{m}.cue.Signal.SAM, 'Frequency_Hz')
                        V1(1,m) = Data{m}.cue.Signal.SAM.Frequency_Hz;
                    elseif isfield(Data{m}.cue.Signal.SAM, 'Rate_Hz')
                        V1(1,m) = Data{m}.cue.Signal.SAM.Rate_Hz;
                    else
                        error('Rate not found for SAM. Ask 2P Slack for help.')
                    end
                else
                    V1(1,m) = NaN; %No rate data
                end
                
                %Modulation depth
                if isfield(Data{m}.cue.Signal.SAM, 'Depth_0_minus1')
                    V2(1,m) = Data{m}.cue.Signal.SAM.Depth_0_minus1;
                else
                    V2(1,m) = NaN; %No depth data
                end
                
                %Intensity
                if isfield(Data{m}.cue.Signal.SAM, 'Level')
                    V3(1,m) = Data{m}.cue.Signal.SAM.Level.dB_SPL;
                else
                    V3(1,m) = Data{m}.cue.Signal.Level.dB_SPL;
                end
                
                %StimLength
                if isfield(Params, 'Output_States')
                    stimLength(1,m) = Params.Output_States(2).StimChans.Burst.Duration_ms;
                elseif isfield(Params, 'Tosca')
                    if  isfield(Params.Tosca.Flowchart(3).State.SigMan.Channels,'Gate')
                        stimLength(1,m) = Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
                    else
                        stimLength(1,m) = Params.Tosca.Flowchart(2).State.SigMan.Channels(1).Gate.Duration_ms;
                    end 
                else
                    error('StimLength not found for SAM. Ask 2P Slack for help.')
                end
                
            else
                %blank trials...
                V1(1,m) = NaN;
                V2(1,m) = NaN;
                V3(1,m) = Data{m}.cue.Signal.Level.dB_SPL;
                stimLength(1,m) = 0;
            end
            
        end
    end      
end

%% 6 = SAM freq

if setup.stim_protocol == 6 %SAM freq
    
    [V1, V2, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        
        rateFound = true;
        
        if isfield(Data{m}.Sound, 'Signal')
            if isfield(Data{m}.Sound.Signal, 'Tone')
                V1(1,m)  = Data{m}.Sound.Signal.Tone.Frequency_kHz;
            elseif isfield(Data{m}.Sound.Signal, 'Waveform')
                V1(1,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
            else
                rateFound = false;
            end
        end
        
        if rateFound
            %stim trials
            V2(1,m) = Data{m}.Sound.Signal.SAM.Depth_0_minus1;
            if isfield(Params, 'Output_States')
                stimLength(1,m) = Params.Output_States(2).StimChans.Burst.Duration_ms;
            elseif isfield(Params, 'Tosca')
                stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels.Gate.Duration_ms;
            else
                error('StimLength not found for SAM Freq. Ask 2P Slack for help.')
            end
        else
            %blank trials
            V1(1,m)  = nan;
            V2(1,m)  = nan;
            stimLength(1,m) = 0;
        end
        
    end
end
        
%% 7 = Behavior go/nogo freq. disc.

if setup.stim_protocol == 7
    
    [V1, V2, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        if isfield(Data{m}.cue.Signal,'Waveform') == true
            V1(1,m)  = Data{m}.cue.Signal.Waveform.Frequency_kHz;
            V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;
            stimLength(1,m) = Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
        else
            %V1(1,m)  = Params.Output_States(2).StimChans.Stimulus.Waveform.Tone.Frequency_kHz;
            %V2(1,m)  = Params.Output_States(2).StimChans.Stimulus.Level.Level;
            V1(1,m)  = Data{m}.cue.Signal.Tone.Frequency_kHz;
            V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;
            stimLength(1,m) =  Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
        end
    end
end
    
%% 8 = ABI
    
if setup.stim_protocol == 8
    
    [V1, V2, V3, stimLength] = deal(nan(1,length(Data)));
    if Params.Info.Version > 2366
        statename= Params.Tosca.Flowchart(3).State.Name;%Addition of ABI for VAD (09/22/23)
    else 
        statename=Params.Flowchart(3).Name;%Addition for old ABI cohort for VAD (11/28/23)
    end
        
    %Figure out statename
%     statename = [];
%     while isempty(statename)
%         for m = 1:length(Data)
%             statenames = string({Data{m}.states.name});
%             stateOI = find(statenames == 'cue' | statenames == 'ABI', 1, 'first');
%             statename = statenames(stateOI);
%         end
%     end
%     
    for m = 1:length(Data)
        %Stim Length and check that ABI stim was Pulse train
        if isfield(Params, 'Output_States')
            stimLength(1,m) = Params.Output_States(2).StimChans.Burst.Duration_ms;
            Stim_type = Params.Output_States(2).StimChans.Waveform.Type;
        elseif isfield(Params, 'Tosca')
            stimLength(1,m) =  Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
            Stim_type = Params.Tosca.Flowchart(3).State.SigMan.Channels.Waveform.Class;
        else
            error('StimLength not found for ABI. Ask 2P Slack for help.')
        end
        
        
        if strcmp(Stim_type, 'Pulse Train')
            
            % check and extract varied stim parameters V1 = Amplitude, V2= Rate, V3= Duration
            if isfield(Params, 'Output_States')
                V1(1,m) = Params.Output_States(2).StimChans.Level.Level;
            elseif isfield(Data{m}.(statename).Current.Level,'dB_Vrms')% Pulse amplitude varied
                V1(1,m)  = Data{m}.(statename).Current.Level.dB_Vrms;
            elseif isfield(Data{m}.(statename).Current.Level, 'dB_re_1_Vrms')
                V1(1,m)  = Data{m}.cue.Current.Level.dB_re_1_Vrms;% (old file)
            elseif isfield(Data{m}.(statename).Current.Level,'Volts')
                V1(1,m)  = Data{m}.(statename).Current.Level.Volts;
            else
                V1(1,m) =  Params.Tosca.Flowchart(3).State.SigMan.Channels.Level.Level;
            end

            %if isfield(Data{m}.cue.Current,'Waveform') == true % Pulserate varied (old Files)
            if isfield(Data{m}.(statename).Current,'Pulse')% Pulse rate varied
                V2(1,m)= Data{m}.(statename).Current.Pulse.Rate_Hz;
            elseif isfield(Data{m}.cue.Current,'Waveform') == true
                V2(1,m)= Data{m}.cue.Current.Waveform.Pulse_rate_Hz;% (old Files)
            else
                if isfield(Params, 'Output_States')
                V2(1,m)=Params.Output_States(2).StimChans.Waveform.Pulse_Train.Pulse_rate_Hz;% (old Files)
                else
                V2(1,m)= Params.Tosca.Flowchart(3).State.SigMan.Channels.Waveform.Rate_Hz;
                end
            end
            
            if isfield(Data{m}.(statename).Current,'Gate')% Burst Duration varied
                V3(1,m)= Data{m}.(statename).Current.Gate.Duration_ms;
            elseif isfield(Params.Output_States(2).StimChans.Burst, 'Duration_ms')
                V3(1,m)=Params.Output_States(2).StimChans.Burst.Duration_ms;%(old Files)
            else
                V3(1,m)=Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
            end
            
        else
            error('Not ABI Pulse Train. Ask 2P Slack for help.')
        end
    end
end

%% 16 = ABI SAM
if setup.stim_protocol == 16
    
    [V1, V2, V3, stimLength] = deal(nan(1,length(Data)));
    %Figure out statename (From Maryse)
    %     statename = [];
    %     while isempty(statename)
    %         for m = 1:length(Data)
    %             statenames = string({Data{m}.states.name});
    %             stateOI = find(statenames == 'cue' | statenames == 'ABI', 1, 'first');
    %             statename = statenames(stateOI);
    %         end
    %     end
    if isfield(Params, 'Output_States')
        for m = 1:length(Data)
                statenames = string({Data{m}.states.name});
                pre_stateOI = find(statenames == 'cue' | statenames == 'ABI', 1, 'first');
                stateOI = statenames(pre_stateOI);
        end
    else
        stateOI = Params.Tosca.Target_state; %Addition of ABI for VAD (09/22/23)
    end
    
    for m = 1:length(Data)
        if ~strcmp(Data{1,m}.Group, 'Group1')  %blank trials
            V1(1,m)  = NaN;
            V2(1,m)  = NaN;
            V3(1,m)  = NaN;
            
        else
            %Stim Length
            if isfield(Params, 'Output_States')
                stimLength(1,m) = Params.Output_States(2).StimChans.Burst.Duration_ms;
            elseif isfield(Params, 'Tosca')
                stimLength(1,m) =  Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
            else
                error('StimLength not found for ABI SAM. Ask 2P Slack for help.')
            end
            
            %Stim params varied (V1=AM freq, V2= Mod. Depth, V3= CurrentLVL)
            %Frequency/Rate
            if isfield(Data{1,m}.(stateOI).Current.SAM, 'Frequency_Hz') || isfield(Data{1,m}.(stateOI).Current.SAM, 'Rate_Hz')|| isfield(Data{1,m}.(stateOI).Current.Pulse, 'Rate_Hz')
                %stim trials...
                if isfield(Data{1,m}.(stateOI).Current.SAM, 'Frequency_Hz')
                    V1(1,m) = Data{1,m}.(stateOI).Current.SAM.Frequency_Hz;
%                     V1(1,m)  = Data{m}.(stateOI).Current.Pulse.Frequency_Hz;%(Old Files)
                elseif isfield(Data{1,m}.(stateOI).Current.SAM, 'Rate_Hz') ||isfield(Data{1,m}.(stateOI).Current.Pulse, 'Rate_Hz')
                    V1(1,m) = Data{1,m}.(stateOI).Current.SAM.Rate_Hz;
%                     V1(1,m) = Data{1,m}.(stateOI).Current.Pulse.Rate_Hz;%(Old Files)
                else
                    error('Rate not found for ABI SAM. Ask 2P Slack for help.')
                end
            else
                V1(1,m) = NaN; %No rate data
                stimLength(1,m) =0;
            end
            
            %Modulation depth
            if isfield(Data{1,m}.(stateOI).Current.SAM, 'Depth_0_minus1')
                V2(1,m) = Data{1,m}.(stateOI).Current.SAM.Depth_0_minus1;
            else
                V2(1,m) = NaN; %No depth data
            end
            
            %Intensity
            if isfield(Data{1,m}.(stateOI).Current.Level, 'dB_Vrms')
                V3(1,m)  = Data{1,m}.(stateOI).Current.Level.dB_Vrms;
            elseif isfield(Data{1,m}.(stateOI).Current.Level, 'dB_re_1_Vrms')
                V3(1,m)  = Data{1,m}.(stateOI).Current.Level.dB_re_1_Vrms;
            elseif isfield(Data{1,m}.(stateOI).Current.Level, 'Volts')
                V3(1,m)  = Data{1,m}.(stateOI).Current.Level.Volts;
            end
            
            %        %Accommodate files with different statenames (From Maryse)
            %         switch statename
            %             case 'cue'
            %                 %V1(1,m)  = Data{m}.cue.Current.SAM.(field); (Old Files)
            %                 V1(1,m)  = Data{m}.(statename).Current.Pulse.(field);
            %                 V2(1,m)  = Data{m}.(statename).Current.SAM.Depth_0_minus1;
            %
            %             case 'ABI'
            %                 V1(1,m)  = Data{m}.(statename).Current.SAM.Frequency_Hz;
            %                 V2(1,m)  = Data{m}.(statename).Current.SAM.Depth_0_minus1;
            %         end
            
            %V3(1,m)  = Params.Output_States(2).StimChans.Waveform.Pulse_Train.Pulse_rate_Hz;
        end
    end
end


%% 9, 11, 111 = Random H20 or Air Puffs

if setup.stim_protocol == 9 || setup.stim_protocol == 11 || setup.stim_protocol == 111
    
    V1 = zeros(1,length(Data)); %V1 indicates if air/water was on (1) or off (0)
    V2 = zeros(1,length(Data)); %V2 is not used for anything for these stim
    stimLength = zeros(1,length(Data)); %stimLength indicates duration that air/water solenoid were open (may be calibrated per experiment)
    paramaterFound = false;
    
    for m = 1:length(Data)
          
        %Check to see if H2O is varied within each trial (e.g. user varied H2O duration within schedule)
        %If it is, use this to determine V1
        if isfield(Data{m}, 'Start')
            if isfield(Data{m}.Start.Timeout,'H2O')
                stimLength(1,m) = Data{m}.Start.Timeout.H2O.Duration_ms;
                paramaterFound = true;
                
                if stimLength(1,m) > 0
                    V1(1,m) = 1;
                end
            
                %Continue to next trial
                continue;
            end
        elseif isfield(Data{m}, 'Sound')
            if isfield(Data{m}.Sound.Timeout,'H2O')
                stimLength(1,m) = Data{m}.Sound.Timeout.H2O.Duration_ms;
                paramaterFound = true;
                
                if stimLength(1,m) > 0
                    V1(1,m) = 1;
                end
            
                %Continue to next trial
                continue;
            end
        end
        
        %If H2O is not varied within trial, program may have used CS+ vs CS- to determine water delivery
        %The exception is Christine's data, where there was only CS-. This data had Target_state == End
        if isfield(Params, 'Tosca') %New Tosca
            if strcmp(Params.Tosca.Target_state, 'End')
                %Christine data
                V1(1,m) = 1;
                paramaterFound = true;
            end
        else
            %Old Tosca: Carolyn and Maryse data
            if strcmp(Data{m}.Type, 'CS+')
                V1(1,m) = 1;
                paramaterFound = true;
            elseif strcmp(Data{m}.Type, 'CS-')
                V1(1,m) = 0;
                paramaterFound = true;
            end
        end

        %Check to see whether paramter was found in any of the above conditions
        if ~paramaterFound
            error('Water/Air value not found. Ask 2P Slack for help.')
        end
        
        %Now get StimLength from Params and apply to all trials where there was water
        if V1(1,m) == 1
            if isfield(Params, 'Tosca')
                stimLength(1,m) = Params.Tosca.Flowchart(1).State.Timeouts(1).TTL.Stop.Expr;
            elseif isfield(Params, 'Flowchart')
                stimLength(1,m) = Params.Flowchart(3).Timeouts(2).TTL.Stop.Expr;
            else
                error('StimLength not found for random air/water. Ask 2P Slack for help.')
            end
        else
            stimLength(1,m) = 0;
        end
    end
end
    
%% 10 = Noiseburst_ITI

if setup.stim_protocol == 10
    
    V2 = zeros(1,length(Data)); %We used to store ITI here, but we were not using it
    [V1, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        
        %Intensity
        if isfield(Data{m},'Sound')
            if isfield(Data{m}.Sound,'Signal')
                V1(1,m)  = Data{m}.Sound.Signal.Level.dB_SPL; %0dB for no stim, 70dB for stim
            elseif isfield(Data{m}.Sound,'Silence')
                V1(1,m)  = Data{m}.Sound.Silence.Level.dB_SPL; %0dB for no stim, 70dB for stim
            else
                error('Stim level not found for NoiseITI. Ask 2P Slack for help.')
            end
        elseif isfield(Data{m},'cue')
            V1(1,m)  = Data{m}.cue.Signal.Level.dB_SPL; %0dB for no stim, 70dB for stim
        else
            V1(1,m) = Params.Output_States(2).StimChans.Level.Level; % programmed without blank trials
        end
        
        %StimLength
        if V1(1,m) <= 0
            stimLength(1,m) = 0;
        else
            if isfield(Params,'Output_States')
                stimLength(1,m) = Params.Output_States(2).StimChans.Burst.Duration_ms;
            elseif isfield(Params, 'Tosca') 
                %I can't figure out how to get rid of this Try/Catch block yet
                try
                    stimLength(1,m) =  Params.Tosca.Flowchart(1).State.SigMan.Channels.Gate.Duration_ms;
                catch
                    stimLength(1,m) =  Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
                end
            else
                error('StimLength not found for NoiseITI. Ask 2P Slack for help.')
            end
        end
    end
end
    
%% 12 = Spontaneous

if setup.stim_protocol == 12 
    V1 = 0;
    V2 = 0;
    stimLength = nan(1,length(Data));
    
    for m = 1:length(Data)
        stimLength(1,m) = 0;
    end
end

%% 13 = Maryse behavior

if setup.stim_protocol == 13
    
    [V1, V2, V3] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        
        %Stim frequency
        if isfield(Data{m}, 'Target_kHz') && isfield(Data{m}, 'Standard_kHz')
            V1(1,m) = Data{1,m}.Standard_kHz;
            V2(1,m) = Data{1,m}.Target_kHz;
        elseif isfield(Data{m}, 'FreqDiscrimTrainHold')
            V1(1,m) = nan;
            V2(1,m) = Data{1,m}.FreqDiscrimTrainHold;
        else
            error('Target frequency not found for Maryse behavior stim. Ask 2P Slack for help.')
        end
    end
    
    %Store stim intensity in V3        
    if isfield(Params, 'Output_States')
        try
            stim_level = Params.Output_States(2).StimChans(1).Stimulus.Level.Level;
        catch
            stim_level = Params.Output_States(2).StimChans(1).Level.Level;
        end
    elseif isfield(Params, 'Tosca')
        try
            stim_level = Params.Tosca.Schedule.Families.Vars(3).expr;
        catch
            stim_level = Params.Tosca.Flowchart(3).State.SigMan.Channels(1).Channel.Level.Level;
        end
    end
    V3(:) = stim_level;
    
    targetFreq = mode(V2);
    if any(isnan(V1))
        V1(:) = targetFreq;
    end
end

%% 15 = Auditory Ripple

if setup.stim_protocol == 15 % auditory ripple
    
    [V1, V2, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        if isfield(Data{1,m}, 'Sound')
            if isfield(Data{1,m}.Sound.Signal, 'Ripple')
                V1(1,m) = Data{1,m}.Sound.Signal.Ripple.Density_cyc_oct;
                V2(1,m) = Data{1,m}.Sound.Signal.Level.dB_SPL;
                stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels.Gate.Duration_ms;
            else
                V1(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels.Waveform.Density_cyc_oct;
                V2(1,m) = Data{1,m}.Sound.Signal.Level.dB_SPL;
                stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels.Gate.Duration_ms;
            end
        else
            if exist('Data{1,m}.cue.Signal.Ripple')
                V1(1,m) = Data{1,m}.cue.Signal.Ripple.Density_cyc_oct;
            elseif isfield(Params.Tosca.Flowchart(3).State.SigMan.Channels,'Waveform')
                V1(1,m) = Params.Tosca.Flowchart(3).State.SigMan.Channels.Waveform.Density_cyc_oct;
            else % new Ripple with opto -- signal waveform info started in the 2nd flowchart
                V1(1,m) = Params.Tosca.Flowchart(2).State.SigMan.Channels(1).Waveform.Density_cyc_oct;
            end
            V2(1,m) = Data{1,m}.cue.Signal.Level.dB_SPL;
            
            if isfield(Params.Tosca.Flowchart(3).State.SigMan.Channels,'Gate')
                stimLength(1,m) =  Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;
            else % new Ripple with opto -- signal waveform info started in the 2nd flowchart
                stimLength(1,m) =  Params.Tosca.Flowchart(2).State.SigMan.Channels(1).Gate.Duration_ms;
            end
        end
    end
end

%% 18 = vocalizations

if setup.stim_protocol == 18
    
    [V1, V2, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        if isfield(Data{1,m}, 'Sound')
            if isfield(Data{1,m}.Sound.Signal, 'File')
                if isfield(Data{1,m}.Sound.Signal.File, 'Token')
                    V1(1,m) = Data{1,m}.Sound.Signal.File.Token;
                    % go digging for level at which sounds were presented
                    V2(1, m) = Params.Tosca.Flowchart(1).State.SigMan.Channels(1).Level.Level;
                    stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels(1).Gate.Duration_ms ;
                else
                    error('Vocalization Token not found.')
                end
            else
                %This is an addition for the blank trials in vocal stimuli
                V1(1,m) = 0;
                V2(1,m) = 0;
                stimLength = 0;
            end
        else 
            error('Vocalization stimulus parameters incorrectly identified.')
        end
    end
end

%% 1018 = vocalizations + noiseburst
% ZG-defined new parameter file for compile block which includes
% presentation of blanks, vocalizations and noise all in one tosca
% parameter file. In stim_parameter 1018, V1 instead of 'noise' is being
% set equal to 101810 to avoid teh issue of having a string as a V1
% variable. Since this cases a problem for this error: Unable to perform assignment because 
% the size of the left side is 1-by-1 and the size of the right side is 1-by-5.

if setup.stim_protocol == 1018

    [V1, V2, stimLength] = deal(nan(1,length(Data)));

    %V1 = string(V1)

    for m =1:length(Data)
        if matches(Data{1,m}.Group, 'Token')
            V1(1,m) = Data{1,m}.Sound.Signal.File.Token; % gets V1 (stim ID) from where the vocal ID is stored
            V2(1, m) = Params.Tosca.Flowchart(1).State.SigMan.Channels(1).Level.Level;% gets V2 (sound level) from where the SPL output is stored
            stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels(1).Gate.Duration_ms; %extracts the stim length

        elseif matches(Data{1,m}.Group, 'noise')
            V1(1,m) = 101810;
            %gave noise a value to resolve the string issue
            V2(1, m) = Params.Tosca.Flowchart(1).State.SigMan.Channels(2).Level.Level;
            stimLength(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels(2).Gate.Duration_ms;

        elseif matches(Data{1,m}.Group, 'blanks')
            V1(1,m) = 0;
            V2(1,m) = 0;
            stimLength = 0;

        else
            error('stimulus group not represented')
        end
    end
end
%% 19 = AM detection behavior (modeled after Maryse Behavior #13)

if setup.stim_protocol == 19
    
    [V1, V2, V3] = deal(nan(1,length(Data)));
    stateOI = Params.Tosca.Target_state; %Addition of ABI for VAD (10/06/23)
    for m = 1:length(Data)
        V1(1,m) = 0; %Standard depth = 0%
        V2(1,m) = Data{1,m}.ModDepth; %Target AM depth
        if isfield(Data{1,m}.Any.Signal.Level, 'dB_Vrms') %accomodate for different units and accoustic / ABI
            V3(1,m) = Data{1,m}.Any.Signal.Level.dB_Vrms;
        elseif isfield(Data{1,m}.Any.Signal.Level, 'Volts')
            V3(1,m) = Data{1,m}.Any.Signal.Level.Volts;
        else
            V3(1,m) = Data{1,m}.Any.Signal.Level.dB_SPL; %Stim intensity
        end
    end
end
%% 21 = ABC/XYZ pattern behavior (21) and randomized ABCXYZ stim (211)
if any(setup.stim_protocol == [21 211])
    clear V1 V2 V3
 
    for m = 1:length(Data)
        groupsall(m) = convertCharsToStrings(Data{1,m}.Group);
        check(m,1) = Data{1,m}.Sound1.Signal.Sweep.Rate;
        check(m,2) = Data{1,m}.Sound2.Signal.Sweep.Rate;
        check(m,3) = Data{1,m}.Sound3.Signal.Sweep.Rate;
        if Data{1,m}.Group =='C';
            if check(m,3) == 40
                V1(1,m) = string('C1');
            elseif check(m,3) == -10
                V1(1,m) = string('C2');
            else
                error('missing C variable')
            end
        elseif Data{1,m}.Group == 'Z';
            if check(m,3) == -40
                V1(1,m) = string('Z1');
            elseif check(m,3) == 40
                V1(1,m) = string('Z2');
            else
                error('missing Z variable')
            end
        else
            V1(1,m) = convertCharsToStrings(Data{1,m}.Group);
            
        end
        V2(1,m) = convertCharsToStrings(strcat(num2str(check(m,1)),'_',num2str(check(m,2)),'_',num2str(check(m,3)))); % put the sweeps here, string together.
        V3(1,m) = Params.Tosca.Flowchart(1).State.SigMan.Channels.Level.Level;
        
        if Data{1,m}.Sound1.Signal.Gate.Duration_ms == 0 % in early versions the z stimulus doesnt have a gate length recorded.
            stimLength(1,m) = 1;
        else
            stimLength(1,m) = (Data{1,m}.Sound1.Signal.Gate.Duration_ms+Data{1,m}.Sound2.Signal.Gate.Duration_ms+Data{1,m}.Sound3.Signal.Gate.Duration_ms); 
        end
    end
end

%% 33 = Fear conditioning with FM sweeps

if setup.stim_protocol == 33
    
    [V1, V2, V3, stimLength] = deal(nan(1,length(Data)));
    
    for m = 1:length(Data)
        if isfield(Data{m}.cue.Signal,'Sweep') == true
            V1(1,m)  = Data{m}.cue.Signal.Sweep.Rate;
            V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;
            stimLength(1,m) = Params.Tosca.Flowchart(2).State.SigMan.Channels.Gate.Duration_ms;
            if Data{m}.Type == 'CS-'
                V3(m) = 0; 
            elseif Data{m}.Type == 'CS+'
               V3(m) = 1;
            else
                error('no cs+ or cs- in data')
            end
%         else
%             %V1(1,m)  = Params.Output_States(2).StimChans.Stimulus.Waveform.Tone.Frequency_kHz;
%             %V2(1,m)  = Params.Output_States(2).StimChans.Stimulus.Level.Level;
%             V1(1,m)  = Data{m}.cue.Signal.Tone.Frequency_kHz;
%             V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;
%             stimLength(1,m) =  Params.Tosca.Flowchart(3).State.SigMan.Channels.Gate.Duration_ms;;
        end
    end
end

%% None of the above

if isempty(V1) && isempty(V2)
    warning(['stim_protocol ' num2str(setup.stim_protocol) ' does not exist yet'])
end

if isempty(stimLength) && ~any(setup.stim_protocol == [13 19]) %For Freq/AM detection behavior, stimLength is defined outside of this script
    warning('StimLength not defined. You may run into errors later if missing this variable. Ask 2P Slack for help.')
end

end %End function
