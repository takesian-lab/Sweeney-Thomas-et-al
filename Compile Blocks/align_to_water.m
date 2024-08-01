function [wTime,block] = align_to_water(block, bout_detection_method)
% for random water presentation, often, in the lack of a cue, a mouse will
% lick the water drop after the initial presentation of the drop. Here, we
% find the times associated with lick bouts rather than the times that a
% trial was initiated. 

% Inputs
% - bout_detection_method (optional) 1 to use moving mean of licks, 2 to use counts

% Version history:
% - V1 = See archived align_to_water_v1
% - V2 = Archived 4/27/2023: alignment was not working well
% - V3 = Current version corrects V2

%% find which trials give water vs no water.

if nargin < 2
    bout_detection_method = 1;
end

if block.setup.stim_protocol == 9 % random/uncued water (9)
    %For protocol 9, this function will return a new Sound_Time where hit and
    %false alarm trials are realigned to licks
    %v1 will remain unchanged
    %v2 will contain the trial outcome
    %Sound_Time will be updated
    v1 = block.parameters.variable1;
elseif block.setup.stim_protocol == 7 % carolyn's frequency stim paradigm
    %For protocol 7, this function will compute a wTime that corresponds
    %to the time of the first lick after each trial (no lick = nan)
    %Sound_Time will remain unchanged
    v1 = block.Outcome == 1;
end

Sound_Time = block.Sound_Time; 
Sound_Time_orig = block.Sound_Time; %save for plotting

%when realigning trials to lick, minimum duration of time required before next trial
minTrialDuration = 2.5; %seconds

%% find the lick bouts [TWO METHODS]
licks = block.concat_licks;
timetrace = block.concat_times;


if bout_detection_method == 1
    % change licks to a continuous trace
    fs = round(1/((timetrace(end)-timetrace(1))/length(timetrace)));
    lickrate = movmean(licks,fs/5); %Perform moving mean over each 200ms of data
    lickthreshold = 0.1; %Arbitrary threshold for saying whether mouse was licking or not
    isLicking = lickrate > lickthreshold;
    lick_timestamp = timetrace(isLicking);
    lickIDX = find(isLicking==1)';
else
    % set the licks into bouts 
    lick_timestamp = timetrace(find(licks==1));
    lickIDX = find(licks==1);

    count = 1;
    difflick = diff(lickIDX);

    for i = 1:length(difflick)
        if i ==1
            bout_start(count,1) = lick_timestamp(1);
            bout_start(count,2) = lickIDX(1);
        end

        if i == length(difflick)
            bout_end(count,1) = lick_timestamp(i);
            bout_end(count,2) = lickIDX(i);
        end

        if difflick(i)>50
            if count == 1 % get the end of the previous bout
                bout_end(count,1) = lick_timestamp(i+1);
                bout_end(count,2) = lickIDX(i+1);
            end

            count = count+1;
            bout_start(count,1) = lick_timestamp(i+1);
            bout_start(count,2) = lickIDX(i+1);

            bout_end(count-1,1) = lick_timestamp(i);
            bout_end(count-1,2) = lickIDX(i);

        end
    end
end

%% Get time of first lick after each trial

wTime = nan(size(Sound_Time));
if ~isempty(lick_timestamp)
nextTrialTimes = [Sound_Time(2:end),lick_timestamp(end)];
for i = 1:length(Sound_Time)

    %Find the time of the first lick between this trial and the next
    T = nextTrialTimes(i); %Time of the following trial or end of block
    t = find((lick_timestamp>Sound_Time(i)) + (lick_timestamp<T) == 2,1,'First'); 

    if ~isempty(t)
        wTime(i) = lick_timestamp(t);
    end
end
else
    warning('no lick timestamps')
wTime = nan(size(Sound_Time));
end


%% Determine all conditions we want the trials to meet

%HIT TRIAL = V1 = 1, V2 = 1
%water trial where mouse licks
%we will move the start time to be aligned with licking
% -if we cannot do this and make sure the trial is at least 2.5s long -> NaN trial

%MISS TRIAL = V1 = 1, V2 = 0
%water trial where mouse didn't lick
% -don't do anything

%WITHHOLD TRIAL = V1 = 0, V2 = 3
%sham trial where mouse didn't lick
% -don't do anything 

%FALSE ALARM TRIAL = V1 = 0, V2 = 4
%sham trial where mouse licked
% - confirm there was no water (e.g. from prior miss trial) -> keep as motor sham trial
% - if there was water from prior miss trial -> NaN trial

%NAN TRIAL = V1 = 0/1, V2 = NaN
%If they don't meet any of the conditions above

%FOR ANALYSIS, WE WILL BE PRIMARILY INTERESTED IN THE HIT AND WITHHOLD TRIALS

%% Step through each trial and determine its outcome

if block.setup.stim_protocol == 9
    
    outcome = nan(size(v1)); %Record Hit, Miss, FP, Withold in outcome
    
    waterOnTheSpout = 0; %Keep track of whether water is on the spout
    
    for i = 1:length(outcome)
        
        %Water trials
        if v1(i) == 1
            
            if ~isnan(wTime(i))
                %HIT -> now make sure the trial is long enough
                if nextTrialTimes(i) - wTime(i) > minTrialDuration
                    outcome(i) = 1; %successful hit
                    Sound_Time(i) = wTime(i); %update sound_time
                else
                    outcome(i) = nan; %Not enough time to do realignment
                end
                
                waterOnTheSpout = 0; %Mouse licked water off
            else
                %MISS
                outcome(i) = 0; %Unlicked water
                
                waterOnTheSpout = 1; %Mouse missed water
            end

        %Sham trials
        else 
            if ~isnan(wTime(i))
                %FALSE ALARM
                
                %first check if there was water on the spout when the false alarm happened
                if waterOnTheSpout
                    %Cannot be a true false alarm, there was water
                    outcome(i) = nan;
                    waterOnTheSpout = 0; %Mouse licked it off here
                else
                    FAtime = wTime(i) - Sound_Time(i);
                    if FAtime > block.setup.constant.after_stim
                        %You can still use this as a sham trial, lick was after needed trial duration
                        outcome(i) = 3; %WITHHOLD
                    else
                        %See if trial is long enough to realign and use as motor sham
                        if nextTrialTimes(i) - wTime(i) > minTrialDuration
                            outcome(i) = 4; %FALSE ALARM
                            Sound_Time(i) = wTime(i); %update sound_time
                        else
                            outcome(i) = nan; %Not enough time to do realignment
                        end
                    end
                end
            else
                %WITHHOLD
                outcome(i) = 3; %Unlicked sham
                %Don't change waterOnTheSpout (identity stays with previous trial)
            end
        end
    end

    %v1 will stay the same and represent solenoid being on or off
    %v2 will represent trials with intended outcome (i.e. 1 = true hit or sham, 0 = something else)
    v2 = zeros(size(v1)); 
    v2(outcome == 1) = 1;
    v2(outcome == 3) = 1;
    
    wTime = Sound_Time; %Replace wTime with corrected Sound_Time
    block.Sound_Time = Sound_Time;
    block.parameters.variable2 = v2; 
    block.parameters.trialsToIgnore = isnan(outcome); %Trials with unspecified outcome
    block.Outcome = outcome; %FA and Miss info is stored in outcome
end


%% look at data

plot_figures = 0;

if plot_figures

    figure;
        
    subplot(2,1,1);
    %Set waterYtime and waterNtime
    waterYtime = Sound_Time_orig(v1 == 1);
    waterNtime = Sound_Time_orig(v1 == 0);
    
    plot(timetrace,lickrate); hold on
    scatter(lick_timestamp, ones(size(lick_timestamp)), 'filled')
    vline(waterYtime,'g-')
    vline(waterNtime,'r-')
    title('Original trials')
    ylabel('Lick rate')
    
    subplot(2,1,2); hold on
    %Set trial identities
    hit_time = Sound_Time(block.Outcome == 1);
    miss_time = Sound_Time(block.Outcome == 0);
    withhold_time = Sound_Time(block.Outcome == 3);
    FA_time = Sound_Time(block.Outcome == 4);
    nan_time = Sound_Time(isnan(block.Outcome)); 
    
    lickrate = movmean(licks,fs/5);
    plot(timetrace,lickrate); hold on
    if ~isempty(hit_time);      vline(hit_time,'g-');       end
    if ~isempty(withhold_time); vline(withhold_time,'r-');  end
    if ~isempty(miss_time);     vline(miss_time,'c-');      end
    if ~isempty(FA_time);       vline(FA_time,'m-');        end
    if ~isempty(nan_time);      vline(nan_time,'k-');       end
    title('Realigned')
    ylabel('Lick rate')

end
