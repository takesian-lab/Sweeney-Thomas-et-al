function Sound_Time = get_Sound_Time_behavior(block, IncludeErrorData)
%Get Sound_Time from block (e.g. if no Bruker or FibPhot data)
%IncludeErrorData = 0 by default
%Set IncludeErrorData to 1 to include Sound Times for error trials
% -- these are approximated because in some error trials the sound state is never reached
% -- that happens in in define_behavior_singleblock and is read from the errorData here

%% 
if nargin < 2
    IncludeErrorData = 0;
end

data = block.errorData.Tosca_times;
starts = block.errorData.start_time;
errors = block.errorData.error_trials;
sounds = block.errorData.New_sound_times;

if ~IncludeErrorData
    %remove error trials
    data(errors) = [];
    starts(errors) = [];
    sounds(errors) = [];
end

catsound_idx = [];
for i = 1:length(data)
    addzeros = zeros(length(data{i}),1);
    zerotrial = data{i} - starts(i);
    Sound_idx = find(zerotrial == sounds(i),1,'first');
    addzeros(Sound_idx) = 1;
    catsound_idx = [catsound_idx; addzeros];
end

Sound_Time = block.concat_times(catsound_idx == 1)';

if length(Sound_Time) ~= length(sounds)
    error('These variables should be the same length... Investigate')
end

end %function