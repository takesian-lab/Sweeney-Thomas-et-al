function [VocalData, fig1] = simple_compute_vocalization(isResponsive, PeakData, StimInfo, plot_figures)
% Analyze responses to vocalizations
%
% Argument(s):
%   isResponsive - array of 0s and 1s indicating whether the
%   cell is significantly responsive to that stim combination
%   PeakData - from simple_check_if_responsive
%   StimInfo - from simple_prepare_variables
%   plot_figures - 0 or 1 to plot figures
%
% Returns:
%   VocalData table:
%     - BestVocal = best vocalization according to highest amplitude response
%     - VI = valence index (preference for pos vs. negative valence)
%     - VI_p = p value for test against shuffled distribution
%     - VI_p_bs = p value for comparing bootstrapped distribution to 0
%     - FI = frequency index (preference for high vs. low frequency)
%     - FI_p = p value for test against shuffled distribution
%     - FI_p_bs = p value for comparing bootstrapped distribution to 0
%     - FRA_area = % of RF == 1
%     - Raw_BestVocal = computed ignoring isRF (best values not restricted to resposive & reliable)


% Search 'TODO'
%% Initial parameters

%% Initial user-defined parameters


nShuffles = 1000;
alpha = 0.05;

%Figures will return as [] if plot_figures is 0
fig1 = [];

%VocalData will return as NaN if there are no isResponsive stim conditions
VocalData = table;
VocalData.BestVocal = single(nan);
VocalData.VI = single(nan);
VocalData.VI_p = single(nan);
VocalData.VI_p_bs = single(nan);
VocalData.FI = single(nan);
VocalData.FI_p = single(nan);
VocalData.FI_p_bs = single(nan);
VocalData.FRA_area = single(nan);
VocalData.Raw_BestVocal = single(nan);
VocalData.Response = single(nan(1, size(StimInfo.V1_Label,1)));

%Check RF for responsive stim conditions and make sure there is >1 stim condition:
if ~any(isResponsive) || isequal(size(isResponsive),[1,1])
    return
end

%% compute FRA area using isResponsive

%Add NaNs to IsRF where there were no trial to analyze (e.g. in Running vs. NotRunning)
type = PeakData.ResponseType;
IsRF_withNans = single(isResponsive);
IsRF_withNans(strcmp(type, "undetermined")) = NaN;

%Compute FRA on non-NaN values only
IsRF_removeNans = IsRF_withNans;
IsRF_removeNans(isnan(IsRF_removeNans)) = [];

VocalData.FRA_area = sum(IsRF_removeNans)/numel(IsRF_removeNans);

%% Define vocalizations and valences
vocalizations = StimInfo.Combined(:,1);
valences = [1 1 1 1 1 1 -1 -1 -1 0 1 1 1 1 1 1 -1 -1 -1 0];% 1 is positive, -1 is negative and 0 is neutral
frequencies = [1 1 1 1 1 -1 -1 -1 -1 0 1 1 1 1 1 -1 -1 -1 -1 0];%1 is high, -1 is low and 0 is mid

%Which PeakData to use for Response value?
%--> Use AUC, but make suppressed AUCs negative for proper inhibitory/excitatory RF detection
%--> If ALL significant responses are suppressed, flip the signs so that we can measure suppressed RFs
RF_type = determine_RF_response_type(type, isResponsive);

%Generate RF from AUC, assuming excitatory or mixed RF for now
response = nan(size(PeakData.Peak_AUC));
for i = 1:length(response)
    response_type = PeakData.ResponseType(i);
    if strcmp(response_type, 'suppressed')
        response(i) = -PeakData.Trough_AUC(i);
    else
        response(i) = PeakData.Peak_AUC(i);
    end
end

%If all significant responses are suppressed, flip the sign of RF
response_orig = response; %Save for storing in VocalData without sign flip
if strcmp(RF_type, 'inhibitory')
    response = -response;
end

VocalData.Response = single(response_orig');

%% d prime sensitivity index (from Romero, Polley et al Cerebral Cortex 2019)
%This is the only part of the code where we keep the responses in a matrix

%Get RF mapping from StimInfo
RF_Map = StimInfo.RF_Map;
nanmat = nan(max(RF_Map));
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

%Assign RF values to RF map
[RF, thresholded_RF] = deal(nanmat);
RF(RF_ind) = response_orig;
temp_response_thresholded = response_orig;
temp_response_thresholded(~isResponsive) = nan;
thresholded_RF(RF_ind) = temp_response_thresholded;
[M,I] = max(thresholded_RF(:));
[best_row, best_col] = ind2sub(size(RF),I);

%% Eliminate stim that aren't represented (e.g. due to loco or changes in stim protocol) or have NaN response
missing_ind = strcmp(type, 'undetermined');
nan_ind = isnan(response);
remove_ind = (missing_ind + nan_ind) > 0;
isResponsive(remove_ind) = [];
if ~isempty(remove_ind)
    valences(remove_ind) = [];
    vocalizations(remove_ind) = [];
    response(remove_ind) = [];
end

%% Find preferred vocalization [Two ways]
% Use vocalization with greatest response

%Raw values
if ~all(isnan(response))
    [~, max_ind] = max(response);

    VocalData.Raw_BestVocal = single(vocalizations(max_ind));
else
    raw_best_ind = [];
end

%Thresholded values
response_thresholded = response;
response_thresholded(~isResponsive) = nan;

if ~all(isnan(response_thresholded))
    [~, max_ind] = max(response_thresholded);

    VocalData.BestVocal = single(vocalizations(max_ind));
else
    best_ind = []; %For fitting
end

%% Valence Index
% Resp_Pos - Resp_Neg / Resp_Pos + Resp_Neg
% Asymmetry index betwen -1 and 1 indicating whether cell responds more to
% pos or negative valence vocalizations
% If value is closer to -1, cell prefers neg valence, if value is closer to
% +1, cell prefers pos valence

pos_ind = valences == 1;
neg_ind = valences == -1;
VocalData.VI = tak_compute_asymmetry_index(response, pos_ind, neg_ind);

% %% Test Valence Index (Two ways)
% 
% if ~isnan(VocalData.VI) %VI will be NaN if there weren't both pos and neg values
% 
%     VI_test = struct;
% 
%     %% FIRST WAY: Create shuffled distribution by shuffling stim IDs
%     %Generate list of shuffled IDs
%     N1 = size(valences,2); % number of stim
% 
%     shuffled_responses = nan(nShuffles,N1);
%     for u = 1:nShuffles
%         permIDs = randperm(N1);
% 
%         %Make sure the original order isn't included
%         while isequal(permIDs,1:N1)
%            % print(permIDs);
%             permIDs = randperm(N1);
%         end
% 
%         %Threshold negative values to be 0
%         tempResponse = response(permIDs);
%         tempResponse(tempResponse < 0) = 0;
% 
%         shuffled_responses(u,:) = tempResponse;
%     end
% 
%     %Compute VI on all shuffled IDs
%     VI_Shuffled = tak_compute_asymmetry_index(shuffled_responses, pos_ind, neg_ind);
% 
%     %Perform two-tailed z test
%     [VI_test.Shuffle.Test, VI_test.Shuffle.P, VI_test.Shuffle.CI, VI_test.Shuffle.R] = ztest(VocalData.VI, mean(VI_Shuffled), std(VI_Shuffled), 'Alpha', alpha);
% 
%     %Store p value
%     VocalData.VI_p = single(VI_test.Shuffle.P);
% 
%     %% SECOND WAY: Compare bootstrapped MI distribution to 0
% 
%     %Calculate the average of r subsamples nShuffle number of times with replacement
%     N1 = size(valences,2); % number of stim
% 
%     if N1 > 2
%         r = N1-1; %Subsample (keep as much of the data as possible to make sure we have both + and - sweeps)
% 
%         shuffle_inds = sort(repmat((1:nShuffles)',[r,1]));
%         bootstrap_ind = randsample(1:N1,r*nShuffles,'true')';
% 
%         VI_bootstrap = nan(nShuffles,1);
%         for u = 1:nShuffles
%             %Define sample
%             sample_ind = bootstrap_ind(shuffle_inds == u);
%             sample_valences = valences(sample_ind);
% 
%             %Make sure we don't have all downs or all ups
%             if all(sample_valences == 0) || all(sample_valences == 1)
%                 hasPosAndNeg = 0;
%                 while ~hasPosAndNeg
%                     sample_ind = randsample(1:N1,r,'true')';
%                     sample_valences = valences(sample_ind);
%                     if ~(all(sample_valences == 0) || all(sample_valences == 1))
%                         hasPosAndNeg = 1;
%                     end
%                 end
%             end
% 
%             %Compute VI
%             sample_response = response(sample_ind);
%             sample_pos_ind = sample_valences == 1;
%             sample_neg_ind = sample_valences == -1;
%             VI_bootstrap(u) = tak_compute_asymmetry_index(sample_response, sample_pos_ind, sample_neg_ind);
% 
% %             if isnan(VI_bootstrap(u))
% %                 if sample_pos_ind == 0 && sample_neg_ind == 0
% %                     %If value is NaN because numerator and denominator are 0, set MI to 0
% %                     VI_bootstrap(u) = 0;
% %                 else
% %                     error('Correct NaNs')
% %                 end
% %             end
%         end
% 
%         %Perform Z-test asking whether our bootstrapped MI distribution is significantly different from 0
%         if VocalData.VI > 0; tail = 'left'; elseif VocalData.VI < 0; tail = 'right'; else; tail = 'both'; end %Determine whether to do left or right-tailed test
%         [VI_test.Bootstrap.Test, VI_test.Bootstrap.P, VI_test.Bootstrap.CI, VI_test.Bootstrap.R] = ztest(0, mean(VI_bootstrap), std(VI_bootstrap), 'Alpha', alpha, 'tail', tail);
% 
%         %Store p value
%         VocalData.VI_p_bs = single(VI_test.Bootstrap.P);
%     else
%         %Cannot compute
%         VI_test.Bootstrap.Test = NaN;
%         VI_test.Bootstrap.P = NaN;
%         VI_test.Bootstrap.CI = [NaN, NaN];
%         VI_test.Bootstrap.R = NaN;
%         VocalData.VI_p_bs = NaN;
%     end
% end


%% Spectral Index
% Resp_Pos - Resp_Neg / Resp_Pos + Resp_Neg
% Asymmetry index betwen -1 and 1 indicating whether cell responds more to
% hi or low frequency vocalizations
% If value is closer to -1, cell prefers low freq, if value is closer to
% +1, cell prefers high frequency

pos_ind = frequencies == 1;
neg_ind = frequencies == -1;
VocalData.FI = tak_compute_asymmetry_index(response, pos_ind, neg_ind);
        
% %% Test Spectral Index (Two ways)
% 
% if ~isnan(VocalData.FI) %VI will be NaN if there weren't both pos and neg values
%     
%     FI_test = struct;
% 
%     %% FIRST WAY: Create shuffled distribution by shuffling stim IDs
%     %Generate list of shuffled IDs
%     N1 = size(frequencies,2); % number of stim
% 
%     shuffled_responses = nan(nShuffles,N1);
%     for u = 1:nShuffles
%         permIDs = randperm(N1);
% 
%         %Make sure the original order isn't included
%         while isequal(permIDs,1:N1)
%             permIDs = randperm(N1);
%         end
% 
%         %Threshold negative values to be 0
%         tempResponse = response(permIDs);
%         tempResponse(tempResponse < 0) = 0;
%         
%         shuffled_responses(u,:) = tempResponse;
%     end
% 
%     %Compute FI on all shuffled IDs
%     FI_Shuffled = tak_compute_asymmetry_index(shuffled_responses, pos_ind, neg_ind);
% 
%     %Perform two-tailed z test
%     [FI_test.Shuffle.Test, FI_test.Shuffle.P, FI_test.Shuffle.CI, FI_test.Shuffle.R] = ztest(VocalData.FI, mean(FI_Shuffled), std(FI_Shuffled), 'Alpha', alpha);
%     
%     %Store p value
%     VocalData.FI_p = single(FI_test.Shuffle.P);
%     
%     %% SECOND WAY: Compare bootstrapped MI distribution to 0
% 
%     %Calculate the average of r subsamples nShuffle number of times with replacement
%     N1 = size(frequencies,2); % number of stim
%     
%     if N1 > 2
%         r = N1-1; %Subsample (keep as much of the data as possible to make sure we have both + and - sweeps)
% 
%         shuffle_inds = sort(repmat((1:nShuffles)',[r,1]));
%         bootstrap_ind = randsample(1:N1,r*nShuffles,'true')';
% 
%         FI_bootstrap = nan(nShuffles,1);
%         for u = 1:nShuffles
%             %Define sample
%             sample_ind = bootstrap_ind(shuffle_inds == u);
%             sample_frequencies = frequencies(sample_ind);
% 
%             %Make sure we don't have all downs or all ups
%             if all(sample_frequencies == 0) || all(sample_frequencies == 1)
%                 hasPosAndNeg = 0;
%                 while ~hasPosAndNeg
%                     sample_ind = randsample(1:N1,r,'true')';
%                     sample_frequencies = frequencies(sample_ind);
%                     if ~(all(sample_frequencies == 0) || all(sample_frequencies == 1))
%                         hasPosAndNeg = 1;
%                     end
%                 end
%             end
% 
%             %Compute VI
%             sample_response = response(sample_ind);
%             sample_pos_ind = sample_frequencies == 1;
%             sample_neg_ind = sample_frequencies == -1;
%             FI_bootstrap(u) = tak_compute_asymmetry_index(sample_response, sample_pos_ind, sample_neg_ind);
% % 
% %             if isnan(FI_bootstrap(u))
% %                 if sample_Resp_pos == 0 && sample_Resp_neg == 0
% %                     %If value is NaN because numerator and denominator are 0, set MI to 0
% %                     FI_bootstrap(u) = 0;
% %                 else
% %                     error('Correct NaNs')
% %                 end
% %             end
%         end
% 
%         %Perform Z-test asking whether our bootstrapped MI distribution is significantly different from 0
%         if VocalData.FI > 0; tail = 'left'; elseif VocalData.FI < 0; tail = 'right'; else; tail = 'both'; end %Determine whether to do left or right-tailed test
%         [FI_test.Bootstrap.Test, FI_test.Bootstrap.P, FI_test.Bootstrap.CI, FI_test.Bootstrap.R] = ztest(0, mean(FI_bootstrap), std(FI_bootstrap), 'Alpha', alpha, 'tail', tail);
% 
%         %Store p value
%         VocalData.FI_p_bs = single(FI_test.Bootstrap.P);
%     else
%         %Cannot compute
%         FI_test.Bootstrap.Test = NaN;
%         FI_test.Bootstrap.P = NaN;
%         FI_test.Bootstrap.CI = [NaN, NaN];
%         FI_test.Bootstrap.R = NaN;
%         VocalData.FI_p_bs = NaN;
%     end
% end
% %% Plot figure
% 
% if plot_figures
%     fig1 = figure;
% 
%     subplot(3,3,1); hold on
%     plot(vocalizations, response, 'Linewidth', 2);
%     xlim([vocalizations(1) vocalizations(end)])
%     set(gca, 'XTick', vocalizations)
%     set(gca, 'XTickLabel', num2str(vocalizations))
% end
% xlabel('Vocalizations')
% ylabel('Response');
% title(['Best Vocalization: ' num2str(VocalData.BestVocal)])
% 
% %Plot VI tests
% if ~isnan(VocalData.VI)
%     subplot(3,3,4); hold on
%     histogram(VI_Shuffled, 'facecolor', 'g', 'facealpha', 0.2)
%     vline(VocalData.VI)
%     ylabel('Shuffles')
%     title(['VI (pos vs. neg): ' num2str(VocalData.VI)])
%     if VI_test.Shuffle.Test; textcolour = ' {\color{green}p=} '; else; textcolour = ' {\color{black}p=} '; end
%     subtitle(strcat('Z=', sprintf('%.2f', VI_test.Shuffle.R), textcolour, sprintf('%.2f', VI_test.Shuffle.P)))
% 
%     subplot(3,3,7); hold on
%     histogram(VI_bootstrap, 'facecolor', 'b', 'facealpha', 0.2)
%     vline(VocalData.VI)
%     vline(0, 'k')
%     ylabel('Bootstraps')
%     xlabel('Valence Index')
%     if VI_test.Bootstrap.Test; textcolour = ' {\color{green}p=} '; else; textcolour = ' {\color{black}p=} '; end
%     subtitle(strcat('Z=', sprintf('%.2f', VI_test.Bootstrap.R), textcolour, sprintf('%.2f', VI_test.Bootstrap.P)))
% else
%     subplot(3,3,4);
%     title('VI could not be computed')
% end
% 
% 
% %Plot FI tests
% if ~isnan(VocalData.FI)
%     subplot(3,3,5); hold on
%     histogram(FI_Shuffled, 'facecolor', 'g', 'facealpha', 0.2)
%     vline(VocalData.FI)
%     ylabel('Shuffles')
%     title(['FI (low vs. high): ' num2str(VocalData.FI)])
%     if FI_test.Shuffle.Test; textcolour = ' {\color{green}p=} '; else; textcolour = ' {\color{black}p=} '; end
%     subtitle(strcat('Z=', sprintf('%.2f', FI_test.Shuffle.R), textcolour, sprintf('%.2f', FI_test.Shuffle.P)))
% 
%     subplot(3,3,8); hold on
%     histogram(FI_bootstrap, 'facecolor', 'b', 'facealpha', 0.2)
%     vline(VocalData.FI)
%     vline(0, 'k')
%     ylabel('Bootstraps')
%     xlabel('Frequency Index')
%     if FI_test.Bootstrap.Test; textcolour = ' {\color{green}p=} '; else; textcolour = ' {\color{black}p=} '; end
%     subtitle(strcat('Z=', sprintf('%.2f', FI_test.Bootstrap.R), textcolour, sprintf('%.2f', FI_test.Bootstrap.P)))
% else
%     subplot(3,3,5);
%     title('FI could not be computed')
% end
% 
 end
