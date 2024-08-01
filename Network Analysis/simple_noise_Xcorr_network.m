function [NoiseCorr, blanks]  ...
    = simple_noise_Xcorr_network(ExtractedData, motorType, cell_list, NoiseCorrOps)

% This function determines the stimulus noise correlations (correlations of
% residuals after subtracting each trial from its mean)
% using pairwise correlations (XCorr) between every cell in the block. 
 
% Code adapted from Maryse Thomas, Carolyn Sweeney, Rahul Brito and Anne Takesian
% modified by Anne Takesian in April 2022

% version 3 September 2022: made compatible with simple_Extracted_Data and
% modified based on M. Thomas simple_xcorr code

% Argument(s): 

% use_responsive_cells - only loop through cells with responsive stimuli
% match_RFs - only use stimuli which both neurons respond to
% ExtractedData
% motor Type :All, Loco, none - from Extracted Data
% cell_rows: cells to analyze in Extracted Data 

% Returns:
% Returns NoiseCorr table with:
% noiseCorr - correlations of residuals measured by Pearson's correlation
% signalCorr - correlations of signals measured by Pearson's correlation
%   " " _z = Z-test result (0 or 1) for comparing each of the above values against control distribution
%   " " _p = p-value of the test
%   " " _r = z-score from the test

% Blanks ('yes blanks' or 'no blanks') - were there blanks to use for
% control distributions

%% Pull out data from Extracted Data for Specific Block
% Define red cells 
% if using red cells only, make list of red cells
if NoiseCorrOps.only_red_cells
    red_cell_list = ExtractedData.CellList.RedCell;
    cell_list = intersect(find(red_cell_list == 1),cell_list); 
end

RFdata =  ExtractedData.CellDataByStim.([motorType])(cell_list);

% find RF trial index that represents conditions that arent presented to
% the cell:
for j = 1:length(RFdata(1).StimTraces)
    if isempty(RFdata(1).StimTraces(j))
        tn(j) = 0;
    else
        tn(j) = 1;
    end
end
trialnan = find(tn ==0);

%% Control mode for signal correlations
% 0 = use blanks, 1 = use shifted data, 2 = do both
if NoiseCorrOps.use_correlation_shift == 0 %|| use_correlation_shift == 2
    %If no blank trials or less blank trials than real trials, default to shifted data instead
    if isempty(RFdata(1).BlankTraces) % the first cell will have same blanks as others in set so r_cellIDX(1)
        NoiseCorrOps.use_correlation_shift = 1;
        blanks = 'no blanks'; % report if no blanks
    else
        blanks = 'yes blanks';
    end
end

mode1 = NoiseCorrOps.use_correlation_shift == 0; %|| use_correlation_shift == 2; % using blank trials
mode2 = NoiseCorrOps.use_correlation_shift == 1; %|| use_correlation_shift == 2; % using shifted traces

%% Initialize correlation matrices
Signal_Correlation = nan(length(cell_list)); % matrix size of all cells in block
Noise_Correlation = nan(length(cell_list));
Signal_Z_test = nan(length(cell_list));
Signal_Z_p = nan(length(cell_list));
Signal_Z_z = nan(length(cell_list));
Noise_Z_test = nan(length(cell_list));
Noise_Z_p = nan(length(cell_list));
Noise_Z_z = nan(length(cell_list));

%% Figure out how many unique pairs there are
pairs = [];
count = 1;
for n1 = 1:size(cell_list,1)
    for n2 = 1:size(cell_list,1)
        pairs(count,1:2) = [n1, n2];
        count = count + 1;
    end
end
unique_pairs = unique(sort(pairs')', 'rows');
unique_pairs(unique_pairs(:,1) == unique_pairs(:,2),:) = [];
nPairs = size(unique_pairs,1); %equivalent to nchoosek(size(traces,1),2)

unique_rows = nan(nPairs,1);
for i = 1:nPairs
    unique_rows(i) = find(sum((pairs - unique_pairs(i,:)) == 0, 2) == 2);
end

%% Make template tables that will be used to store NoiseCorr results 
if size(cell_list,1) > 1
    npairs = nchoosek(size(cell_list,1),2); % number of pairs
else
    npairs = 1;
end
template = table;
template.noiseCorr = nan(npairs,1);  
template.signalCorr = nan(npairs,1); 
template.noise_z = nan(npairs,1);
template.noise_p = nan(npairs,1);
template.noise_r = nan(npairs,1);
template.signal_z = nan(npairs,1);
template.signal_p = nan(npairs,1);
template.signal_r = nan(npairs,1);

%Return empty template if only 1 cell found (can't compute correlations)
if size(cell_list,1) <= 1
    NoiseCorr = template;
    return;
end

%% Loop through cell pairs to measure signal and noise correlations
for c1 = 1:length(cell_list)-1 % loop through each cell (cell #1)
    for c2 = c1+1:length(cell_list) % loop through each neighboring cell for pairwise comparison (cell #2)
        cell1 = c1;
        cell2 = c2;

        cell1_RF = RFdata(c1).RF;
        cell1_RF(trialnan)=[]; % eliminate rows with nan

        cell2_RF = RFdata(c2).RF; % cell 2 RF
        cell2_RF(trialnan)=[]; % eliminate rows with nan

        % find matching responses within RFs of both neurons
        if NoiseCorrOps.match_RFs == 1 % if only looking at responses to stimuli that both neurons are responsive to...
            A =  intersect(find(cell1_RF),find(cell2_RF)); % matching RFs
            [ii] = A;
        else
            A = find((cell1_RF+cell2_RF)>=1);
            [ii] = A; % use stimuli which either cell 1 or cell 2 are responsive
        end

        %index = 1; % for testing Pearson's function
        %Find rows corresponding to all matched stimuli
        if ~isempty(ii)
            cell_1_residual = [];
            cell_2_residual = [];
            cell_1_signal = [];
            cell_2_signal = [];
            XCorr_S_Total = []; 

            if NoiseCorrOps.plot_traces; sz = ceil(sqrt(length(ii))); fig1 = figure;  fig2 = figure; fig3 = figure; end % if plot traces for testing

            for n=1:length(ii) % loop through stimuli that elicit responses in both neurons or all stimuli if not matching RFs
                row = ii(n);
                if size(RFdata(c1).StimTraces{row},1)>NoiseCorrOps.trial_lim % have to have enough trials for noise analysis
                    F1 = RFdata(c1).StimTraces{row};
                    F2 = RFdata(c2).StimTraces{row};
                    trial_n = size(F1,1);

                    cell_1_mean = RFdata(c1).StimTracesAveraged(row,:); % mean response of cell 1 to the given stimulus
                    cell_2_mean = RFdata(c2).StimTracesAveraged(row,:); % mean response of cell 2 to the given stimulus

                    cell_1_signal = [cell_1_signal; F1]; % building vector of signals for all responses to all stimuli of cell 1
                    cell_2_signal = [cell_2_signal; F2]; % building vector of signals for all responses to all stimuli of cell 2
                    cell_1_residual = [cell_1_residual; F1-cell_1_mean]; % residuals for all responses to stimulus of cell 1
                    cell_2_residual = [cell_2_residual; F2-cell_2_mean]; % residuals for all responses to stimulus of cell 2

                    % Signal Correlation - XCorr correlation across all
                    % responses to the same-stimuli but unmatched
                    % (non-simultaneous) trials
                    [XCorr_A,lags] = xcorr([F1; F2]', NoiseCorrOps.maxlag*NoiseCorrOps.fs, 'coef'); % XCorr across all trials from 2 cells
                    max_XCorr = max(XCorr_A,[],1)';

                    % index of unmatched pairs in cross correlation matrix
                    XCorr_matrix = zeros(trial_n*2,trial_n*2);
                    XCorr_matrix(1:trial_n,trial_n+1:2*trial_n) = 1; % set right quadrant to 1
                    XCorr_matrix(trial_n*2*trial_n+1:trial_n*2+1:end) = 0; % set diagonal elements (matched stimuli) to 0
                    XCorr_ind = find(XCorr_matrix);
                    XCorr_S = max_XCorr(XCorr_ind);
                    XCorr_S_Total = [XCorr_S_Total; XCorr_S];

                    %loop XCorr (slow way) to double-check indexing
                    count = 1;
                    for i=1:size(F2,1)
                        trace1 = F2(i,:);
                        for j = 1:size(F1,1)
                            if i ~=j
                            trace2 = F1(j,:);
                              [CorrX, lag] = xcorr(trace1,trace2, NoiseCorrOps.maxlag*NoiseCorrOps.fs, 'coeff');
                              XCorr_max_test(count,:) = max(CorrX,[],2);
                              count = count +1;
                            end
                        end
                    end


                    if NoiseCorrOps.plot_traces % if plot traces of comparisons
                        [stim_units, V1_label, V2_label] = assign_stim(block.setup.stim_protocol);
                        figure(fig1);
                        for s = 1:size(F1,1) % plot residuals of each cell
                            subplot(length(i),size(F1,1),(n-1)*size(F1,1)+s);
                            plot((F1(s,:)-cell_1_mean)','b'); hold on; plot((F2(s,:)-cell_2_mean)','r');
                            a = ['Residuals of Cells ' num2str(block.cell_number(cell1))  ' & ' num2str(block.cell_number(cell2))];
                            b = [num2str(unique_stim_v1(i(n))) ' ' stim_units{1}  ' ' num2str(unique_stim_v2(ii(n))) ' ' stim_units{2}];
                            title(a);
                            subtitle(b);
                        end

                        figure(fig2); % plot traces for each cell
                        for s = 1:size(F1,1)
                            subplot(length(i),size(F1,1),(n-1)*size(F1,1)+s);
                            plot(F1(s,:)','b'); hold on; plot(F2(s,:)','r');
                            a = ['Traces of Cells ' num2str(block.cell_number(cell1))  ' & ' num2str(block.cell_number(cell2))];
                            b = [num2str(unique_stim_v1(i(n))) ' ' stim_units{1}  ' ' num2str(unique_stim_v2(ii(n))) ' ' stim_units{2}];
                            title(a);
                            subtitle(b);
                        end

                        figure(fig3);
                        subplot(sz,sz,n);
                        plot(F1','b'); hold on; plot(F2','r');
                        title(['Traces of Cell ' num2str(block.cell_number(cell1)) ' & ' num2str(block.cell_number(cell2)) ' ' num2str(unique_stim_v1(i(n))) ' ' V1_label ' ' num2str(unique_stim_v2(ii(n))) ' ' V2_label]);
                    end
                end

                % can use the following to observe or test correlations
                %                     for r=1:5
                %                         x = corrcoef(cell_1_signal(r,:), cell_2_signal(r,:))
                %                         Pearsons_signal(index) = x(1,2);
                %                         index = index +1
                %                     end
                %                     assignin('base','Pearsons_signal',Pearsons_signal)
            end

            % Perform Noise and Signal Correlations across All Trials
            row_num = size(cell_1_signal,1); % number of trials compared (mixed stimuli)
            if row_num>NoiseCorrOps.trial_lim % if total # of comparisons greater than some limit

                % Noise Correlation - XCorr across all residuals from 2
                % cells, pull out the relevant comparisons, sum and assign to
                % correlation matrix
                [XCorr_A,lags] = xcorr([cell_1_residual; cell_2_residual]', NoiseCorrOps.maxlag*NoiseCorrOps.fs, 'coef'); % XCorr across all trials from 2 cells
                max_XCorr = max(XCorr_A,[],1)';

                % index of unmatched pairs in cross correlation matrix
                XCorr_matrix_un = zeros(row_num*2,row_num*2);
                XCorr_matrix_un(1:row_num,row_num+1:2*row_num) = 1; % set right quadrant to 1
                XCorr_matrix_un(row_num*2*row_num+1:row_num*2+1:end) = 0; % set diagonal elements (matched stimuli) to 0
                XCorr_unmatched_ind = find(XCorr_matrix_un);
                
                % index of matched pairs in cross correlation matrix
                XCorr_matrix_matched = zeros(row_num*2,row_num*2);
                XCorr_matrix_matched(row_num*2*row_num+1:row_num*2+1:end) = 1; % set diagonal elements (matched stimuli) to 1
                XCorr_matched_ind = find(XCorr_matrix_matched);
                
                XCorr_unmatched = max_XCorr(XCorr_unmatched_ind); % unmatched stimuli
                XCorr_N1_2 = max_XCorr(XCorr_matched_ind); % matched stimuli only
                Noise_Correlation(cell1,cell2) = mean(mean(XCorr_N1_2));
                Noise_Correlation(cell2,cell1) = mean(mean(XCorr_N1_2));

                % Noise XCorrelation Validation: Perform shuffles on unmatched trials for control
                % distribution
                for u = 1:NoiseCorrOps.shuffles
                    % Pull random correlations from unmatched pairs, n =
                    % #stimuli pairs compared above (row_num)
                    rand_trial_indices = randsample(length(XCorr_unmatched),row_num,true); % n rand numbers from 1:length of non-matched correlations list w replacement
                    rand_trial_noise_comparisons = XCorr_unmatched(rand_trial_indices); % pull out rand noise correlations from unmatched trials
                    noise_correlations_unmatched(u) = mean(rand_trial_noise_comparisons); % unmatched noise correlations for control distribution
                end

                % perform z test to determine if noise correlation is from same
                % distribution as shuffled trials
                [z_test_noise, z_P_noise, z_CI_noise,z_R_noise]...
                    = ztest(Noise_Correlation(cell1,cell2), mean(noise_correlations_unmatched), std(noise_correlations_unmatched),"Alpha",NoiseCorrOps.p_value);
                Noise_Z_test(cell1,cell2) = z_test_noise; Noise_Z_test(cell2,cell1) = z_test_noise;
                Noise_Z_p(cell1,cell2) = z_P_noise; Noise_Z_p(cell2,cell1) = z_P_noise;
                Noise_Z_z(cell1,cell2) = z_R_noise; Noise_Z_z(cell2,cell1) = z_R_noise;

                % Signal XCorrelations
                Signal_Correlation(cell1,cell2) = mean(XCorr_S_Total); % signal correlation matrix
                Signal_Correlation(cell2,cell1) = mean(XCorr_S_Total);
                n_sig_corr = length(XCorr_S_Total);
                

                %% Signal Correlation Validation: Perform shuffles on unmatched blanks or time-shifted data for controls

                %Blank trials
                blank_cell1 = RFdata(c1).BlankTraces;
                blank_cell2 = RFdata(c2).BlankTraces;
                NB = size(blank_cell1,1);

                %check if enough blank trials to sample
                if ~isempty(blank_cell1) && NB>NoiseCorrOps.trial_lim
                    mode1 = 1;
                    mode2 = 0;
                else
                    mode1 = 0;
                    mode2 = 1;
                end

                if mode1 % if using blanks as signal control

                    % Compute XCorr of blanks for all trials except
                    % matched
                    [XCorr_blanks,lags] = xcorr([blank_cell1; blank_cell2]', NoiseCorrOps.maxlag*NoiseCorrOps.fs, 'coef'); % XCorr across all trials from 2 cells
                    max_XCorr_blanks = max(XCorr_blanks,[],1)';

                    % index of unmatched pairs in cross correlation matrix
                    XCorr_matrix = zeros(NB*2,NB*2);
                    XCorr_matrix(1:NB,NB+1:2*NB) = 1; % set right quadrant to 1
                    XCorr_matrix(NB*2*NB+1:NB*2+1:end) = 0; % set diagonal elements (matched stimuli) to 0
                    XCorr_ind = find(XCorr_matrix);
                    Blanks = max_XCorr_blanks(XCorr_ind);

                    % Boostrap blank or signal corrrelations to compare
                    % distributions of same number
                    if length(Blanks)>n_sig_corr % if bootstrapping blanks
                        signal_correlation_blanks = nan(NoiseCorrOps.shuffles,1);
                        for u = 1:NoiseCorrOps.shuffles
                            % Pull number of correlations from matched
                            % blank pairs = number of stim pairs
                            rand_blank_indices= randsample(length(Blanks),n_sig_corr, true); % n rand numbers from 1:length of matched blank trials with replacement
                            rand_trial_blanks = Blanks(rand_blank_indices); % pull out signal correlations from matched blank trials
                            signal_correlation_blanks(u) = mean(rand_trial_blanks);
                        end

                        % perform z test to determine if signal correlations is from same
                        % distribution as blank trials
                        [z_test_signal,z_P_signal,z_CI_signal,z_R_signal]...
                            = ztest(Signal_Correlation(cell1,cell2), mean(signal_correlation_blanks), std(signal_correlation_blanks),"Alpha",NoiseCorrOps.p_value);
                        Signal_Z_test(cell1,cell2) = z_test_signal; Signal_Z_test(cell2,cell1) = z_test_signal;
                        Signal_Z_p(cell1,cell2) = z_P_signal; Signal_Z_p(cell2,cell1) = z_P_signal;
                        Signal_Z_z(cell1,cell2) = z_R_signal; Signal_Z_z(cell2,cell1) = z_R_signal;

                    else % if bootstrapping signal
                        correlation_signal = nan(NoiseCorrOps.shuffles,1);
                        for u = 1:NoiseCorrOps.shuffles
                            % Pull number of blank correlations from
                            % matched signal pairs = number of blank pairs
                            rand_sig_indices= randsample(n_sig_corr,length(Blanks),true); % n rand numbers from 1:length of matched blank trials with replacement
                            rand_trial_sig = XCorr_S_Total(rand_sig_indices); % pull out signal correlations
                            signal_correlations(u) = mean(rand_trial_sig);
                        end

                        % perform z test to determine if signal correlations is from same
                        % distribution as blank trials
                        [z_test_signal,z_P_signal,z_CI_signal,z_R_signal]...
                            = ztest(mean(Blanks), mean(signal_correlations), std(signal_correlations),"Alpha",NoiseCorrOps.p_value);
                        Signal_Z_test(cell1,cell2) = z_test_signal; Signal_Z_test(cell2,cell1) = z_test_signal;
                        Signal_Z_p(cell1,cell2) = z_P_signal; Signal_Z_p(cell2,cell1) = z_P_signal;
                        Signal_Z_z(cell1,cell2) = z_R_signal; Signal_Z_z(cell2,cell1) = z_R_signal;
                    end

                end

                %% Measure signal correlation of shifted traces as control and determine Z-score
                N2 = size(cell_1_signal,2); % number of frames in each trial (same for both cells)
                if mode2  %if isempty(r_blanks)
                    %Pull out comparisons from shifted trials to compute the mean Pearsons & do this shuffles number of times
                    signal_correlation_shifted = nan(NoiseCorrOps.shuffles,1);

                    for u = 1:NoiseCorrOps.shuffles
                        % Pull number of correlations from unmatched shifted pairs =
                        % number of stimuli pairs compared above (row_num)
                        rand_shifting= randsample(1:N2,row_num,true); % random frame for which to start shift
                        shifted_c1_trials = cell2mat(arrayfun(@(x) circshift(cell_1_signal(x,:),[1 rand_shifting(x)]),(1:numel(rand_shifting))','un',0));

                        [XCorr_shifted,lags] = xcorr([shifted_c1_trials;cell_2_signal ]', NoiseCorrOps.maxlag*NoiseCorrOps.fs, 'coef'); % XCorr across all trials from 2 cells
                        max_XCorr_shifted = max(XCorr_shifted,[],1)';

                        % index of unmatched pairs in cross correlation matrix
                        XCorr_matrix = zeros(row_num*2,row_num*2);
                        XCorr_matrix(1:row_num,row_num+1:2*row_num) = 1; % set right quadrant to 1
                        XCorr_matrix(row_num*2*row_num+1:row_num*2+1:end) = 0; % set diagonal elements (matched stimuli) to 0
                        XCorr_ind = find(XCorr_matrix);
                        Shifted = max_XCorr_shifted(XCorr_ind);
                        signal_correlation_shifted(u) = mean(Shifted);
                    end

                    % perform z test to determine if signal correlations is from same
                    % distribution as shuffled trials
                    [z_test_signal,z_P_signal,z_CI_signal,z_R_signal]...
                        = ztest(Signal_Correlation(cell1,cell2), mean(signal_correlation_shifted), std(signal_correlation_shifted),"Alpha", NoiseCorrOps.p_value);
                    Signal_Z_test(cell1,cell2) = z_test_signal; Signal_Z_test(cell2,cell1) = z_test_signal;
                    Signal_Z_p(cell1,cell2) = z_P_signal; Signal_Z_p(cell2,cell1) = z_P_signal;
                    Signal_Z_z(cell1,cell2) = z_R_signal; Signal_Z_z(cell2,cell1) = z_R_signal;
                end

                if NoiseCorrOps.plot_histograms

                    % histogram of signal correlations
                    figure;

                    if mode1
                        if length(Blanks)>n_sig_corr % if shuffled blanks
                            min_g = min(min(signal_correlation_blanks), Signal_Correlation(cell1,cell2))-0.1;
                            max_g = max(max(signal_correlation_blanks), Signal_Correlation(cell1,cell2))+0.1;

                            fig4 = histogram(signal_correlation_blanks,min_g:.01:max_g,'facecolor','b','facealpha',0.2,'edgecolor','none');
                            [maxcount, whichbin] = max(fig4.Values); hold on;
                            plot([Signal_Correlation(cell1,cell2) Signal_Correlation(cell1,cell2)],[0 maxcount], '--r');
                        else % if shuffled signal
                            min_g = min(min(signal_correlations), mean(Blanks))-0.1;
                            max_g = max(max(signal_correlations), mean(Blanks))+0.1;

                            fig4 = histogram(signal_correlations,min_g:.01:max_g,'facecolor','b','facealpha',0.2,'edgecolor','none');
                            [maxcount, whichbin] = max(fig4.Values); hold on;
                            plot([mean(Blanks) mean(Blanks)],[0 maxcount], '--r');
                        end
                        title(['Signal Correlation of Cell ' num2str(block.cell_number(cell1))  ' and Cell ' num2str(block.cell_number(cell2)) ]);
                    end

                    if mode2 % if using shifted control distribution
                        min_g = min(min(signal_correlation_shifted), Signal_Correlation(cell1,cell2))-0.1;
                        max_g = max(max(signal_correlation_shifted),Signal_Correlation(cell1,cell2))+0.1;

                        fig4 = histogram(signal_correlation_shifted,min_g:.01:max_g,'facecolor','g','facealpha',0.2,'edgecolor','none');
                        [maxcount, whichbin] = max(fig4.Values); hold on;
                        plot([Signal_Correlation(cell1,cell2) Signal_Correlation(cell1,cell2)],[0 maxcount], '--r');
                        title(['Signal Correlation of Cell ' num2str(block.cell_number(cell1))  ' and Cell ' num2str(block.cell_number(cell2)) ]);
                    end

                    Signal_Z_test(cell1,cell2)
                    if Signal_Z_test(cell1,cell2) ==1
                        subtitle(strcat('Corr= ', sprintf('%.2f',Signal_Correlation(cell1,cell2)),...
                            ' {\color{green}Zs=} ', sprintf('%.2f', z_R_noise)), 'FontSize',10);
                    else
                        subtitle(strcat('Corr= ', sprintf('%.2f',Signal_Correlation(cell1,cell2)),...
                            ' {\color{black}Zs=} ', sprintf('%.2f', z_R_noise)), 'FontSize',10);
                    end

                    % histogram of noise correlations
                    figure;
                    min_g = min(min(noise_correlations_unmatched), Noise_Correlation(cell1,cell2))-0.1;
                    max_g = max(max(noise_correlations_unmatched), Noise_Correlation(cell1,cell2))+0.1;

                    fig5 = histogram(noise_correlations_unmatched,min_g:.01:max_g,'facecolor','b','facealpha',0.2,'edgecolor','none');
                    [maxcount, whichbin] = max(fig4.Values); hold on;
                    plot([Noise_Correlation(cell1,cell2) Noise_Correlation(cell1,cell2)],[0 maxcount], '--r');

                    Noise_Z_test(cell1,cell2)
                    if Noise_Z_test(cell1,cell2) ==1
                        subtitle(strcat('Corr= ', sprintf('%.2f',Noise_Correlation(cell1,cell2)),...
                            ' {\color{green}Zs=} ', sprintf('%.2f', z_R_noise)), 'FontSize',10);
                    else
                        subtitle(strcat('Corr= ', sprintf('%.2f',Noise_Correlation(cell1,cell2)),...
                            ' {\color{black}Zs=} ', sprintf('%.2f', z_R_noise)), 'FontSize',10);
                    end
                    title(['Noise Correlation of Cell ' num2str(block.cell_number(cell1))  ' and Cell ' num2str(block.cell_number(cell2)) ]);
                    pause;
                end
            else
              Noise_Correlation(cell1,cell2) = nan; Noise_Correlation(cell2,cell1) = nan;
              Noise_Z_test(cell1,cell2) = nan; Noise_Z_test(cell2,cell1) = nan;
              Noise_Z_p(cell1,cell2) = nan; Noise_Z_p(cell2,cell1) = nan; 
              Noise_Z_z(cell1,cell2) = nan; Noise_Z_z(cell2,cell1) = nan;
              Signal_Correlation(cell1,cell2) = nan; Signal_Correlation(cell2,cell1) = nan;
              Signal_Z_test(cell1,cell2) = nan; Signal_Z_test(cell2,cell1) = nan;
              Signal_Z_p(cell1,cell2) = nan; Signal_Z_p(cell2,cell1) = nan; 
              Signal_Z_z(cell1,cell2) = nan; Signal_Z_z(cell2,cell1) = nan;
            end
        
        else
              Noise_Correlation(cell1,cell2) = nan; Noise_Correlation(cell2,cell1) = nan;
              Noise_Z_test(cell1,cell2) = nan; Noise_Z_test(cell2,cell1) = nan;
              Noise_Z_p(cell1,cell2) = nan; Noise_Z_p(cell2,cell1) = nan; 
              Noise_Z_z(cell1,cell2) = nan; Noise_Z_z(cell2,cell1) = nan;
              Signal_Correlation(cell1,cell2) = nan; Signal_Correlation(cell2,cell1) = nan;
              Signal_Z_test(cell1,cell2) = nan; Signal_Z_test(cell2,cell1) = nan;
              Signal_Z_p(cell1,cell2) = nan; Signal_Z_p(cell2,cell1) = nan; 
              Signal_Z_z(cell1,cell2) = nan; Signal_Z_z(cell2,cell1) = nan;

        end
    end
end

%% Build NoiseCorr table

NoiseCorr = template;
variables = {'noise', 'signal'};
variable2 = {'Noise', 'Signal'};

for v = 1:length(variables)

M1 = eval(strcat(variable2{v}, '_Correlation'));
M2 = eval(strcat(variable2{v}, '_Z_test'));
M3 = eval(strcat(variable2{v}, '_Z_p'));
M4 = eval(strcat(variable2{v}, '_Z_z'));

% remove zeros, inf means tested pairs that were non-significant
M1(M1==0) = inf; M2(M2==0) = inf; M3(M3==0) = inf; M4(M4==0) = inf;
M1 = M1(find(tril(M1,-1))); M2 = M2(find(tril(M2,-1))); % Reshape and keep values below the diagonal
M3 = M3(find(tril(M3,-1))); M4 = M4(find(tril(M4,-1)));
M1(M1==inf) = 0; M2(M2==inf) = 0; M3(M3==inf) = 0; M4(M4==inf) = 0;

NoiseCorr.(strcat(variables{v}, 'Corr')) = M1;
NoiseCorr.(strcat(variables{v}, '_z')) = M2;
NoiseCorr.(strcat(variables{v}, '_p')) = M3;
NoiseCorr.(strcat(variables{v}, '_r')) = M4;

end






