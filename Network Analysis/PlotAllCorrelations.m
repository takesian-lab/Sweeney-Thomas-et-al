function PlotAllCorrelations(Correlations_Matrix_All,save_path, Ops)
% Plot output of noise/signal correlation data

% Argument(s):
% Correlations_Table - output from Network Analysis Pipeline
% Save_path  - user-defined save path
% Ops - parameters defined in Network Analysis Pipeline

% Output
% Various graphs of correlation histograms and correlations by distance 
% by correlation type for a given defined cell type

% ToDo - add graphs to compare cell types 

%% Initialize Figures

if contains(Ops.stim,'pontaneous')
    figures = cell(1,12); %Store figures
else
    figures = cell(1,34);
end

%% Sort Data According to Loco, Cell Type and Stim Type

% Sort by Cell Types
cell_type = find(strcmp(Correlations_Matrix_All.CellType,Ops.celltype));

% Sort by Stim Type
stim_type = find(strcmp(Correlations_Matrix_All.StimName,Ops.stim));

% Sort by Mouse_Name
Correlations_Matrix_All.MouseName = categorical(Correlations_Matrix_All.MouseName);
mice = unique(Correlations_Matrix_All.MouseName);
for n = 1:size(mice,1)
    mouse_indices{n} = find(Correlations_Matrix_All.MouseName == mice(n));
end

% Sort by Block_Name
Correlations_Matrix_All.BlockName = categorical(Correlations_Matrix_All.BlockName);
blocks = unique(Correlations_Matrix_All.BlockName);
for n = 1:size(blocks,1)
    block_indices{n} = find(Correlations_Matrix_All.BlockName == blocks(n));
end

% Sort by Locomotion
if any(strcmp(Correlations_Matrix_All.LocoType, 'AllLoc')) || any(strcmp(Correlations_Matrix_All.LocoType, 'NoLoc')) || any(strcmp(Correlations_Matrix_All.LocoType, 'YesLoc'))  
    all_loco_ind = find(strcmp(cellstr(Correlations_Matrix_All.LocoType),'AllLoc'));
    no_loco_ind = find(strcmp(cellstr(Correlations_Matrix_All.LocoType),'NoLoc'));
    yes_loco_ind = find(strcmp(cellstr(Correlations_Matrix_All.LocoType),'YesLoc'));
else
    all_loco_ind = find(strcmp(cellstr(Correlations_Matrix_All.LocoType),'All'));
    no_loco_ind = find(strcmp(cellstr(Correlations_Matrix_All.LocoType),'NotRunning'));
    yes_loco_ind = find(strcmp(cellstr(Correlations_Matrix_All.LocoType),'Running'));
end

% Significant Correlations
TraceMax_Sig = find(Correlations_Matrix_All.TraceMaxZTest == 1);
TraceZero_Sig = find(Correlations_Matrix_All.TraceZeroZTest == 1);

if ~contains(Ops.stim,'pontaneous')
    Noise_Sig = find(Correlations_Matrix_All.NoiseZTest == 1);
    Signal_Sig = find(Correlations_Matrix_All.SignalZTest == 1);
    NoiseX_Sig = find(Correlations_Matrix_All.NoiseXZTest == 1);
    SignalX_Sig = find(Correlations_Matrix_All.SignalXZTest == 1);
end

if isempty(cell_type) && isempty(stim_type) % all correlations
    all = all_loco_ind;
    no = no_loco_ind;
    yes = yes_loco_ind;
elseif ~isempty(cell_type) && isempty(stim_type)  % sort correlations by Cell Type
    all = intersect(cell_type,all_loco_ind);
    no = intersect(cell_type,no_loco_ind);
    yes = intersect(cell_type,yes_loco_ind);
elseif isempty(cell_type) && ~isempty(stim_type) % sort correlations by Stim type
    all = intersect(stim_type,all_loco_ind);
    no = intersect(stim_type,no_loco_ind);
    yes = intersect(stim_type,yes_loco_ind);
elseif ~isempty(cell_type) && ~isempty(stim_type)
    all = intersect(cell_type,intersect(stim_type,all_loco_ind)); % sort correlations by Cell and Stim types
    no = intersect(cell_type,intersect(stim_type,no_loco_ind));
    yes = intersect(cell_type, intersect(stim_type,yes_loco_ind));
end

%% Plot graphs

if contains(Ops.stim,'pontaneous'); loops = 2; % if spont, only trace data
else; loops = 6; end

fig_count = 1;
color = ['k' 'b' 'c'];

% pull out distances for all, no loco, and loco data
Distance{1} = Correlations_Matrix_All.Distance(all);
Distance{2} = Correlations_Matrix_All.Distance(no);
Distance{3} = Correlations_Matrix_All.Distance(yes);

titles = {'Trace Max Correlations', 'Trace Zero Correlations', ...
    'Noise Correlations', 'Signal Correlations',...
    'Noise XCorrelations', 'Signal XCorrelations'};
sig_titles = {'Significant Trace Max Correlations', 'Significant Trace Zero Correlations',...
    'Significant Noise Correlations', 'Significant Signal Correlations',...
    'Significant NoiseX Correlations', 'SignificantX Signal Correlations',};

correlation_type = {'TraceMax', 'TraceZero','Noise', 'Signal','NoiseX','SignalX'};
loco_type = {'all', 'no loco','loco'};

for r = 1:loops % loop through correlation types

Sig_Correlations{1} = eval(strcat('Correlations_Matrix_All.',correlation_type{r}, 'Correlations(intersect(all,', correlation_type{r}, '_Sig))'));
Sig_Correlations{2} = eval(strcat('Correlations_Matrix_All.',correlation_type{r}, 'Correlations(intersect(no,', correlation_type{r}, '_Sig))'));
Sig_Correlations{3} = eval(strcat('Correlations_Matrix_All.',correlation_type{r}, 'Correlations(intersect(yes,', correlation_type{r}, '_Sig))'));
Sig_Distance{1} = eval(strcat('Correlations_Matrix_All.Distance(intersect(all,', correlation_type{r}, '_Sig))'));
Sig_Distance{2} = eval(strcat('Correlations_Matrix_All.Distance(intersect(no,', correlation_type{r}, '_Sig))'));
Sig_Distance{3} = eval(strcat('Correlations_Matrix_All.Distance(intersect(yes,', correlation_type{r}, '_Sig))'));
Correlations{1} = eval(strcat('Correlations_Matrix_All.',correlation_type{r}, 'Correlations(all)'));
Correlations{2} = eval(strcat('Correlations_Matrix_All.',correlation_type{r}, 'Correlations(no)'));
Correlations{3} = eval(strcat('Correlations_Matrix_All.',correlation_type{r}, 'Correlations(yes)'));

    % Histograms of Loco, Non Loco and All - Correlations
    figures{1,fig_count} = figure;
    figure(figures{1,fig_count});
    for n=1:3
        if Ops.histogram_line == 1   % if opted for figure with outlined histograms
            histogram(Correlations{n},Ops.Ops.bins, 'Normalization',Ops.histogram_type,'DisplayStyle','stairs','LineWidth', 2, 'EdgeColor',color(n)); hold on
        else % standard histogram
            histogram(Correlations{n},Ops.bins, 'Normalization',Ops.histogram_type,'EdgeColor',color(n),'FaceColor', color(n)); hold on;
        end
        title(strcat(titles(r), ' ', Ops.celltype));
        xlabel('correlations'); ylabel(Ops.histogram_type);
    end
    fig_count = fig_count + 1;
    legend('All', 'No Loco','Loco');

    % Histograms of Loco, Non Loco and All - Only Significant Correlations
    figures{1,fig_count} = figure;
    figure(figures{1,fig_count});
    for n=1:3
        if Ops.histogram_line == 1   % if opted for figure with outlined histograms
            histogram(Sig_Correlations{n},Ops.bins, 'Normalization',Ops.histogram_type,'DisplayStyle','stairs','LineWidth', 2, 'EdgeColor',color(n)); hold on
        else % standard histogram
            histogram(Sig_Correlations{n},Ops.bins, 'Normalization',Ops.histogram_type,'EdgeColor',color(n),'FaceColor', color(n)); hold on;
        end
        title(strcat(sig_titles(r), ' ', Ops.celltype));
        xlabel('correlations'); ylabel(Ops.histogram_type);
        num_C(n) = length(Correlations{n}); % number of total cells
        num_sig(n) = length(Sig_Correlations{n}); % number of significantly correlated cells (with motor)
        ratio(n) = round((num_sig(n)/num_C(n))*100,2); % percentage of sig correlated cells

        indices_neg_correlations = find(Sig_Correlations{n}<0);
        indices_pos_correlations = find(Sig_Correlations{n}>0);

        Neg_Correlations = Sig_Correlations{n}(indices_neg_correlations);
        Pos_Correlations = Sig_Correlations{n}(indices_pos_correlations);
        ratio_neg(n) = round((length(Neg_Correlations)/num_C(n))*100,2); % percentage of sig correlated cells
        ratio_pos(n) = round((length(Pos_Correlations)/num_C(n))*100,2); % percentage of sig correlated cells

    end
    fig_count = fig_count + 1;

    subtitle({ 
        ['All loco = ' num2str(num_sig(1))  '/', num2str(num_C(1)) ' =', num2str(ratio(1)) '% sig, ' num2str(ratio_neg(1)) '% neg, ',...
        num2str(ratio_pos(1)) '% pos']
        ['No loco = ' num2str(num_sig(2))  '/', num2str(num_C(2)) ' =', num2str(ratio(2)) '% sig, ' num2str(ratio_neg(2)) '% neg, ',...
        num2str(ratio_pos(2)) '% pos']
        ['Loco = ' num2str(num_sig(3))  '/', num2str(num_C(3)) ' =', num2str(ratio(3)) '% sig, ' num2str(ratio_neg(3)) '% neg, ',...
        num2str(ratio_pos(3)) '% pos']
        });
    legend('All', 'No Loco','Loco');

    % Plot Correlations as a Function of Distance
    for n=1:3
        indices_neg_correlations = find(Sig_Correlations{n}<0);
        indices_pos_correlations = find(Sig_Correlations{n}>0);

        Neg_Correlations = Sig_Correlations{n}(indices_neg_correlations);
        Pos_Correlations = Sig_Correlations{n}(indices_pos_correlations);

        pairwise_distances{n} = Sig_Distance{n}; % distances of significantly correlated

        distance_neg = pairwise_distances{n}(indices_neg_correlations);
        distance_pos = pairwise_distances{n}(indices_pos_correlations);

        Neg_Fit = fitlm(distance_neg, Neg_Correlations, 'linear'); % fit linear regression model to negative correlations
        coefs_neg = Neg_Fit.Coefficients.Estimate;
        Pos_Fit = fitlm(distance_pos, Pos_Correlations, 'linear'); % fit linear regression model to positive correlations
        coefs_pos = Pos_Fit.Coefficients.Estimate;

        figures{1,fig_count} = figure;
        figure(figures{1,fig_count});
        scatter(Distance{n},Correlations{n},'k'); hold on;
        scatter(distance_pos, Pos_Correlations,'r'); hold on; % plot sig correlated in red
        scatter(distance_neg, Neg_Correlations,'c'); hold on; % plot sig correlated in blue
        hline1 = refline(coefs_neg(2),coefs_neg(1));
        hline2 =  refline(coefs_pos(2),coefs_pos(1));
        hline1.Color = 'c';
        hline2.Color = 'r';
        xlabel('distance (microns)'); ylabel('correlations');
        title([correlation_type{r} ' Correlations by Distance '  loco_type{n} 'For cell type ' Ops.celltype])
        legend('All', ['Sig Positive ' num2str(Pos_Fit.Rsquared.Adjusted)], ['Sig Negative ', num2str(Neg_Fit.Rsquared.Adjusted)])
        fig_count = fig_count + 1;
    end
end

%% Plot Stability of Correlations Across Recording Sessions

if Ops.plot_stability

    if contains(Ops.stim,'pontaneous'); loops = 2; % if spont, only trace data
    else; loops = 4; end

    for r=1:loops % trace, noise and signal correlations

        if r == 1; Correlations = Correlations_Matrix_All.TraceMaxCorrelations; % trace max
        elseif r == 2; Correlations = Correlations_Matrix_All.TraceZeroCorrelations; % trace zero
        elseif r == 3; Correlations = Correlations_Matrix_All.NoiseCorrelations; % noise
        elseif r == 4; Correlations = Correlations_Matrix_All.SignalCorrelations; % signal
        end

        titles = {'Trace Max', 'Trace Zero', 'Noise','Signal'};

        ind1 = 1;
        ind2 = 1;
        for i = 1:size(mice) % for each mouse
            % for ii = 1:size(blocks)
            %    indices = intersect(blocks{ii}, mouse_indices{i});
            matched_cells = unique(Correlations_Matrix_All.Cell1(mouse_indices{i})); % unique cells per mouse
            for ii = 1:size(matched_cells)-1
                cell1 = matched_cells(ii);
                cell1_indices = [find(Correlations_Matrix_All.Cell1==cell1); find(Correlations_Matrix_All.Cell2==cell1)];
                for iii = ii+1:size(matched_cells)
                    cell2 = matched_cells(iii);
                    cell2_indices = [find(Correlations_Matrix_All.Cell1==cell2); find(Correlations_Matrix_All.Cell2==cell2)];
                    AA = intersect(cell1_indices,cell2_indices,'rows');
                    BB = intersect(mouse_indices{i},AA,'rows');
                    CC = intersect(BB, all_loco_ind,'rows');
                    DD = intersect(CC, stim_type,'rows');
                    %EE = intersect(DD,Noise_Sig,'rows');

                    matched_FOV_correlations = rmmissing(Correlations(DD));
                    n_FOVs = size(matched_FOV_correlations,1); % number of FOVs with comparisons
                    if n_FOVs == Ops.timepoints
                        var_correlations(ind1) = var(rmmissing(Correlations(DD)));
                        mice_all{ind1} = mice(i);
                        cell1_all(ind1) = cell1;
                        cell2_all(ind1) = cell2;
                        ind1=ind1+1;
                    end

                end
            end

            HH = intersect(mouse_indices{i}, all_loco_ind, 'rows');
            II = intersect(HH, find(~isnan(Correlations)));

            % generate control distribution with random comparisons
            % between different cells in the same mouse
            for iiii = 1:Ops.nshuffles
                %GG = intersect(Noise_Sig,,'rows');
                gen_rand_cells = randsample(size(II,1),3);
                JJ = Correlations(II(gen_rand_cells)); % random cells within that mouse FOVs from any day
                var_correlation_control(ind2) = var(JJ);
                ind2=ind2+1;
            end
        end
        p_value = 0.05;
        [z_test, z_P, z_CI, z_R] = ztest(mean(var_correlations), mean(var_correlation_control), var(var_correlation_control), "Alpha", p_value);
        Ops.histogram_type = 'probability'
        Ops.bins = [-0.4:0.01:0.4];
        figure(figures{1,fig_count});
        histogram(var_correlations,Ops.bins, 'Normalization',Ops.histogram_type,'EdgeColor','r','FaceColor', 'r'); hold on;
        histogram(var_correlation_control,Ops.bins, 'Normalization',Ops.histogram_type,'EdgeColor','k','FaceColor', 'k'); hold on;
        vline(mean(var_correlations), 'r');
        vline(mean(var_correlation_control),'k');
        title([titles(r) 'Correlation var vs Shuffled']);
        subtitle(['ztest=' num2str(z_test) ' P =' num2str(z_P)]);
        xlabel('var of correlations'); ylabel(Ops.histogram_type);
        legend('Matched Cells', 'UnMatched Cell Control');
        fig_count = fig_count + 1;
    end

    % Table of Variances
    Variance_across_FOVs = table;
    Variance_across_FOVs = table(var_correlations, mice_all, cell1_all, cell2_all);
    cd(save_path)
    writetable(Variance_across_FOVs,'Variance_across_FOVs.csv','Delimiter',',');
    save('Variance_Table.mat', 'Variance_across_FOVs') ;

end

%% Save Figures
if Ops.save_figures
    if ~isempty(save_path)
        cd(save_path)
        if contains(Ops.stim,'pontaneous')
            figureNames = {'Hist_TraceMaxCorrelations'...
                'Hist_SigTraceMaxCorrelations'...
                'TraceMax_Corr_Distance_All'...
                'TraceMax_Corr_Distance_No_Loco '....
                'TraceMax_Corr_Distance_Loco '....
                'Hist_TraceZeroCorrelations'...
                'Hist_SigTraceZeroCorrelations'...
                'TraceZero_Corr_Distance_All'...
                'TraceZero_Corr_Distance_No_Loco '....
                'TraceZero_Corr_Distance_Loco '....
                'TraceMax_Corr_Stability'...
                'TraceZero_Corr_Stability'...
                };
        else
            figureNames = {'Hist_TraceMaxCorrelations'...
                'Hist_SigTraceMaxCorrelations'...
                'TraceMax_Corr_Distance_All'...
                'TraceMax_Corr_Distance_No_Loco'....
                'TraceMax_Corr_Distance_Loco'....
                'Hist_TraceZeroCorrelations'...
                'Hist_SigTraceZeroCorrelations'...
                'TraceZero_Corr_Distance_All'...
                'TraceZero_Corr_Distance_No_Loco '....
                'TraceZero_Corr_Distance_Loco '.... 
                'Hist_NoiseCorrelations'...
                'Hist_SigNoiseCorrelations' ...
                'Noise_Corr_Distance_All'...
                'Noise_Corr_Distance_No_Loco '....
                'Noise_Corr_Distance_Loco '....
                'Hist_SignalCorrelations' ...
                'Hist_SigSignalCorrelations' ...
                'Signal_Corr_Distance_ All'....
                'Signal_Corr_Distance_No_Loco '....
                'Signal_Corr_Distance_Loco '...
                'Hist_NoiseXCorrelations'...
                'Hist_SigNoiseXCorrelations' ...
                'NoiseX_Corr_Distance_All'...
                'NoiseX_Corr_Distance_No_Loco '....
                'NoiseX_Corr_Distance_Loco '....
                'Hist_SignalXCorrelations' ...
                'Hist_SigSignalXCorrelations' ...
                'SignalX_Corr_Distance_ All'....
                'SignalX_Corr_Distance_No_Loco '....
                'SignalX_Corr_Distance_Loco '...
                'TraceMax_Corr_Stability'...
                'TraceZero_Corr_Stability'...
                'Noise_Corr_Stability'...
                'Signal_Corr_Stability'
                };
        end

        for f = 1:size(figures,2)
            if isempty(figures{1,f}); continue; end
            h_name = ['Fig', num2str(f), '_', figureNames{f}];
            saveas(figures{1,f}, h_name, 'fig');
        end

    end
end