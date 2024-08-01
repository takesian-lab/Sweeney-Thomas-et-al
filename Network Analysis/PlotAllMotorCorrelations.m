function PlotAllMotorCorrelations(Correlations_Matrix_All,save_path, Ops)
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
    figures = cell(1,24);
end

%% Sort Data According to Loco, Cell Type and Stim Type

% Sort by Cell Types
cell_type = find(strcmp(Correlations_Matrix_All.MotorCellType,Ops.celltype));

% Sort by Stim Type
stim_type = find(strcmp(Correlations_Matrix_All.MotorStimName,Ops.stim));

% Sort by Mouse_Name
Correlations_Matrix_All.MotorMouseName = categorical(Correlations_Matrix_All.MotorMouseName);
mice = unique(Correlations_Matrix_All.MotorMouseName);
for n = 1:size(mice,1)
    mouse_indices{n} = find(Correlations_Matrix_All.MotorMouseName == mice(n));
end

% Sort by Block_Name
Correlations_Matrix_All.MotorBlockName = categorical(Correlations_Matrix_All.MotorBlockName);
blocks = unique(Correlations_Matrix_All.MotorBlockName);
for n = 1:size(blocks,1)
    block_indices{n} = find(Correlations_Matrix_All.MotorBlockName == blocks(n));
end

% Significant Correlations
Trace_Sig.Motor = find(Correlations_Matrix_All.MotorZTest == 1);

% Determine if xyz correlations exist for controls
existXYCorr = any(strcmp('XYCorrelations',Correlations_Matrix_All.Properties.VariableNames));
existZCorr = any(strcmp('ZCorrelations',Correlations_Matrix_All.Properties.VariableNames));

if existXYCorr  % if XYCorrelations calculated
    if any(strcmp('XYZTest',Correlations_Matrix_All.Properties.VariableNames))
        Trace_Sig.XY = find(Correlations_Matrix_All.XYZTest == 1);
    elseif any(strcmp('XY_ZTest',Correlations_Matrix_All.Properties.VariableNames)) %OLD VARIABLE NAME (discontinue one day)
        Trace_Sig.XY = find(Correlations_Matrix_All.XY_ZTest == 1);
    end
end

if  existZCorr % if ZCorrelations calculated
    if any(strcmp('ZZTest',Correlations_Matrix_All.Properties.VariableNames))
        Trace_Sig.Z = find(Correlations_Matrix_All.ZZTest == 1);
    elseif any(strcmp('Z_ZTest',Correlations_Matrix_All.Properties.VariableNames)) %OLD VARIABLE NAME (discontinue one day)
        Trace_Sig.Z = find(Correlations_Matrix_All.Z_ZTest == 1);
    end
end


%% Plot graphs

if existXYCorr && existZCorr; loops = 3; % if XYZ controls exist
elseif existXYCorr && ~existZCorr; loops = 2;
else loops = 1;
end

fig_count = 1;
color = ['k' 'b' 'c'];

titles = {'Motor Correlation', 'XY Correlation', 'Z Correlation'};
sig_titles = {'Significant Motor Correlations', 'Significant XY Correlations', 'Significant Z Correlations'};
correlation_type = {'MotorCorrelations','XYCorrelations','ZCorrelations'};
sig_type = {'Motor', 'XY', 'Z'};

for r = 1:loops
    All_Correlations = eval(strcat('Correlations_Matrix_All.',correlation_type{r}));
    Sig_ind = eval(strcat('Trace_Sig.',sig_type{r}));

    if isempty(cell_type)
        Sig_Correlations = All_Correlations(Sig_ind);
        Correlations  = All_Correlations;
    else
        Sig_Correlations = All_Correlations(intersect(cell_type,Sig_ind));
        Correlations  = All_Correlations(cell_type);
    end

    % Histograms of Correlations
    figures{1,fig_count} = figure;
    figure(figures{1,fig_count});
    if Ops.histogram_line == 1   % if opted for figure with outlined histograms
        histogram(Correlations,Ops.Ops.bins, 'Normalization',Ops.histogram_type,'DisplayStyle','stairs','LineWidth', 2, 'EdgeColor',color(1)); hold on
    else % standard histogram
        histogram(Correlations,Ops.bins, 'Normalization',Ops.histogram_type,'EdgeColor',color(1),'FaceColor', color(1)); hold on;
    end
    title(strcat(titles{r}, ' ', '  ', Ops.celltype, ' cells'));
    num_C = length(Correlations); % number of total cells
    num_sig = length(Sig_Correlations); % number of significantly correlated cells (with motor)
    ratio = round((num_sig/num_C)*100,2); % percentage of sig correlated cells

    indices_neg_correlations = find(Sig_Correlations<0);
    indices_pos_correlations = find(Sig_Correlations>0);

    Neg_Correlations = Sig_Correlations(indices_neg_correlations);
    Pos_Correlations = Sig_Correlations(indices_pos_correlations);
    ratio_neg = round((length(Neg_Correlations)/num_C)*100,2); % percentage of sig correlated cells
    ratio_pos = round((length(Pos_Correlations)/num_C)*100,2); % percentage of sig correlated cells

    subtitle(strcat(num2str(num_sig), '/', num2str(num_C), ' =', num2str(ratio), '% sig, ', num2str(ratio_neg),'% neg, ',...
        num2str(ratio_pos),'% pos'));
    xlabel('correlations'); ylabel(Ops.histogram_type);
hold on;

% Histograms of Loco, Non Loco and All - Only Significant Correlations
if Ops.histogram_line == 1   % if opted for figure with outlined histograms
    histogram(Sig_Correlations,Ops.bins, 'Normalization',Ops.histogram_type,'DisplayStyle','stairs','LineWidth', 2, 'EdgeColor','c'); hold on
else % standard histogram
    histogram(Sig_Correlations,Ops.bins, 'Normalization',Ops.histogram_type,'EdgeColor',color(1),'FaceColor', 'c', 'FaceAlpha',1); hold on;
end
xlabel('correlations'); ylabel(Ops.histogram_type);

fig_count = fig_count + 1;
end


%% Report % of Sig Correlations for both Locomotor and XY 

% plot Motor Correlations versus XY Correlations
figures{1,fig_count} = figure;
figure(figures{1,fig_count});
scatter(Correlations_Matrix_All.XYCorrelations, Correlations_Matrix_All.MotorCorrelations); hold on;
title('XY versus Motor Correlations')
xlabel('XY Correlations'); ylabel('Motor Correlations');

% fit line and report rsquare
mdl = fitlm(Correlations_Matrix_All.XYCorrelations, Correlations_Matrix_All.MotorCorrelations);
coefs = mdl.Coefficients.Estimate;
refline(coefs(2),coefs(1));

% report % of correlations both significant for motor and XY
both_sig = length(intersect(Trace_Sig.XY, Trace_Sig.Motor)); % both sig
motor_sig = length(Trace_Sig.Motor); % total motor correlations
percent_share = (both_sig./motor_sig)*100;
subtitle(strcat('Percentage both motor and XY correlations = ', num2str(percent_share),' ', 'Rsquare adjusted = ', num2str(mdl.Rsquared.adjusted)));
fig_count = fig_count + 1;

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
            figureNames = {'Hist_Motor_Correlations'...
                'Hist_XY_Correlations'...
                'Hist_Z_Correlations'...
                'Motor versus XY Correlations'
                };

        for f = 1:size(figures,2)
            if isempty(figures{1,f}); continue; end
            h_name = ['Fig', num2str(f), '_', figureNames{f}];
            saveas(figures{1,f}, h_name, 'fig');
        end

    end
end