function [Summary, RowResultsTable, ColumnResultsTable] = plot_simple_ExtractedData_tuning(ExtractedData, UserOps)
% This code takes the row or column corresponding to the Best Response in the receptive field and averages for all cells

%% Options

normalizeToBestResponse = 1; %1 to normalize to BestResponse (BF), 0 to normalize to max value

%% Subset data from ExtractedData (check parameters at the top of function for options!)

%Return if no stim analysis
if isempty(fields(ExtractedData.StimAnalysis))
    return;
end

%I RECOMMEND NOT CHANGING THESE AND USING USEROPS INSTEAD
if nargin < 2
    UserOps = struct;
    UserOps.Loco                = 'All'; %All, Running, or NotRunning
    UserOps.ResponseType        = ''; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
    UserOps.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
    UserOps.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
    UserOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1)
    UserOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
    UserOps.sortbyCondition     = 0; %0 to ignore, 1 for first part, 2 for second part [e.g. Passive_60dB]
    UserOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
    UserOps.plotCorrelations	= 0; %0 = don't plot dependent variable correlations, 1 = plot
    UserOps.smooth_rasters      = 1; %0 = don't smooth, 1 = smooth
    UserOps.FibPhotChannel      = 2; %1 = blue, 2 = green, 3 = red, 4 = all of them
    UserOps.SuppressOutput      = 0; %0 to print Ops to command line, 1 to suppress
end

disp('Plotting tuning data')

StimData = simple_subset_ExtractedData(ExtractedData, UserOps);

Ops = StimData.Ops;
StimInfo = StimData.StimInfo;
Summary = StimData.Summary;
StimAnalysis = StimData.StimAnalysis;
CellDataByStim = StimData.CellDataByStim;
GroupList = Summary.Group;
Groups = unique(GroupList);

if isempty(Summary)
    disp('No responsive cells found to plot')
    return
end

%Make groups label without underscores for plots
GroupsLabel = Groups;
for g = 1:length(GroupsLabel)
    GroupsLabel{g} = regexprep(Groups{g},'_',' ','emptymatch');
end

%% Get read to plot average best row and best column

%Store "Best Row" and "Best Column" values

response_mat = StimAnalysis.Response;

%Switch row or column variable based on stim type
switch ExtractedData.Ops.stim_protocol
    case 2 %RF
        BestCol = StimAnalysis.BF;
        BestRow = StimAnalysis.BF_I;
        
    case 5 %SAM
        BestCol = StimAnalysis.BestRate;
        BestRow = StimAnalysis.BestDepth;
        
    case 6 %SAMfreq
        BestCol = StimAnalysis.BF;
        BestRow = StimAnalysis.BestDepth;
        
    otherwise
        error('This stim protocol is not coded for yet')
end

%Preallocate
[binaryAll_row_mat, binaryAct_row_mat, binaryPro_row_mat, binarySup_row_mat, binaryExc_row_mat] = deal(nan(height(StimAnalysis),length(StimInfo.V1)));
[binaryAll_col_mat, binaryAct_col_mat, binaryPro_col_mat, binarySup_col_mat, binaryExc_col_mat] = deal(nan(height(StimAnalysis),length(StimInfo.V2)));
[padded_binaryAll_row_mat, padded_binaryAct_row_mat, padded_binaryPro_row_mat, padded_binarySup_row_mat, padded_binaryExc_row_mat] = deal(nan(height(StimAnalysis),length(StimInfo.V1)*2-1));
[padded_binaryAll_col_mat, padded_binaryAct_col_mat, padded_binaryPro_col_mat, padded_binarySup_col_mat, padded_binaryExc_col_mat] = deal(nan(height(StimAnalysis),length(StimInfo.V2)*2-1));
[row_mat, row_mat_norm] = deal(nan(height(StimAnalysis),length(StimInfo.V1)));
[col_mat, col_mat_norm] = deal(nan(height(StimAnalysis),length(StimInfo.V2)));
[padded_row_mat, padded_row_norm] = deal(nan(height(StimAnalysis),length(StimInfo.V1)*2-1));
[padded_col_mat, padded_col_norm] = deal(nan(height(StimAnalysis),length(StimInfo.V2)*2-1));
padded_row_center = length(StimInfo.V1);
padded_col_center = length(StimInfo.V2);
Summary.BestRow_Slope = nan(height(StimAnalysis),1);
Summary.BestColumn_Slope = nan(height(StimAnalysis),1);

%Determine if Summary already contains best sparseness columns
compute_sparseness = 0;
if ~any(ismember(StimAnalysis.Properties.VariableNames, 'V1_BestSparseness'))
    compute_sparseness = 1;
    Summary.V1_BestSparseness = nan(height(StimAnalysis),1);
    Summary.V2_BestSparseness = nan(height(StimAnalysis),1);
    Summary.Binary_V1_BestSparseness = nan(height(StimAnalysis),1);
    Summary.Binary_V2_BestSparseness = nan(height(StimAnalysis),1);
end

%Loop through each cell to calculate normalized response in best row and column
for c = 1:height(StimAnalysis)
    
    %get index corresponding to BestRow and BestCol
    BestCol_ind = find(StimInfo.V1 == BestCol(c));
    BestRow_ind = find(StimInfo.V2 == BestRow(c));

    %compute and add best sparseness columns to Summary if not already there
    if compute_sparseness
        Summary.V1_BestSparseness(c) = StimAnalysis.V1_Sparseness{c}(BestRow_ind);
        Summary.V2_BestSparseness(c) = StimAnalysis.V2_Sparseness{c}(BestCol_ind);
        Summary.Binary_V1_BestSparseness(c) = StimAnalysis.Binary_V1_Sparseness{c}(BestRow_ind);
        Summary.Binary_V2_BestSparseness(c) = StimAnalysis.Binary_V2_Sparseness{c}(BestCol_ind);
    end
    
    %Get row corresponding to BestResponse (V2)
    row_mat(c,:) = response_mat(c,StimInfo.Combined(:,2) == BestRow(c));
    row_mat(c,isinf(row_mat(c,:))) = nan; %Fix infinity
    
    %Get column corresponding to BestResponse (V1)
    col_mat(c,:) = response_mat(c,StimInfo.Combined(:,1) == BestCol(c));
    col_mat(c,isinf(col_mat(c,:))) = nan; %Fix infinity
    
    %Get best response
    if ~isequal(col_mat(c,BestRow_ind), row_mat(c,BestCol_ind))
        error('Troubleshoot') %These should be the same
    else
        BestResponse = row_mat(c,BestCol_ind);
    end

    %STORE BINARY VALUES
    IsRF = CellDataByStim(c).RF';
    ResponseType = CellDataByStim(c).PeakData.ResponseType';
    
    %Four types of binary mats
    binaryAll = IsRF;
    binaryAct = strcmp(ResponseType, 'activated');
    binaryPro = strcmp(ResponseType, 'prolonged');
    binarySup = strcmp(ResponseType, 'suppressed');
    binaryExc = binaryAct + binaryPro;
    
    binaryAll_row_mat(c,:) = binaryAll(StimInfo.Combined(:,2) == BestRow(c)); 
    binaryAct_row_mat(c,:) = binaryAct(StimInfo.Combined(:,2) == BestRow(c)); 
    binaryPro_row_mat(c,:) = binaryPro(StimInfo.Combined(:,2) == BestRow(c)); 
    binarySup_row_mat(c,:) = binarySup(StimInfo.Combined(:,2) == BestRow(c)); 
    binaryExc_row_mat(c,:) = binaryExc(StimInfo.Combined(:,2) == BestRow(c)); 
    binaryAll_col_mat(c,:) = binaryAll(StimInfo.Combined(:,1) == BestCol(c)); 
    binaryAct_col_mat(c,:) = binaryAct(StimInfo.Combined(:,1) == BestCol(c)); 
    binaryPro_col_mat(c,:) = binaryPro(StimInfo.Combined(:,1) == BestCol(c)); 
    binarySup_col_mat(c,:) = binarySup(StimInfo.Combined(:,1) == BestCol(c)); 
    binaryExc_col_mat(c,:) = binaryExc(StimInfo.Combined(:,1) == BestCol(c));

    %NORMALIZE
    temp_row = row_mat(c,:);
    temp_column = col_mat(c,:);
    
    %If inhibitory RF, flip sign before normalizing (use raw mat to see negative results)
    if strcmp(Summary.RF_Type(c), 'inhibitory')
        temp_row = -temp_row;
        temp_column = -temp_column;
        BestResponse = -BestResponse;
    elseif BestResponse < 0 || max(temp_row) < 0 || max(temp_column) < 0
        error('Troubleshoot') %max response shouldn't be negative unless RF is inhibitory
    end

    if normalizeToBestResponse 
        %Set BestResponse as max positive and negative value
        temp_row(temp_row > BestResponse) = BestResponse;
        temp_column(temp_column > BestResponse) = BestResponse;
        temp_row(temp_row < -BestResponse) = -BestResponse;
        temp_column(temp_column < -BestResponse) = -BestResponse;
    else
        %Set max positive response as max positive and negative value
        temp_row(temp_row < -max(temp_row)) = -max(temp_row);
        temp_column(temp_column < -max(temp_column)) = -max(temp_column);
    end
    
    %Perform final normalization and store result
    row_mat_norm(c,:) = temp_row./max(temp_row);
    col_mat_norm(c,:) = temp_column./max(temp_column);
  
    %Shift data to store in padded matrices
    row_ind_to_shift = padded_row_center - BestCol_ind;
    padded_row_mat(c,row_ind_to_shift+(1:length(StimInfo.V1))) = row_mat(c,:);
    padded_row_norm(c,row_ind_to_shift+(1:length(StimInfo.V1))) = row_mat_norm(c,:);
    padded_binaryAll_row_mat(c,row_ind_to_shift+(1:length(StimInfo.V1))) = binaryAll_row_mat(c,:);
    padded_binaryAct_row_mat(c,row_ind_to_shift+(1:length(StimInfo.V1))) = binaryAct_row_mat(c,:);
    padded_binaryPro_row_mat(c,row_ind_to_shift+(1:length(StimInfo.V1))) = binaryPro_row_mat(c,:);
    padded_binarySup_row_mat(c,row_ind_to_shift+(1:length(StimInfo.V1))) = binarySup_row_mat(c,:);
    padded_binaryExc_row_mat(c,row_ind_to_shift+(1:length(StimInfo.V1))) = binaryExc_row_mat(c,:);
    
    col_ind_to_shift = padded_col_center - BestRow_ind;
    padded_col_mat(c,col_ind_to_shift+(1:length(StimInfo.V2))) = col_mat(c,:);
    padded_col_norm(c,col_ind_to_shift+(1:length(StimInfo.V2))) = col_mat_norm(c,:);
    padded_binaryAll_col_mat(c,col_ind_to_shift+(1:length(StimInfo.V2))) = binaryAll_col_mat(c,:);
    padded_binaryAct_col_mat(c,col_ind_to_shift+(1:length(StimInfo.V2))) = binaryAct_col_mat(c,:);
    padded_binaryPro_col_mat(c,col_ind_to_shift+(1:length(StimInfo.V2))) = binaryPro_col_mat(c,:);
    padded_binarySup_col_mat(c,col_ind_to_shift+(1:length(StimInfo.V2))) = binarySup_col_mat(c,:);
    padded_binaryExc_col_mat(c,col_ind_to_shift+(1:length(StimInfo.V2))) = binaryExc_col_mat(c,:);
    
    %Compute slope on normalized data
    row_coefficients = polyfit(1:length(row_mat_norm(c,:)), row_mat_norm(c,:), 1); % Get coefficients of a line fit through the data.
    Summary.BestRow_Slope(c) = row_coefficients(1);
    
    col_coefficients = polyfit(1:length(col_mat_norm(c,:)), col_mat_norm(c,:), 1); % Get coefficients of a line fit through the data.
    Summary.BestColumn_Slope(c) = col_coefficients(1);
end

%% Filter activity mats by binary mats?

filter_activity_mats_by_binary_mats = 0;

if filter_activity_mats_by_binary_mats
    padded_row_norm(~(padded_binaryAll_row_mat == 1)) = 0;  
end

%% Plot figure
   
figure; %('units','normalized','outerposition',[0 0 1 1]); hold on

for i = 1:4
    if i == 1
        data = padded_row_norm;
        plotTitle = 'Best Row Padded: Normalized';
        plotXTickLabel = (-length(StimInfo.V1)+1):1:(length(StimInfo.V1)-1);
        unit_ind = 1;
    elseif i == 2
        data = col_mat_norm;
        plotTitle = 'Best Column: Normalized';
        plotXTickLabel = StimInfo.V2;
        unit_ind = 2;
    elseif i == 3
        data = padded_row_mat;
        plotTitle = 'Best Row Padded: Raw';
        plotXTickLabel = (-length(StimInfo.V1)+1):1:(length(StimInfo.V1)-1);
        unit_ind = 1;
    elseif i == 4
        data = col_mat;
        plotTitle = 'Best Column: Raw';
        plotXTickLabel = StimInfo.V2;
        unit_ind = 2;
    end
    
    subplot(2,2,i); hold on
    for g = 1:length(Groups)
        isGroup = strcmp(Groups{g},GroupList);
        plotmat = data(isGroup,:);
        y = mean(plotmat,1, 'omitnan');
        nValsForSEM = sum(~isnan(plotmat),1) - 1;
        nValsForSEM(nValsForSEM < 0) = 0;
        y_sem = std(plotmat,[],1,'omitnan')./sqrt(nValsForSEM);
        errorbar(y, y_sem)
        set(gca, 'XTick', 1:length(y))
        set(gca, 'XTickLabel', plotXTickLabel)
        xlabel(StimInfo.Units{unit_ind})
    end
    title(plotTitle)
    legend(GroupsLabel)
end

tak_suptitle(StimInfo.StimType)
drawnow

%% Plot binary figure
   
figure; %('units','normalized','outerposition',[0 0 1 1]); hold on

for i = 1:10
    if i == 1
        data = padded_binaryAll_row_mat;
        plotTitle = 'Padded Best Row: Binary All';
        plotXTickLabel = (-length(StimInfo.V1)+1):1:(length(StimInfo.V1)-1);
        unit_ind = 1;
    elseif i == 2
        data = binaryAll_col_mat;
        plotTitle = 'Best Column: Binary All';
        plotXTickLabel = StimInfo.V2;
        unit_ind = 2;
    elseif i == 3
        data = padded_binaryAct_row_mat;
        plotTitle = 'Padded Best Row: Binary Act';
        plotXTickLabel = (-length(StimInfo.V1)+1):1:(length(StimInfo.V1)-1);
        unit_ind = 1;
    elseif i == 4
        data = binaryAct_col_mat;
        plotTitle = 'Best Column: Binary Act';
        plotXTickLabel = StimInfo.V2;
        unit_ind = 2;
    elseif i == 5
        data = padded_binaryPro_row_mat;
        plotTitle = 'Padded Best Row: Binary Pro';
        plotXTickLabel = (-length(StimInfo.V1)+1):1:(length(StimInfo.V1)-1);
        unit_ind = 1;
    elseif i == 6
        data = binaryPro_col_mat;
        plotTitle = 'Best Column: Binary Pro';
        plotXTickLabel = StimInfo.V2;
        unit_ind = 2;
    elseif i == 7
        data = padded_binarySup_row_mat;
        plotTitle = 'Padded Best Row: Binary Sup';
        plotXTickLabel = (-length(StimInfo.V1)+1):1:(length(StimInfo.V1)-1);
        unit_ind = 1;
    elseif i == 8
        data = binarySup_col_mat;
        plotTitle = 'Best Column: Binary Sup';
        plotXTickLabel = StimInfo.V2;
        unit_ind = 2;
    elseif i == 9
        data = padded_binaryExc_row_mat;
        plotTitle = 'Padded Best Row: Binary Exc';
        plotXTickLabel = (-length(StimInfo.V1)+1):1:(length(StimInfo.V1)-1);
        unit_ind = 1;
    elseif i == 10
        data = binaryExc_col_mat;
        plotTitle = 'Best Column: Binary Exc';
        plotXTickLabel = StimInfo.V2;
        unit_ind = 2;
    end
    
    subplot(5,2,i); hold on
    for g = 1:length(Groups)
        isGroup = strcmp(Groups{g},GroupList);
        plotmat = data(isGroup,:);
        y = mean(plotmat,1, 'omitnan');
        nValsForSEM = sum(~isnan(plotmat),1) - 1;
        nValsForSEM(nValsForSEM < 0) = 0;
        y_sem = std(plotmat,[],1,'omitnan')./sqrt(nValsForSEM);
        errorbar(y, y_sem)
        set(gca, 'XTick', 1:length(y))
        set(gca, 'XTickLabel', plotXTickLabel)
        xlabel(StimInfo.Units{unit_ind})
    end
    title(plotTitle)
    legend(GroupsLabel)
end

tak_suptitle(StimInfo.StimType)
drawnow

%% Plot Best Row by columns

figure; hold on
tiledlayout(7, length(StimInfo.V1))

for f = 1:7
    if f == 1
        data = row_mat_norm;
        ylab = 'Norm. Activity';
    elseif f == 2
        data = row_mat;
        ylab = 'Raw Activity';
    elseif f == 3
        data = binaryAll_row_mat;
        ylab = 'Binary All';
    elseif f == 4
        data = binaryAct_row_mat;
        ylab = 'Binary Act';
    elseif f == 5
        data = binaryPro_row_mat;
        ylab = 'Binary Pro';
    elseif f == 6
        data = binarySup_row_mat;
        ylab = 'Binary Sup';
    elseif f == 7
        data = binaryExc_row_mat;
        ylab = 'Binary Exc';
    end
    
    for v = 1:length(StimInfo.V1)
        nexttile; hold on
        currentV1 = StimInfo.V1(v);

        for g = 1:length(Groups)
            isGroup = strcmp(Groups{g},GroupList);
            temp_BestCol = BestCol(isGroup);
            temp_ind = temp_BestCol == currentV1;

            plotmat = data(isGroup,:);
            plotmat(~temp_ind,:) = [];
            y = mean(plotmat,1, 'omitnan');
            nValsForSEM = sum(~isnan(plotmat),1) - 1;
            nValsForSEM(nValsForSEM < 0) = 0;
            y_sem = std(plotmat,[],1,'omitnan')./sqrt(nValsForSEM);
            errorbar(y, y_sem)
            set(gca, 'XTick', 1:length(y))
            set(gca, 'XTickLabel', StimInfo.V2)
            xlabel(StimInfo.Units{2})
            if v == 1
                ylabel(ylab)
            end
        end
        title([num2str(currentV1) StimInfo.Units{1}])
        if f == 1
            legend(GroupsLabel)
        end
    end
end

tak_suptitle(StimInfo.StimType)
drawnow

%% Plot Best Column by rows
% 
% figure; hold on
% tiledlayout(7, length(StimInfo.V2))
% 
% for f = 1:7
%     if f == 1
%         data = col_mat_norm;
%         ylab = 'Norm. Activity';
%     elseif f == 2
%         data = col_mat;
%         ylab = 'Raw Activity';
%     elseif f == 3
%         data = binaryAll_col_mat;
%         ylab = 'Binary All';
%     elseif f == 4
%         data = binaryAct_col_mat;
%         ylab = 'Binary Act';
%     elseif f == 5
%         data = binaryPro_col_mat;
%         ylab = 'Binary Pro';
%     elseif f == 6
%         data = binarySup_col_mat;
%         ylab = 'Binary Sup';
%     elseif f == 7
%         data = binaryExc_col_mat;
%         ylab = 'Binary Exc';
%     end
%     
%     for v = 1:length(StimInfo.V2)
%         nexttile; hold on
%         currentV2 = StimInfo.V2(v);
% 
%         for g = 1:length(Groups)
%             isGroup = strcmp(Groups{g},GroupList);
%             temp_BestRow = BestRow(isGroup);
%             temp_ind = temp_BestRow == currentV2;
% 
%             plotmat = data(isGroup,:);
%             plotmat(~temp_ind,:) = [];
%             y = mean(plotmat,1, 'omitnan');
%             nValsForSEM = sum(~isnan(plotmat),1) - 1;
%             nValsForSEM(nValsForSEM < 0) = 0;
%             y_sem = std(plotmat,[],1,'omitnan')./sqrt(nValsForSEM);
%             errorbar(y, y_sem)
%             set(gca, 'XTick', 1:length(y))
%             set(gca, 'XTickLabel', StimInfo.V1)
%             xlabel(StimInfo.Units{1})
%             if v == 1
%                 ylabel(ylab)
%             end
%         end
%         title([num2str(currentV2) StimInfo.Units{2}])
%         if f == 1
%             legend(GroupsLabel)
%         end
%     end
% end
% 
% tak_suptitle(StimInfo.StimType)
% drawnow

%% Make Results Table

for i = 1:height(Summary)
    
    TempResultsTable = repmat(Summary(i,:), length(StimInfo.V1), 1);
    TempResultsTable.Variable = StimInfo.V1;
    TempResultsTable.Raw = row_mat(i,:)';
    TempResultsTable.Norm = row_mat_norm(i,:)';
    TempResultsTable.BinaryAll = binaryAll_row_mat(i,:)';
    TempResultsTable.BinaryAct = binaryAct_row_mat(i,:)';
    TempResultsTable.BinaryPro = binaryPro_row_mat(i,:)';
    TempResultsTable.BinarySup = binarySup_row_mat(i,:)';
    TempResultsTable.BinaryExc = binaryExc_row_mat(i,:)';
        
    if i == 1
        RowResultsTable = TempResultsTable;
    else
        RowResultsTable = [RowResultsTable; TempResultsTable];
    end
end

for i = 1:height(Summary)
    
    TempResultsTable = repmat(Summary(i,:), length(StimInfo.V2), 1);
    TempResultsTable.Variable = StimInfo.V2;
    TempResultsTable.Raw = col_mat(i,:)';
    TempResultsTable.Norm = col_mat_norm(i,:)';
    TempResultsTable.BinaryAll = binaryAll_col_mat(i,:)';
    TempResultsTable.BinaryAct = binaryAct_col_mat(i,:)';
    TempResultsTable.BinaryPro = binaryPro_col_mat(i,:)';
    TempResultsTable.BinarySup = binarySup_col_mat(i,:)';
    TempResultsTable.BinaryExc = binaryExc_col_mat(i,:)';
    
    if i == 1
        ColumnResultsTable = TempResultsTable;
    else
        ColumnResultsTable = [ColumnResultsTable; TempResultsTable];
    end
end
