%% Make an average padded RF for each group

%Version History:
%V1 = archived 8/29/2023 MET
%V2 = this version, adds more normalization options and figures and reorganizes code
%% User options

%Load ExtractedData first
clearvars -except ExtractedData
 
%OPTIONS (applied to all plots)
save_stats_table = 1; %excel spreadsheet with data for stats
save_figures = 1; %using save_plots
center_RF_dimension = 1; %1 for horizontal, 2 for vertical
control_group = "PYR"; %"" for pvalue plots
sort_by_GCAMP = 0; %0 for groups, 1 for gcamp, 2 to combine
normalize_to_raw_BF = 0; %normalize and center to BF computed from raw data
use_CF_instead_of_BF = 0; %applies to stim protocol 2 only
use_mask = 1; %apply binary mask to raw data
use_reliability = 0; %Plot reliability R value in heat map instead of AUC
recompute_bandwidth = 0; %recompute gaussian bandwidth at BF for each cell (WARNING: this adds computation time)
nanpad_matrices = 0; %0 for zero padding, 1 for nan-padding
save_path = pwd;

%FILTERS (Script will plot figures for each combination of these elements)
%Do not leave {} empty -> must be {''}
LOCO = {'All'}; %'All', 'Running', 'NotRunning'
PLOTBINARY = [0 1]; %[0 1]
RESPONSETYPE = {'activated'}; %{activated, 'prolonged', 'suppressed'}; {'', 'excitatory'}; OR {''}
RFTYPE = {''}; %, 'excitatory', 'inhibitory'}; OR {''}
RFMASK = {''}; %'excitatory', 'suppressed' plot only RF responses that are excitatory(activated + prolonged) or suppressed

cd(save_path)

%% Stim-specific setup

%Specify best frequency columns per stimulus
switch ExtractedData.Ops.stim_protocol
    case 2 %RF
        
        if use_CF_instead_of_BF
            v1_column = 'CF'; 
            v2_column = 'CF_I';
            if normalize_to_raw_BF
                %Do not try to compute raw CF without thinking more about it...
                %E.g. how do we determine tip of V-shaped tuning curve without R&R?
                warning('Normalize to raw BF not applicable to CF. Using R&R CF')
            end
        else
            v1_column = 'BF'; 
            v2_column = 'BF_I';
        end
        
    case 3 %FM
        v1_column = 'BestSpeed';
        v2_column = 'BestSpeed';
         
    case 5 %SAM
        v1_column = 'BestRate'; 
        v2_column = 'BestDepth';
        
    case 6 %SAM freq
        v1_column = 'BF'; 
        v2_column = 'BestDepth';
        
    otherwise
        error('Add info for this stim procotol')
end

%Setup for RF Mapping
Parameters = ExtractedData.StimInfo.Parameters;
V1 = ExtractedData.StimInfo.V1;
V2 = ExtractedData.StimInfo.V2;
Order = [1 2];
Combined = ExtractedData.StimInfo.Combined; %Original combined values
if center_RF_dimension == 2
    [v1_column, v2_column] = deal(v2_column, v1_column); %Swap v1_column and v2_column
    [V1, V2] = deal(V2, V1); %Swap V1 and V2
    Parameters = fliplr(Parameters);
    Order = [2 1];
end

%Remove 0% row from SAM
if ExtractedData.Ops.stim_protocol == 5 
    V1(V1 == 0) = [];
    V2(V2 == 0) = [];
end

%Generate map indices (see add_RF_mapping_to_StimInfo)
RF_Map = nan(length(V1)*length(V2),2);
OriginalOrder = nan(length(V1)*length(V2),1);
count = 1;
for vv = 1:length(V2)
    for v = 1:length(V1)
        RF_Map(count,1) = vv; %Y
        RF_Map(count,2) = v; %X
        OriginalOrder(count,1) = find(((Combined(:,Order(1)) == V1(v)) + (Combined(:,Order(2)) == V2(vv))) == 2);
        count = count + 1;
    end
end
if nanpad_matrices
    nanmat = nan(length(V2), length(V1));
else
    nanmat = zeros(length(V2), length(V1));
end
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

%Pad and center matrix
dims = size(nanmat);
if nanpad_matrices
    padded_nanmat = nan(dims(1),dims(2)+dims(2)-1); % pad matrix
else
    padded_nanmat = zeros(dims(1),dims(2)+dims(2)-1); % pad matrix
end
centerColumn = ceil(size(padded_nanmat,2)./2);
padded_sparseness_mat = nan(1,dims(2)+dims(2)-1); %Does not need to switch between 0/nan padding

%Create new StimInfo to use in simple_compute_tuning_sparseness script
StimInfo = struct;
StimInfo.Parameters = Parameters;
StimInfo.V1 = V1;
StimInfo.V2 = V2;
StimInfo.RF_Map = RF_Map;

%% Loop through all filtering options

Filters = struct; %Store filtered data and filter criteria here

count = 1;
for L = 1:length(LOCO)
    for B = 1:length(PLOTBINARY)
        for M = 1:length(RFMASK)
            for R = 1:length(RESPONSETYPE)
                for F = 1:length(RFTYPE)
                    %Record filters
                    Filters(count).DATA = {};
                    Filters(count).LOCO = LOCO{L};
                    Filters(count).PLOTBINARY = PLOTBINARY(B);
                    Filters(count).RFMASK = RFMASK(M);
                    Filters(count).RESPONSETYPE = RESPONSETYPE(R);
                    Filters(count).RFTYPE = RFTYPE(F);
                    Filters(count).TITLE = strcat({'Loco='}, LOCO{L}, {' PlotBinary='}, num2str(PLOTBINARY(B)), {' ResponseMask='}, string(RFMASK{M}),...
                    {' ResponseType='}, string(RESPONSETYPE{R}), {' RFType='}, string(RFTYPE{F}));

                    %Subset and store data
                    SubsetOps = struct;
                    SubsetOps.Loco                = LOCO{L};
                    SubsetOps.ResponseType        = RESPONSETYPE{R}; 
                    SubsetOps.RF_Type             = RFTYPE{F};
                    SubsetOps.IsResponsive        = 1; %Only look at responsive cells
                    SubsetOps.sortbyGCAMP         = sort_by_GCAMP;
                    SubsetOps.SuppressOutput      = 1;
                    Filters(count).DATA = simple_subset_ExtractedData(ExtractedData, SubsetOps); %This function filters the data based on the above ops

                    count = count + 1;
                end
            end
        end
    end
end

%% Loop through filters and center & normalize RFs

%Take group list from full extracted data in case not all groups are represented in each filter
GroupList = unique(ExtractedData.Summary.(LOCO{1}).Group);

for F = 1:size(Filters,2)
    
    %Get data for current filter set
    SubsetData = Filters(F).DATA;
    Summary = SubsetData.Summary;
    CellDataByStim = SubsetData.CellDataByStim;
    StimAnalysis = SubsetData.StimAnalysis;
    Groups = Summary.Group;

    %Initialize
    [Padded_RF, Padded_RF_BFRow, Padded_Sparseness, Bandwidth, Bandwidth_BF] = deal([]);

    %Go through all cells in Summary and store in padded mat
    for c = 1:height(Summary)
        
        %Initialize
        temp_padded_RF = padded_nanmat;
        temp_padded_sparseness = padded_sparseness_mat;
        
        %Get RF data (OriginalOrder reorders if needed for padding in vertical dimension)
        PeakData = StimAnalysis.Response(c,OriginalOrder);
        RFData = CellDataByStim(c).RF(OriginalOrder);
        RTData = CellDataByStim(c).PeakData.ResponseType(OriginalOrder);
        ReliabilityData = [CellDataByStim(c).ReliabilityData.R];
        ReliabilityData = ReliabilityData(OriginalOrder);
        %ReliabilityData(~RFData) = 0;
        
        %Put into RF dimensions
        [vmat, IsRF, vmatR] = deal(nanmat);
        vmat(RF_ind) = PeakData;
        IsRF(RF_ind) = RFData;
        vmatR(RF_ind) = ReliabilityData;
        IsRF(isnan(vmat)) = nan; %Replace IsRF with NaN wherever there weren't enough trials to compute a response (important for Running/NotRunning)
        
        %Correct for inf in data
        vmat(isinf(vmat)) = nan;
        
        %Assign Response Mask
        %--------------------
        if ~isempty(Filters(F).RFMASK{:})
            responseMASK = zeros(size(nanmat));
            tempmask = responseMASK(:);
            switch Filters(F).RFMASK{:}
                case 'activated'
                    tempmask(strcmp(RTData,'activated')) = 1;
                case 'prolonged'
                    tempmask(strcmp(RTData,'prolonged')) = 1;
                case 'excitatory'
                    tempmask(strcmp(RTData,'activated')) = 1;
                    tempmask(strcmp(RTData,'prolonged')) = 1;
                case 'suppressed'
                    tempmask(strcmp(RTData,'suppressed')) = 1;
                otherwise
                    error('Check for typos in RFMASK')
            end
            tempmask(~RFData) = 0; %Make sure we aren't including data for non-responsive stim
            responseMASK(RF_ind) = tempmask;
            IsRF = responseMASK;
        end
        %--------------------------
        
        %Find BF coord
        if normalize_to_raw_BF
            [~, maxind] = max(vmat(:),[],1);
            [BF_dim2, BF_dim1] = ind2sub(size(vmat), maxind);
        else
            BF_dim1 = find(V1 == StimAnalysis.(v1_column)(c));
            BF_dim2 = find(V2 == StimAnalysis.(v2_column)(c));
            
            %Stim specific computations
            if ExtractedData.Ops.stim_protocol == 5 && (isempty(BF_dim1) || isempty(BF_dim2))%For SAM when BF == 0%
                BF_dim1 = 1;
                BF_dim2 = 1;
                vmat(:) = NaN; %Replace data with NaNs for these cells
            elseif ExtractedData.Ops.stim_protocol == 3 %FM
                if center_RF_dimension == 1
                    BF_dim1 = find(V1 == abs(StimAnalysis.(v1_column)(c)));
                    BF_dim2 = find(V2 == sign(StimAnalysis.(v2_column)(c)));
                else
                    BF_dim1 = find(V1 == sign(StimAnalysis.(v1_column)(c)));
                    BF_dim2 = find(V2 == abs(StimAnalysis.(v2_column)(c)));
                end
            end
        end
        
        %NORMALIZE
        if Filters(F).PLOTBINARY
            vmat = IsRF;
        else
            BF_val = abs(vmat(BF_dim2, BF_dim1));
            vmat = vmat./BF_val;
            if use_mask
                vmat(IsRF == 0) = 0;
            end
        end
        
        %TEMP
        if use_reliability
            vmat = vmatR;
        end
        
        %PAD
        shift1 = centerColumn - BF_dim1 + 1;
        shift2 = size(vmat,2)+shift1-1;
        %temp_padded_RF(BF_dim2,shift1:shift2) = vmat(BF_dim2,:); %TEMP: test what would happen if we only filled BF row
        temp_padded_RF(:,shift1:shift2) = vmat;
        temp_padded_RF_BFRow = temp_padded_RF(BF_dim2,:);
        
        %Compute sparseness
        [SparsenessData, ~] = simple_compute_tuning_sparseness(vmat, IsRF, StimInfo, 0, 'RF', PeakData);
        temp_padded_sparseness(1,shift1:shift2) = SparsenessData.V2_Sparseness{:}; %V2 because we want to know what V2 sparseness/bandwidth is while RF is being centered in the V1 dimension

        %Compute bandwidth at each row
        temp_bandwidth = nan(1,length(V2));
        temp_bandwidth_est = nan(1,length(V2));
        temp_bandwidth_padded = nan(1,length(V2));
        temp_bandwidth_padded_est = nan(1,length(V2));
        if recompute_bandwidth
            for i = BF_dim2%:length(temp_bandwidth)
                %Flip sign for inhibitory RFs
                if ~Filters(F).PLOTBINARY && (strcmp(Filters(F).RESPONSETYPE, 'suppressed') || strcmp(Filters(F).RFTYPE, 'inhibitory'))
                    filt_vmat = -vmat; %Flip sign before computing inhibitory bandwidths
                    filt_padded = -temp_padded_RF(i,:);
                else
                    filt_vmat = vmat;
                    filt_padded = temp_padded_RF(i,:);
                end
                
                %PADDED
                best_ind = zeros(1,size(temp_padded_RF,2));
                best_ind(centerColumn) = 1;
                paddedLabel = (-(size(temp_padded_RF,2)-1)/2):1:((size(temp_padded_RF,2)-1)/2);
                FitTable = tak_fit_gaussian(paddedLabel, filt_padded, best_ind, 'LogTransform', 0, 'PlotFigure', 0);
                FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit
                if ~isempty(FitTable)
                    temp_bandwidth_padded(i) = FitTable.Width;
                    temp_bandwidth_padded_est(i) = FitTable.EstimatedWidth;
                end
                
                %NOT PADDED
                best_ind = zeros(1,length(V1));
                best_ind(BF_dim1) = 1;
                FitTable = tak_fit_gaussian(1:length(best_ind), filt_vmat(BF_dim2,:), best_ind, 'LogTransform', 0, 'PlotFigure', 0);

                FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit
                if ~isempty(FitTable)
                    temp_bandwidth(i) = FitTable.Width;
                    temp_bandwidth_est(i) = FitTable.EstimatedWidth;
                end
            end
        end
        temp_bandwidth_BF = temp_bandwidth(BF_dim2);
        temp_bandwidth_est_BF = temp_bandwidth_est(BF_dim2);
        temp_bandwidth_padded_BF = temp_bandwidth_padded(BF_dim2);
        temp_bandwidth_padded_est_BF = temp_bandwidth_padded_est(BF_dim2);
        
        %Compute bandwidth of average padded RF
        temp_bandwidth_avg = nan;
        temp_bandwidth_avg_est = nan;
        if recompute_bandwidth
            %Flip sign for inhibitory RFs
            if ~Filters(F).PLOTBINARY && (strcmp(Filters(F).RESPONSETYPE, 'suppressed') || strcmp(Filters(F).RFTYPE, 'inhibitory'))
                filt_padded = -temp_padded_RF;
            else
                filt_padded = temp_padded_RF;
            end
                
            %Average RF
            filt_padded_avg = mean(filt_padded,1,'omitnan');
            
            %Filter
            best_ind = zeros(1,size(temp_padded_RF,2));
            best_ind(centerColumn) = 1;
            paddedLabel = (-(size(temp_padded_RF,2)-1)/2):1:((size(temp_padded_RF,2)-1)/2);
            FitTable = tak_fit_gaussian(paddedLabel, filt_padded_avg, best_ind, 'LogTransform', 0, 'PlotFigure', 0);
            FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit
            if ~isempty(FitTable)
                temp_bandwidth_avg = FitTable.Width;
                temp_bandwidth_avg_est = FitTable.EstimatedWidth;
            end
        end
        
        %Concatenate RF to store
        Padded_RF =  cat(3,Padded_RF,temp_padded_RF); %3D concatenation
        Padded_RF_BFRow = [Padded_RF_BFRow; temp_padded_RF_BFRow]; %2D concatenation
        Padded_Sparseness = [Padded_Sparseness; temp_padded_sparseness]; %2D concatenation
        Bandwidth = [Bandwidth; temp_bandwidth]; %2D concatenation
        Bandwidth_BF = [Bandwidth_BF; temp_bandwidth_BF];
    
        %Save data to do stats
        temp_save = mean(temp_padded_RF,1,'omitnan');
        for i = 1:length(temp_save); Summary.(['Norm' num2str(i)])(c) = temp_save(i); end
        for i = 1:length(temp_save); Summary.(['Norm_BFRow' num2str(i)])(c) = temp_padded_RF_BFRow(i); end
        for i = 1:length(temp_save); Summary.(['Norm_Sparseness' num2str(i)])(c) = temp_padded_sparseness(i); end
        %for i = 1:length(temp_save); Summary.(['Norm_Bandwidth' num2str(i)])(c) = Bandwidth(i); end
        Summary.('Norm_Bandwidth_BF')(c) = temp_bandwidth_BF;
        Summary.('Norm_Bandwidth_BF_est')(c) = temp_bandwidth_est_BF;
        Summary.('Norm_Bandwidth_BF_padded')(c) = temp_bandwidth_padded_BF;
        Summary.('Norm_Bandwidth_BF_padded_est')(c) = temp_bandwidth_padded_est_BF;
        Summary.('Norm_Bandwidth_Avg_padded')(c) = temp_bandwidth_avg;
        Summary.('Norm_Bandwidth_Avg_padded_est')(c) = temp_bandwidth_avg_est;
    end

    %Split and store data by Group
    for G = 1:length(GroupList)
        Filters(F).Padded_RF{G} = Padded_RF(:,:,strcmp(Groups, GroupList{G}));
        Filters(F).Padded_RF_BFRow{G} = Padded_RF_BFRow(strcmp(Groups, GroupList{G}),:);
        Filters(F).Padded_Sparseness{G} = Padded_Sparseness(strcmp(Groups, GroupList{G}),:);
        Filters(F).Bandwidth{G} = Bandwidth(strcmp(Groups, GroupList{G}),:);
        Filters(F).Bandwidth_BF{G} = Bandwidth_BF(strcmp(Groups, GroupList{G}));
    end
    
    %Save data for stats
    if save_stats_table
        Filters(F).DATA.Summary = Summary;
        writetable(Summary,['Tuning_' Filters(F).TITLE{:} '_Summary.csv'],'Delimiter',',');
    end
end

%% Plot figures - HEAT MAPS
%1st row: Average of all RFs - no additional normalization
%2nd row: Line graph version of above
%3rd row: Average of all RFs but each row is normalized to the center value
%4th row: Line graph version of above

for F = 1:length(Filters)
    figure; t = tiledlayout(3,length(GroupList));

    for G = 1:length(GroupList)
        if isempty(Filters(F).Padded_RF)
            continue;
        end
        
        Padded_RF = Filters(F).Padded_RF{G};
        avg_padded_RF = squeeze(mean(Padded_RF,3,'omitnan'));

        for f = 1:3
            if f == 1
                nexttile(G)
                imagesc(avg_padded_RF);
                caxis([0 0.3])
            elseif f == 2
                nexttile(G+length(GroupList))
                surf(flipud(avg_padded_RF));
                %shading interp
            elseif f == 3
                nexttile(G+length(GroupList)*2)
                contour(flipud(avg_padded_RF));
            end
            caxis([0 0.3])
            title(GroupList(G))
            if Filters(F).PLOTBINARY
                colormap('bone')
            else
                colormap(bluewhitered)
            end
            colorbar
            set(gca, 'XTick', 1:size(avg_padded_RF,2))
            set(gca, 'YTick', 1:size(avg_padded_RF,1))
            paddedLabel = (-(size(avg_padded_RF,2)-1)/2):1:((size(avg_padded_RF,2)-1)/2);
            set(gca, 'XTickLabel', num2str(paddedLabel'))
            if f == 1
                set(gca, 'YTickLabel', num2str(V2))
            elseif f == 2
                set(gca, 'YTickLabel', num2str(flipud(V2)))
            end
            xlabel(['Relative ' Parameters{1}])
            ylabel(Parameters{2})
        end
    end
    
    sgtitle(Filters(F).TITLE)
end

%% PLOT FIGURES - PVALUE PLOTS
%1st row: pvalue of every position in the padded nanmat compared to the control group
%2nd row: -1 if position is sig. less than control, 1 if position is sig. greater

pvalue = 0.05; %alpha level
controlInd = strcmp(GroupList, control_group);
expGroups = GroupList(~controlInd);
                
for F = 1:length(Filters)
    figure; tiledlayout(2,length(expGroups));
    
    for G = 1:length(expGroups)
        if isempty(Filters(F).Padded_RF)
            continue;
        end
        
        %Get Padded_RF of experimental and control groups
        expInd = strcmp(GroupList,expGroups{G});
        Exp_Padded_RF = Filters(F).Padded_RF{expInd};
        Ctrl_Padded_RF = Filters(F).Padded_RF{controlInd};

        %Initialize
        pvalueMat = nan(size(padded_nanmat));
        signMat = ones(size(padded_nanmat));
        
        %Loop through every position in the Padded RFs
        for i = 1:numel(pvalueMat)
            [sub1,sub2] = ind2sub(size(pvalueMat),i);
            controlVals = squeeze(Exp_Padded_RF(sub1,sub2,:));
            expVals = squeeze(Ctrl_Padded_RF(sub1,sub2,:));

            %Remove missing data
            if Filters(F).PLOTBINARY == 0
                %In non-binary data, 0s and NaNs represent parts of the padded matrices that had no data
                controlVals(controlVals == 0) = [];
                controlVals(isnan(controlVals)) = [];
                expVals(expVals == 0) = [];
                expVals(isnan(expVals)) = [];
            end
            
            %If there is any data left to compare, do a t-test on the two distributions and record sign
            if ~isempty(controlVals) && ~isempty(expVals)
                [~, pvalueMat(i), ~] = ttest2(controlVals, expVals);
                if mean(expVals,'omitnan') < mean(controlVals,'omitnan')
                    signMat(i) = -1;
                end
            end
        end
        
        nexttile(G)
        imagesc(pvalueMat);
        title(expGroups{G})
        subtitle('pvalue')
        colorbar
        set(gca, 'XTick', 1:size(avg_padded_RF,2))
        set(gca, 'YTick', 1:size(avg_padded_RF,1))
        paddedLabel = (-(size(avg_padded_RF,2)-1)/2):1:((size(avg_padded_RF,2)-1)/2);
        set(gca, 'XTickLabel', num2str(paddedLabel'))
        set(gca, 'YTickLabel', num2str(V2))
        xlabel(['Relative ' Parameters{1}])
        ylabel(Parameters{2})
            
        nexttile(G+length(expGroups))
        binary_pvalue = pvalueMat < pvalue;
        binary_pvalue = binary_pvalue.*signMat;
        imagesc(binary_pvalue)
        subtitle(['pvalue < ' num2str(pvalue)])
        caxis([-1 1])
        colorbar
        set(gca, 'XTick', 1:size(avg_padded_RF,2))
        set(gca, 'YTick', 1:size(avg_padded_RF,1))
        paddedLabel = (-(size(avg_padded_RF,2)-1)/2):1:((size(avg_padded_RF,2)-1)/2);
        set(gca, 'XTickLabel', num2str(paddedLabel'))
        set(gca, 'YTickLabel', num2str(V2))
        xlabel(['Relative ' Parameters{1}])
        ylabel(Parameters{2})
    end
    
    sgtitle(Filters(F).TITLE)
end

%% Plot figures - LINE PLOTS: ROW AND COLUMN AVERAGES AND SPARSENESS
%1st row left: average and SEM of all padded mat rows averaged together
%2nd row left: above graph with each row normalized to the center value
%3rd row left: above graph with each row normalized to the center value but for each cell individually before averaging
%4th row left: average and SEM of BF row only
%1st row right: average and SEM of all padded mat columns averaged together
%2nd row right: above graph with each column normalized to the max value
%3rd row right: above graph with each column normalized to the max value for each cell individually before averaging
%4th row right: average sparseness
 
for F = 1:length(Filters)
    figure; tiledlayout(4,2)
    
    for G = 1:length(GroupList)
        if isempty(Filters(F).Padded_RF)
            continue;
        end
        
        Padded_RF = Filters(F).Padded_RF{G};

        %~~~~~~~~ROWS~~~~~~~~
        
        %Average all rows of each padded RF
        nexttile(1); hold on
        tempRows = squeeze(mean(Padded_RF(:,:,:),1,'omitnan'));
        avgRows = mean(tempRows,2,'omitnan');
        semRows = std(tempRows,[],2,'omitnan')./sqrt(size(tempRows,2)-1);
        errorbar(avgRows, semRows);
        set(gca, 'XTick', 1:size(Padded_RF,2))
        paddedLabel = (-(size(avg_padded_RF,2)-1)/2):1:((size(avg_padded_RF,2)-1)/2);
        set(gca, 'XTickLabel', num2str(paddedLabel'))
        ylabel('Response')
        title('All rows averaged')
        
        %Normalize each row to BF column
        %METHOD 1:
        nexttile(3); hold on
        avgRows = avgRows./avgRows(centerColumn);
        plot(avgRows);
        set(gca, 'XTick', 1:size(Padded_RF,2))
        paddedLabel = (-(size(avg_padded_RF,2)-1)/2):1:((size(avg_padded_RF,2)-1)/2);
        set(gca, 'XTickLabel', num2str(paddedLabel'))
        ylabel('Response')
        title('All rows averaged & renormalized to BF')
        
        %METHOD 2:
        nexttile(5); hold on
        tempRows_norm = tempRows./tempRows(centerColumn,:);
        avgRows = mean(tempRows_norm,2,'omitnan');
        semRows = std(tempRows_norm,[],2,'omitnan')./sqrt(size(tempRows_norm,2)-1);
        errorbar(avgRows, semRows);
        set(gca, 'XTick', 1:size(Padded_RF,2))
        paddedLabel = (-(size(avg_padded_RF,2)-1)/2):1:((size(avg_padded_RF,2)-1)/2);
        set(gca, 'XTickLabel', num2str(paddedLabel'))
        ylabel('Response')
        title('All rows averaged & renormalized to BF (each cell individually)')
        
        %Plot BF row only
        nexttile(7); hold on
        Padded_RF_BFRow = Filters(F).Padded_RF_BFRow{G};
        avgRows = mean(Padded_RF_BFRow,1,'omitnan');
        semRows = std(Padded_RF_BFRow,[],1,'omitnan')./sqrt(size(Padded_RF_BFRow,1)-1);
        errorbar(avgRows, semRows);
        set(gca, 'XTick', 1:size(Padded_RF,2))
        paddedLabel = (-(size(avg_padded_RF,2)-1)/2):1:((size(avg_padded_RF,2)-1)/2);
        set(gca, 'XTickLabel', num2str(paddedLabel'))
        xlabel(['Relative ' Parameters{1}])
        ylabel('Response')
        title('BF row')
        
        %~~~~~~~~COLUMNS~~~~~~~~
        
        %Take center column from all padded RFs and average together
        nexttile(2); hold on
        tempColumns = squeeze(Padded_RF(:,centerColumn,:));
        avgColumns = mean(tempColumns,2,'omitnan');
        semColumns = std(tempColumns,[],2,'omitnan')./sqrt(size(tempColumns,2)-1);
        errorbar(avgColumns, semColumns);
        set(gca, 'XTick', 1:size(Padded_RF,1))
        set(gca, 'XTickLabel', num2str(V2))
        ylabel('Response')
        title('All columns averaged')
        
        %Normalize each column to max value
        %METHOD 1
        nexttile(4); hold on
        avgColumns = avgColumns./max(avgColumns);
        plot(avgColumns);
        set(gca, 'XTick', 1:size(Padded_RF,1))
        set(gca, 'XTickLabel', num2str(V2))
        ylabel('Response')
        title('All columns averaged & renormalized to max value')
        
        %METHOD 2
        nexttile(6); hold on
        tempColumns_norm = tempColumns./max(tempColumns,[],1);
        avgColumns = mean(tempColumns_norm,2,'omitnan');
        semColumns = std(tempColumns_norm,[],2,'omitnan')./sqrt(size(tempColumns_norm,2)-1);
        errorbar(avgColumns, semColumns);
        set(gca, 'XTick', 1:size(Padded_RF,1))
        set(gca, 'XTickLabel', num2str(V2))
        xlabel(Parameters{2})
        ylabel('Response')
        title('All columns averaged & renormalized to max value (each cell individually)')
        
        %~~~~~~~~SPARSENESS~~~~~~~~
        
        nexttile(8); hold on
        Padded_Sparseness = Filters(F).Padded_Sparseness{G};
        avgRows = mean(Padded_Sparseness,1,'omitnan');
        semRows = std(Padded_Sparseness,[],1,'omitnan')./sqrt(size(Padded_Sparseness,1)-1);
        errorbar(avgRows, semRows);
        set(gca, 'XTick', 1:size(Padded_RF,2))
        paddedLabel = (-(size(avg_padded_RF,2)-1)/2):1:((size(avg_padded_RF,2)-1)/2);
        set(gca, 'XTickLabel', num2str(paddedLabel'))
        xlabel(['Relative ' Parameters{1}])
        ylabel('Sparseness')
        title('Sparseness')
        
    end
    
    %Add legends
    for f = 1:8
        nexttile(f)
        legend(GroupList)
    end
    
    sgtitle(Filters(F).TITLE)
end

%% Plot figures - COMPUTE BANDWIDTHS FROM LINE GRAPHS

% 
% 
% for F = 1:length(Filters)
%     
%     if strcmp(Filters(F).RFTYPE, "inhibitory") || strcmp(Filters(F).RESPONSETYPE, "suppressed")
%         warning('Bandwidth computations must be adjusted for inhibitory RFs (flip sign before computing');
%         continue;
%     end
%         
%     figure; tiledlayout(1,3)
%     
%     Filters(F).BW1 = cell(1,length(GroupList));
%     Filters(F).BW2 = cell(1,length(GroupList));
%     Filters(F).BW3 = cell(1,length(GroupList));
% 
%     for G = 1:length(GroupList)
%         if isempty(Filters(F).Padded_RF)
%             continue;
%         end
%         
%         Padded_RF = Filters(F).Padded_RF{G};
% 
%         %~~~~~~~~ROWS~~~~~~~~
%         
%         %Average all rows of each padded RF
%         tempRows = squeeze(mean(Padded_RF(:,:,:),1,'omitnan'));
%         
%         %Scatterplot of the bandwidths
%         estimatedBW = nan(1,size(tempRows,2));
%         for i = 1:size(tempRows,2) %estimate bandwidth using gaussian fit for each value in tempRows
%                 
%             %TODO: should bandwidths be centered on centercolumn?????
%             best_ind = zeros(1,size(tempRows,1));
%             best_ind(centerColumn) = 1;
%             FitTable = tak_fit_gaussian(1:size(tempRows,1), tempRows(:,i), [], 'LogTransform', 0, 'PlotFigure', 0);
%             
%             FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit
%             if ~isempty(FitTable)
%                 estimatedBW(i) = FitTable.Width;
%             end
%         end
%         Filters(F).BW1{G} = estimatedBW;
%         
%         %Plot
%         nexttile(1); hold on
%         scatter(zeros(size(estimatedBW))+G, estimatedBW, 50, 'jitter', 'on')
%         line([G-0.25 G+0.25], [median(estimatedBW,'omitnan') median(estimatedBW,'omitnan')], 'Color', 'k', 'Linewidth', 2)
%         set(gca, 'XTick', 1:length(GroupList))
%         set(gca, 'XTickLabel', GroupList)
%         ylabel('Bandwidth')
%         title('All rows averaged')
%         
%         
%         %Normalize each row to BF column
%         %METHOD 2:
%         tempRows_norm = tempRows./tempRows(centerColumn,:);
%         
%         %Scatterplot of the bandwidths
%         estimatedBW = nan(1,size(tempRows,2));
%         for i = 1:size(tempRows,2) %estimate bandwidth using gaussian fit for each value in tempRows
%                 
%             %TODO: should bandwidths be centered on centercolumn?????
%             best_ind = zeros(1,size(tempRows,1));
%             best_ind(centerColumn) = 1;
%             FitTable = tak_fit_gaussian(1:size(tempRows,1), tempRows_norm(:,i), [], 'LogTransform', 0, 'PlotFigure', 0);
%             
%             FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit
%             if ~isempty(FitTable)
%                 estimatedBW(i) = FitTable.Width;
%             end
%         end
%         Filters(F).BW2{G} = estimatedBW;
%         
%         %Plot
%         nexttile(2); hold on
%         scatter(zeros(size(estimatedBW))+G, estimatedBW, 50, 'jitter', 'on')
%         line([G-0.25 G+0.25], [median(estimatedBW,'omitnan') median(estimatedBW,'omitnan')], 'Color', 'k', 'Linewidth', 2)
%         set(gca, 'XTick', 1:length(GroupList))
%         set(gca, 'XTickLabel', GroupList)
%         ylabel('Bandwidth')
%         title('All rows averaged & renormalized to BF (each cell individually)')
%         
%         %Plot BF row only
%         Padded_RF_BFRow = Filters(F).Padded_RF_BFRow{G};
% 
%         %Scatterplot of the bandwidths
%         estimatedBW = nan(1,size(tempRows,2));
%         for i = 1:size(tempRows,2) %estimate bandwidth using gaussian fit for each value in tempRows
%                 
%             %TODO: should bandwidths be centered on centercolumn?????
%             best_ind = zeros(1,size(tempRows,1));
%             best_ind(centerColumn) = 1;
%             FitTable = tak_fit_gaussian(1:size(tempRows,1), Padded_RF_BFRow(i,:), [], 'LogTransform', 0, 'PlotFigure', 0);
%             
%             FitTable = FitTable(FitTable.BestFit == 1,:); %Only keep best fit
%             if ~isempty(FitTable)
%                 estimatedBW(i) = FitTable.Width;
%             end
%         end
%         Filters(F).BW3{G} = estimatedBW;
%         
%         %Plot
%         nexttile(3); hold on
%         scatter(zeros(size(estimatedBW))+G, estimatedBW, 50, 'jitter', 'on')
%         line([G-0.25 G+0.25], [median(estimatedBW,'omitnan') median(estimatedBW,'omitnan')], 'Color', 'k', 'Linewidth', 2)
%         set(gca, 'XTick', 1:length(GroupList))
%         set(gca, 'XTickLabel', GroupList)
%         ylabel('Bandwidth')
%         title('BF row')
%         
%     end
%     
%     sgtitle(Filters(F).TITLE)
% end

%% Save figures

if save_figures
    save_plots;
end