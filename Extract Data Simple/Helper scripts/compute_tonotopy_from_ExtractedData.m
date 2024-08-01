%% compute_tonotopy_from_ExtractedData

%% Setup

plotfigures = 1;
dprimethreshold = 0; %Only include cells with dprime greater than this
datapath = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData\NDNF vs VIP vs PYR April 11 more mice';
datafile = 'ExtractedData_RF_20230412-031712.mat';

cd(datapath)
load(datafile)

%%
SubsetOps = struct;
SubsetOps.Loco                = 'All'; %All, Running, or NotRunning
SubsetOps.ResponseType        = 'activated'; %'' -> no filtering, 'none', 'activated', 'prolonged', 'suppressed', 'excitatory' (use excitatory to combine activated and prolonged)
SubsetOps.RF_Type             = ''; %'' -> no filtering, 'none', 'excitatory', 'inhibitory', 'mixed'
SubsetOps.IsResponsive        = 1; %1 for responsive (and reliable), 0 for not, 2 for both
SubsetOps.FOV                 = ''; %'' -> no filtering, or add FOV name here to only keep that FOV (e.g. L1), put ~ in front to keep everything excep that FOV
SubsetOps.Group               = ''; %'' -> no filtering, or add group name here to only keep that group, put ~ in front to keep everything excep that Group
SubsetOps.Ensemble            = ''; %'' -> no filtering, or add ensemble name for Markpoints Experiments, put ~ in front to keep everything excep that ensemble
SubsetOps.sortbyGCAMP         = 0; %0 for groups, 1 for gcamp, 2 to combine
SubsetOps.sortbyRedCell       = 0; %0 = don't sort, 1 = red cell only, 2 = green cell only
SubsetOps.SuppressOutput      = 1; %0 to print Ops to command line, 1 to suppress

SubsetData = simple_subset_ExtractedData(ExtractedData, SubsetOps);
BlockInfo = SubsetData.BlockInfo; %Need for refImg and conv_factor
CellList = SubsetData.CellList; %Need for x/y position of cells
Summary = SubsetData.Summary; %Need for BF
StimInfo = SubsetData.StimInfo; %Need for stim frequencies

%Set line color range to match frequency spectrum
linecolors = colormap('jet'); %Set colormap for plot lines at the beginning to not affect other plots
close all %closes figure that automatically pops up when you run colormap
freqs = StimInfo.V1*1000;
minfreq = log2(min(freqs));
maxfreq = log2(max(freqs));
range = minfreq:(maxfreq-minfreq)/size(linecolors,1):maxfreq;
range = range(1:size(linecolors,1)); %Make sure arrays are the same size
freqTicks = 0:1/(length(freqs)-1):1; %For colorbar ticks

%% Loop through each block

for b = 1:size(BlockInfo,1)
    
    blockname = BlockInfo.Block(b);
    conv_factor = BlockInfo.conv_factorX(b);
    refImg = BlockInfo.refImg{b};
    refImg = imadjust(int16(refImg)); 
    
    %Find corresponding cells in CellList
    cellInd = strcmp(CellList.Block,blockname);
    if sum(cellInd)<=1; continue; end
    
    %Threshold by dprime (optionally)
    dprime = Summary.dPrime;
    cellInd(dprime < dprimethreshold) = 0;
    if sum(cellInd)<=1; continue; end
    
    currentCells = CellList(cellInd,:);
    BF = Summary.BF(cellInd)*1000; %All cells will have a BF if we set SubsetOps.IsResponsive = 1
    CF = Summary.CF(cellInd)*1000; 
    
    %Get all cell pairs for this block (code from Network_Analysis_Pipeline)
    nCells = size(currentCells,1);
    c1_list = []; sz = nCells; for mm = 1:sz; new_list = ones(sz-mm,1)*mm; c1_list=[c1_list; new_list]; end
    c2_list = []; sz = nCells; for mm = 1:sz; new_list = (mm+1:sz)'; c2_list=[c2_list; new_list]; end
    nPairs = length(c1_list);
    
    %Initialize PairList
    PairList = table;
    PairList.Cell1 = nan(nPairs,1);
    PairList.Cell2 = nan(nPairs,1);
    PairList.Cell1X = nan(nPairs,1);
    PairList.Cell1Y = nan(nPairs,1);
    PairList.Cell2X = nan(nPairs,1);
    PairList.Cell2Y = nan(nPairs,1);
    PairList.Cell1_BF = nan(nPairs,1);
    PairList.Cell2_BF = nan(nPairs,1);
    PairList.Cell1_CF = nan(nPairs,1);
    PairList.Cell2_CF = nan(nPairs,1);
    PairList.DeltaBF = nan(nPairs,1);
    PairList.DeltaCF = nan(nPairs,1);
    PairList.Distance = nan(nPairs,1);

    %Record information for each pair
    for p = 1:size(PairList,1)
        c1 = c1_list(p); %index
        c2 = c2_list(p); %index
        
        PairList.Cell1(p) = currentCells.CellNumber(c1);
        PairList.Cell2(p) = currentCells.CellNumber(c2);
        PairList.Cell1X(p) = currentCells.PositionX(c1);
        PairList.Cell1Y(p) = currentCells.PositionY(c1);
        PairList.Cell2X(p) = currentCells.PositionX(c2);
        PairList.Cell2Y(p) = currentCells.PositionY(c2);
        PairList.Cell1_BF(p) = BF(c1);
        PairList.Cell2_BF(p) = BF(c2);
        PairList.Cell1_CF(p) = CF(c1);
        PairList.Cell2_CF(p) = CF(c2);
        PairList.DeltaBF(p) = log2(BF(c2)) - log2(BF(c1));
        PairList.DeltaCF(p) = log2(CF(c2)) - log2(CF(c1));
    
        cell1XY = [PairList.Cell1X(p), PairList.Cell1Y(p)];
        cell2XY = [PairList.Cell2X(p), PairList.Cell2Y(p)];
        PairList.Distance(p) = norm(cell2XY - cell1XY);
    end
    
    %COMPUTE VECTORS USING NEAREST NEIGHBORING CELL ONLY
    nearestcelldistance = nan(nCells,1);
    nearestcellvectors_BF = nan(nCells,1);
    nearestcellvectors_CF = nan(nCells,1);
    gradients_BF = nan(nCells,4);
    gradients_CF = nan(nCells,4);
    for c = 1:length(nearestcellvectors_BF)
       c1 = currentCells.CellNumber(c);
       
       %Find nearest neighbor
       ind = find(((PairList.Cell1 == c1) + (PairList.Cell2 == c1)) > 0); %find where cell was either pair cell 1 or 2
       [~,ind2] = min(PairList.Distance(ind));
       closestcell = ind(ind2);
       
       %Is c1 cell1 or cell2?
       if isequal(c1, PairList.Cell1(closestcell))
           c1X = PairList.Cell1X(closestcell);
           c1Y = PairList.Cell1Y(closestcell);
           c2X = PairList.Cell2X(closestcell);
           c2Y = PairList.Cell2Y(closestcell);
           c1BF = PairList.Cell1_BF(closestcell);
           c1CF = PairList.Cell1_CF(closestcell);
           c2BF = PairList.Cell2_BF(closestcell);
           c2CF = PairList.Cell2_CF(closestcell);
       else
           c1X = PairList.Cell2X(closestcell);
           c1Y = PairList.Cell2Y(closestcell);
           c2X = PairList.Cell1X(closestcell);
           c2Y = PairList.Cell1Y(closestcell);
           c1BF = PairList.Cell2_BF(closestcell);
           c1CF = PairList.Cell2_CF(closestcell);
           c2BF = PairList.Cell1_BF(closestcell);
           c2CF = PairList.Cell1_CF(closestcell);
       end
       
       %Compute vector
       D = PairList.Distance(closestcell);
       BFdiff = log2(c2BF) - log2(c1BF);
       CFdiff = log2(c2CF) - log2(c1CF);
       diffx = c2X - c1X;
       diffy = c2Y - c1Y;
       vectorsum_BF = (abs(BFdiff))*([diffx,diffy]);
       vectorsum_CF = (abs(CFdiff))*([diffx,diffy]);
       nearestcellvectors_BF(c) = norm(vectorsum_BF);
       nearestcellvectors_CF(c) = norm(vectorsum_CF);
       nearestcelldistance(c) = D;
       
       gradients_BF(c,:) = [c1X, c1Y, vectorsum_BF(1), vectorsum_BF(2)];
       gradients_CF(c,:) = [c1X, c1Y, vectorsum_CF(1), vectorsum_CF(2)];
    end
    

    %COMPUTE HETEROGENEITY OF LOCAL TUNING
    %Variation in local BF tuning was measured from neuropil- corrected two-photon imaging data.
    %Within a given field of view, all somatic ROIs were identified within a 50 μm radius of the
    %reference cell. Provided that a minimum of 5 cells were identified within this area, the median
    %BF was computed across all cells within this local neighborhood. The absolute value of the BF
    %difference for each cell versus the neighborhood median was calculated before repeating the process
    %with a different reference neuron. The interquartile range of this BF distribution was operationally
    %defined as the local BF heterogeneity [Romero et al 2018]
    
    cellMeanDistance = nan(nCells,1);
    cellMeanTuning_BF = nan(nCells,1);
    cellMeanTuning_CF = nan(nCells,1);
    for c = 1:nCells
        c1 = currentCells.CellNumber(c);
        
        %Find 5 closest cells to average together for heterogeneity of tuning
        [~, sortInd] = sort(PairList.Distance, 'ascend');
        sortedPairList = PairList(sortInd,:);
        ind = find(((sortedPairList.Cell1 == c1) + (sortedPairList.Cell2 == c1)) > 0); %find where cell was either pair cell 1 or 2
        if length(ind)<4
            continue;
        end
        closest_ind = ind(1:4);
        closestPairList = sortedPairList(closest_ind,:);
        
        %Get BFs of all pairs
        closestBF = nan(1,4);
        closestCF = nan(1,4);
        for cc = 1:size(closestPairList,1)
            if closestPairList.Cell1 == c1
                closestBF(cc) = log2(closestPairList.Cell2_BF(cc));
                closestCF(cc) = log2(closestPairList.Cell2_CF(cc));
            else
                closestBF(cc) = log2(closestPairList.Cell1_BF(cc));
                closestCF(cc) = log2(closestPairList.Cell1_CF(cc));
            end
        end
        %Concatenate with BF of current cell
        closestBF = [log2(BF(c)), closestBF];
        closestCF = [log2(CF(c)), closestCF];
        
        %Find median
        BFmedian = median(closestBF,'omitnan');
        CFmedian = median(closestCF,'omitnan');
        cellMeanTuning_BF(c) = mean(abs(closestBF-BFmedian),'omitnan');
        cellMeanTuning_CF(c) = mean(abs(closestCF-CFmedian),'omitnan');
        cellMeanDistance(c) = mean(closestPairList.Distance,'omitnan');
    end
    
    %Record cell information
    CellList.CellVector_BF(cellInd) = nearestcellvectors_BF;
    CellList.CellVector_CF(cellInd) = nearestcellvectors_CF;
    CellList.NearestCellDistance(cellInd) = nearestcelldistance;
    CellList.MeanLocalDistance(cellInd) = cellMeanDistance; %Distance between each cell and its closest 5 neighbors
    CellList.MeanLocalTuning_BF(cellInd) = cellMeanTuning_BF; %DeltaBF between each cell and its closest 5 neighbors
    CellList.MeanLocalTuning_CF(cellInd) = cellMeanTuning_CF; %DeltaCF between each cell and its closest 5 neighbors
    
    %Record block information
    BlockInfo.Vector_BF(b) = mean(nearestcellvectors_BF,'omitnan');
    BlockInfo.Vector_BF_STD(b) = std(nearestcellvectors_BF,'omitnan');
    BlockInfo.Vector_CF(b) = mean(nearestcellvectors_CF,'omitnan');
    BlockInfo.Vector_CF_STD(b) = std(nearestcellvectors_CF,'omitnan');
    BlockInfo.Gradients_BF{b} = gradients_BF;
    BlockInfo.Gradients_CF{b} = gradients_CF;
    BlockInfo.DistanceAllCells(b) = mean(PairList.Distance,'omitnan');
    BlockInfo.DistanceAllCellsSTD(b) = std(PairList.Distance,'omitnan');
    
    %% FIGURE
    if plotfigures

        figure('units','normalized','outerposition',[0 0 1 1])

        %BACKGROUNDS
        %Plot suite2p image
        for i = 2:3%1:3
            ax(i) = subplot(3,3,i); hold on
            imagesc(refImg)
            axis square
            xlim([0 512])
            ylim([0 512])
            colormap(ax(i), 'bone')
            set(ax(i),'YDir','reverse')
            axis off
            if i == 1
                plot(CRlineX, fliplr(CRlineY), 'w')
            end
        end

        %Plot ROI masks only to get separate colorbar
        for i = 5:6
            ax(i) = subplot(3,3,i); hold on
            axis square
            xlim([0 512])
            ylim([0 512])
            colormap(ax(i),'jet')
            set(ax(i),'YDir','reverse')
            axis off
        end

        %CELLS

        %Plot ROI masks and vector arrows
        
        f = [2 3 5 6];
        for i = 1:length(f)
            switch f(i)
                case 2
                    plottitle = 'BF';
                    dataToPlot = BF;
                    vectorsToPlot = gradients_BF;
                    plotcbar = 0;
                case 3
                    plottitle = 'CF';
                    dataToPlot = CF;
                    vectorsToPlot = gradients_CF;
                    plotcbar = 0;
                case 5
                    plottitle = 'BF';
                    dataToPlot = BF;
                    vectorsToPlot = gradients_BF;
                    plotcbar = 1;
                case 6
                    plottitle = 'CF';
                    dataToPlot = CF;
                    vectorsToPlot = gradients_CF;
                    plotcbar = 1;
            end

            colorToPlot = nan(size(dataToPlot));
            for ii = 1:length(dataToPlot)
                [~, colorToPlot(ii)] = min(abs(range - log2(dataToPlot(ii))));
            end

            ax(f(i)) = subplot(3,3,f(i));
            ax(f(i)).ColorOrder = linecolors(colorToPlot,:);
            for c = 1:nCells
                xcirc = currentCells.xcirc{c};
                ycirc = currentCells.ycirc{c};

                plot(xcirc,ycirc,'Linewidth', 1.5);
%                 if ~plotcbar
%                     text(max(xcirc),max(ycirc),num2str(currentCells.CellNumber(c)), 'Color', 'w');
%                 else
%                     %text(max(xcirc),max(ycirc),num2str(currentCells.CellNumber(c)), 'Color', 'k');
%                 end
            end

            if plotcbar
                quiver(vectorsToPlot(:,1),vectorsToPlot(:,2),vectorsToPlot(:,3),vectorsToPlot(:,4),'k')
                cbar = colorbar('EastOutside', 'YTick', freqTicks, 'TickLabels', freqs);
            end

            title(plottitle)
        end

        %VECTORS
        f = [8 9];
        for i = 1:length(f)
            switch f(i)
                case 8
                    vectorsToPlot = gradients_BF;
                case 9
                    vectorsToPlot = gradients_CF;
            end
            ax(f(i)) = subplot(3,3,f(i)); hold on
            xVector = sum(vectorsToPlot(:,3));
            yVector = sum(vectorsToPlot(:,4));
            quiver(xVector,yVector,0,'r');
            xlim([0.85 1.15])
            ylim([0.85, 1.15])
            hline(1,'k')
            vline(1,'k')
            axis square
            axis off
        end

        %TUNING HETEROGENEITY
        ax(4) = subplot(3,3,4); hold on
        scatter(cellMeanDistance, cellMeanTuning_BF)
        xlabel('Mean distance')
        ylabel('BF Tuning Heterogeneity')
        title('BF Tuning heterogeneity between 5 closest cells')

        ax(7) = subplot(3,3,7); hold on
        scatter(cellMeanDistance, cellMeanTuning_CF)
        xlabel('Mean distance')
        ylabel('CF Tuning Heterogeneity')
        title('CF Tuning heterogeneity between 5 closest cells')
        
        sgtitle(regexprep(blockname,'_',' ','emptymatch'))
    end

end

%% Plot group differences

Groups = unique(Summary.Group);

%Distance vs. tuning similarity

figure; hold on

for g = 1:length(Groups)
    subplot(2,length(Groups),g); hold on
    groupInd = strcmp(Summary.Group,Groups{g});
    X = CellList.MeanLocalDistance(groupInd);
    Y = CellList.MeanLocalTuning_BF(groupInd);
    
    scatter(X,Y);
    lsline
    [R,P] = corrcoef(X,Y,'Rows','pairwise');
    subtitle(['Pearsons R ' num2str(round(R(2),3)) ' p ' num2str(round(P(2),3))])
%     %Moving mean
%     [~,ind] = sort(X,'ascend');
%     roll = movmean(Y(ind),10);
%     plot(X(ind),roll,'Linewidth',2,'Color','k')
    title([Groups{g}])
    xlim([0 512])
    ylim([0 1.2])
    %hline(mean(Y,'omitnan'))
    ylabel('BF')
    
    subplot(2,length(Groups),g+length(Groups)); hold on
    groupInd = strcmp(Summary.Group,Groups{g});
    Y = CellList.MeanLocalTuning_CF(groupInd);
    
    scatter(X,Y);
    lsline
    [R,P] = corrcoef(X,Y,'Rows','pairwise');
    subtitle(['Pearsons R ' num2str(round(R(2),3)) ' p ' num2str(round(P(2),3))])
%     %Moving mean
%     [~,ind] = sort(X,'ascend');
%     roll = movmean(Y(ind),10);
%     plot(X(ind),roll,'Linewidth',2,'Color','k')
    title([Groups{g}])
    xlim([0 512])
    ylim([0 1.2])
    %hline(mean(Y,'omitnan'))
    ylabel('CF')
end

sgtitle('Tuning heterogeneity vs. Distance')

%% Plot group differences

Groups = unique(Summary.Group);

%Dprime vs. tuning similarity

figure; hold on

for g = 1:length(Groups)
    subplot(2,length(Groups),g); hold on
    groupInd = strcmp(Summary.Group,Groups{g});
    dprime = Summary.dPrime(groupInd);
    X = dprime;
    Y = CellList.MeanLocalTuning_BF(groupInd);
    
    scatter(X,Y);
    lsline
    [R,P] = corrcoef(X,Y,'Rows','pairwise');
    subtitle(['Pearsons R ' num2str(round(R(2),3)) ' p ' num2str(round(P(2),3))])
%     %Moving mean
%     [~,ind] = sort(X,'ascend');
%     roll = movmean(Y(ind),10);
%     plot(X(ind),roll,'Linewidth',2,'Color','k')
    title([Groups{g}])
    xlim([0 5])
    ylim([0 1.2])
    ylabel('BF')
    
    subplot(2,length(Groups),g+length(Groups)); hold on
    groupInd = strcmp(Summary.Group,Groups{g});
    Y = CellList.MeanLocalTuning_CF(groupInd);
    
    scatter(X,Y);
    lsline
    [R,P] = corrcoef(X,Y,'Rows','pairwise');
    subtitle(['Pearsons R ' num2str(round(R(2),3)) ' p ' num2str(round(P(2),3))])
%     %Moving mean
%     [~,ind] = sort(X,'ascend');
%     roll = movmean(Y(ind),10);
%     plot(X(ind),roll,'Linewidth',2,'Color','k')
    title([Groups{g}])
    xlim([0 5])
    ylim([0 1.2])
    ylabel('CF')
end

sgtitle('Tuning heterogeneity vs. dPrime')

%% Vectors

Groups = unique(Summary.Group);

figure; hold on

subplot(1,5,1); hold on

for g = 1:length(Groups)
    groupInd = strcmp(BlockInfo.Group,Groups{g});
    Y = BlockInfo.Vector_BF(groupInd);
    scatter(zeros(length(Y))+g,Y,[],'jitter', 'off');
    line([g-0.25 g+0.25], [mean(Y,'omitnan') mean(Y,'omitnan')], 'Color', 'k', 'Linewidth',1)
end
set(gca,'Xtick', 1:length(Groups))
set(gca,'Xticklabel', Groups)
xlim([0.5 length(Groups)+0.5])
ylabel('BF Vector')

subplot(1,5,2); hold on

for g = 1:length(Groups)
    groupInd = strcmp(BlockInfo.Group,Groups{g});
    Y = BlockInfo.Vector_CF(groupInd);
    scatter(zeros(length(Y))+g,Y,[],'jitter', 'off');
    line([g-0.25 g+0.25], [mean(Y,'omitnan') mean(Y,'omitnan')], 'Color', 'k', 'Linewidth',1)
end
set(gca,'Xtick', 1:length(Groups))
set(gca,'Xticklabel', Groups)
xlim([0.5 length(Groups)+0.5])
ylabel('CF Vector')

subplot(1,5,3); hold on

for g = 1:length(Groups)
    groupInd = strcmp(BlockInfo.Group,Groups{g});
    Y = BlockInfo.DistanceAllCells(groupInd);
    scatter(zeros(length(Y))+g,Y,[],'jitter', 'off');
    line([g-0.25 g+0.25], [mean(Y,'omitnan') mean(Y,'omitnan')], 'Color', 'k', 'Linewidth',1)
end
set(gca,'Xtick', 1:length(Groups))
set(gca,'Xticklabel', Groups)
xlim([0.5 length(Groups)+0.5])
ylabel('Average cell distance')

subplot(1,5,4); hold on

for g = 1:length(Groups)
    groupInd = strcmp(CellList.Group,Groups{g});
    Y = CellList.MeanLocalTuning_BF(groupInd);
    scatter(zeros(length(Y))+g,Y,[],'jitter', 'off');
    line([g-0.25 g+0.25], [mean(Y,'omitnan') mean(Y,'omitnan')], 'Color', 'k', 'Linewidth',1)
end
set(gca,'Xtick', 1:length(Groups))
set(gca,'Xticklabel', Groups)
xlim([0.5 length(Groups)+0.5])
ylabel('Local deltaBF')

subplot(1,5,5); hold on

for g = 1:length(Groups)
    groupInd = strcmp(CellList.Group,Groups{g});
    Y = CellList.MeanLocalTuning_CF(groupInd);
    scatter(zeros(length(Y))+g,Y,[],'jitter', 'off');
    line([g-0.25 g+0.25], [mean(Y,'omitnan') mean(Y,'omitnan')], 'Color', 'k', 'Linewidth',1)
end
set(gca,'Xtick', 1:length(Groups))
set(gca,'Xticklabel', Groups)
xlim([0.5 length(Groups)+0.5])
ylabel('Local deltaCF')
