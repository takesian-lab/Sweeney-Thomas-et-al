function [fig1] = simple_plot_frequency_tuning(IsResponsive, PeakData, StimInfo)

%% Setup

%Get RF mapping from StimInfo
RF_Map = StimInfo.RF_Map;
nanmat = nan(max(RF_Map));
RF_ind = sub2ind(size(nanmat),RF_Map(:,1),RF_Map(:,2));

%Reconstruct RF maps
[IsRF, RF_AUC, RF_Peak, IsRF_sign] = deal(nanmat);
IsRF(RF_ind) = IsResponsive;

%Which PeakData to use for Response value?
%--> Use AUC, but make suppressed AUCs negative for proper inhibitory/excitatory RF detection
%--> If ALL significant responses are suppressed, flip the signs so that we can measure suppressed RFs

%Generate RF from AUC, assuming excitatory or mixed RF for now
temp_RF_AUC = nan(size(IsResponsive));
temp_RF_Peak = nan(size(IsResponsive));
temp_IsRF_sign = IsResponsive;
for i = 1:length(temp_RF_AUC)
    response_type = PeakData.ResponseType(i);
    if strcmp(response_type, 'suppressed')
        temp_RF_AUC(i) = -PeakData.Trough_AUC(i);
        temp_RF_Peak(i) = PeakData.Trough_3fr(i);
        IsRF_sign(i) = -1;
    else
        temp_RF_AUC(i) = PeakData.Peak_AUC(i);
        temp_RF_Peak(i) = PeakData.Peak_3fr(i);
    end
end

%Assign RF values to RF map
RF_AUC(RF_ind) = temp_RF_AUC;
RF_Peak(RF_ind) = temp_RF_Peak; %For plotting
IsRF_sign(RF_ind) = temp_IsRF_sign;

%convert freqs to octaves
V1 = StimInfo.V1;
V2 = StimInfo.V2;

%Check RF for responsive stim conditions and make sure there is >1 stim condition:
if ~any(any(IsRF)) || isequal(size(IsRF),[1,1])
    return
end

%% Find threshold and BF

%find BF
thresholded_RF = RF_AUC;
thresholded_RF(~IsRF) = nan;
[M,I] = max(thresholded_RF(:));
[BF_row, BF_col] = ind2sub(size(RF_AUC),I);
BF_I = V2(BF_row);
BF = V1(BF_col);

%% Plot RF and FRA

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);

for f = 1:9

    if f == 1
        use_mask = 0;
        use_bluewhitered = 0;
        datatoplot = RF_Peak;
        figtitle = 'Peak';
    elseif f == 2
        use_mask = 1;
        use_bluewhitered = 0;
        datatoplot = RF_Peak;
        figtitle = 'Peak';
    elseif f == 3
        use_mask = 0;
        use_bluewhitered = 0;
        datatoplot = RF_AUC;
        figtitle = 'AUC';
    elseif f == 4
        use_mask = 1;
        use_bluewhitered = 0;
        datatoplot = RF_AUC;
        figtitle = 'AUC';
    elseif f == 5
        use_mask = 0;
        use_bluewhitered = 0;
        datatoplot = IsRF_sign;
        figtitle = 'IsRF Sign';
    elseif f == 6
        use_mask = 0;
        use_bluewhitered = 1;
        datatoplot = RF_Peak;
        figtitle = 'Peak';
    elseif f == 7
        use_mask = 1;
        use_bluewhitered = 1;
        datatoplot = RF_Peak;
        figtitle = 'Peak';
    elseif f == 8
        use_mask = 0;
        use_bluewhitered = 1;
        datatoplot = RF_AUC;
        figtitle = 'AUC';
    elseif f == 9
        use_mask = 1;
        use_bluewhitered = 1;
        datatoplot = RF_AUC;
        figtitle = 'AUC';
    end
    
    subplot(2,5,f); 
    
    RFmask = ones(size(nanmat))*0.1;
    RFmask(IsRF_sign == 1) = 1;
    %RFmask(IsRF == 1) = 1;
    if use_mask
        imagesc(datatoplot, 'AlphaData', RFmask);
    else
        imagesc(datatoplot, 'AlphaData', ones(size(nanmat)));
    end
    title(figtitle)
    xlabel(StimInfo.Parameters{1})
    ylabel(StimInfo.Parameters{2})
    set(gca,'XTick',1:length(V1))
    set(gca,'YTick',1:length(V2))
    set(gca,'XTickLabel',num2str(V1))
    set(gca,'YTickLabel',num2str(V2))
    text(BF_col - 0.25, BF_row, 'BF', 'Color', 'w') %-0.25 centers the text

    if use_bluewhitered
        colormap(gca, bluewhitered(256)); %bluewhitered_TAKlab)
    end
    c = colorbar;
    
    c.Title.String = 'Z';
end

return

caxis([0 4])
