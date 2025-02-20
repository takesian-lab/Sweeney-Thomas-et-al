function [BehaviorData, fig1] = simple_compute_tuning_maryse_behavior(isResponsive, PeakData, StimInfo, plot_figures)
% Compute responses to Maryse behavior stim for a single cell and plot figures
%
% Argument(s): 
%   isResponsive - array of 0s and 1s indicating whether the
%   cell is significantly responsive to that stim combination
%   PeakData - from simple_check_if_resopnsive
%   StimInfo - from simple_prepare_variables
%   plot_figures - 0 or 1 to plot figures
%
% Returns:
%   BehaviorData table:
%     - BestOctave = best octave difference according to highest amplitude response
%     - BestOctave_fit = best octave difference according to gaussian estimate
%     - fit_halfwidth = bandwidth of octave difference preference
%     - MI = modulation index (Resp_alt - Resp_rep / Resp_alt + Resp_rep)
%   fig1 - Gaussian fit + Neurometric curve
% 
% Version History:
%  - V1 = current version made for simple_extract_data
%
% TODO: Find best starting values for fitting
% Search 'TODO'
%% Initial parameters

%Figures will return as [] if plot_figures is 0
fig1 = [];

%BehaviorData will return as NaN if there are no isResponsive stim conditions
BehaviorData = table;
BehaviorData.BestOctave = single(nan);
BehaviorData.BestOctave_fit = single(nan);
BehaviorData.fit_halfwidth = single(nan);
BehaviorData.MI = single(nan);

%Check RF for responsive stim conditions:
if ~any(isResponsive)
    return
end

%Variables will return as NaN if fits cannot be made
[BestOctave_fit, fit_halfwidth] = deal(single(nan));

%% Find preferred octave difference [Two ways]

ovs = StimInfo.V1;
V2 = StimInfo.V2;

%For now, skip analysis if I have FP, Hit, etc in V2
if length(V2) > 1
    return;
end

%Remove undetermined stim from PeakData
PeakData(strcmp(PeakData.ResponseType, 'undetermined'),:) = [];

%Which PeakData to use for Response value?
response = PeakData.Peak_AUC;

%Method #1: Use ov difference with greatest response
response_thresholded = response;
response_thresholded(~isResponsive) = nan;
[~, max_ind] = max(response_thresholded);
BestOctave = ovs(max_ind); %Best Octave

%% Method #2: Fit a gaussian across all octave differences

x = double(ovs);
y = double(response);
options = fitoptions('gauss1');
options.Lower = [0 -inf 0]; %CHECK THIS
options.Upper = [inf max(x) max(x)]; %CHECK THIS

if length(y(~isnan(y))) >= 3
    try
        gauss_fit = fit(x(~isnan(y)), y(~isnan(y)), 'gauss1',options);
        x_ovs = ovs(1):0.1:ovs(length(ovs));
        gauss_curve = gauss_fit(x_ovs);

        %fit_RMS = gauss_fit.c1*2; %gauss_RMS_width
        fit_halfwidth = gauss_fit.c1*2.3548; %full width of the gaussian curve at half the maximum 

        [amplitude, index] = max(gauss_fit(x_ovs)); 
        BestOctave_fit = x_ovs(index);

        if plot_figures
            fig1 = figure;
            subplot(2,1,1); hold on
            plot(ovs,response);
            plot(x_ovs, gauss_curve,'r');
            sz_half_BW_ovs = fit_halfwidth/2;
            hw_y = [amplitude/2 amplitude/2];
            hw_x = [BestOctave_fit-(sz_half_BW_ovs) BestOctave_fit+(sz_half_BW_ovs)];
            plot(hw_x,hw_y,'c');
            xlim([ovs(1)-0.05 ovs(end)+0.05])
            set(gca, 'XTick', ovs)
            set(gca, 'XTickLabel', num2str(ovs))
            title(['Gaussian fit: Best Octave = ' num2str(round(BestOctave_fit,1))])
            xlabel('Octaves');
            ylabel('Response');
            legend({'Raw', 'Fit', 'Bandwidth'})
        end
    catch
        warning('Fit not computed for some reason')
    end
else
    %warning('Not enough points to fit behavior tuning')
end

%% Modulation Index
% Resp_alt - Resp_rep / Resp_alt + Resp_rep

if length(y(~isnan(y))) >= 6

    %Check that ovs are ascending
    if ~issorted(ovs,'ascend')
        error('Ovs should be ascending')
    end

    Resp_rep = mean(response(1:3), 'omitnan'); %Response to the repeating tones (average 3 smallest ov diff)
    Resp_alt = mean(response(end-2:end), 'omitnan'); %Response to the alternating tones (average 3 largest ov diff)

    MI = (Resp_alt - Resp_rep)/(Resp_alt + Resp_rep);

    %Plot neurometric curve
    if plot_figures
        subplot(2,1,2); hold on
        plot(response)
        plot(1:3, zeros(1,3) + Resp_rep, 'r');
        plot(length(response)-2:length(response), zeros(1,3) + Resp_alt, 'r');
        xlim([1-0.05 length(ovs)+0.05])
        set(gca,'XTick',1:length(ovs))
        set(gca,'XTickLabel',num2str(ovs))
        ylabel('Response')
        xlabel('Octaves')
        title(['Modulation index: ' num2str(MI)])
    end
else
    MI = nan;
end

%% Store Data
BehaviorData.BestOctave = single(BestOctave);
BehaviorData.BestOctave_fit = single(BestOctave_fit);
BehaviorData.fit_halfwidth = single(fit_halfwidth);
BehaviorData.MI = single(MI);
