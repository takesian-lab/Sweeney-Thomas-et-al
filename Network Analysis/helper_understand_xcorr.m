%helper_understand_xcorr

load('SampleData.mat')

fs = 30;
maxlag_in_s = 20;
maxlag = floor(fs*maxlag_in_s); 

%Get data
loco_speed = SampleData.loco;
xoff = SampleData.xoff;
yoff = SampleData.yoff;
XY = sqrt(xoff.^2 + yoff.^2);
zdrift = SampleData.zdrift;
timestamp = SampleData.timestamp;

%Zscore traces
loco_speed_zscored = zscore(loco_speed, '', 2);
XY_zscored = zscore(XY, '', 2);
zdrift_zscored = zscore(zdrift, '', 2);

%% Figure 1
figure; hold on;

%Plot locomotor activity
subplot(4,1,1); hold on
title('Locomotor activity')
ylabel('Speed (cm/s)') 
plotMax = max(loco_speed) + 1;
area(timestamp,(loco_speed > 0)*plotMax, 'EdgeColor', 'none', 'Facecolor', [210/255, 248/255, 210/255])
area(timestamp,(loco_speed == 0)*plotMax, 'EdgeColor', 'none', 'Facecolor', [238/255, 144/255, 144/255])
plot(timestamp, loco_speed, 'LineWidth', 0.25, 'Color', 'k')
hline(0.7) %Noise floor of the wheel
ylim([0 plotMax])
xlim([timestamp(1) timestamp(end)])

subplot(4,1,2); hold on
title('X/Y movement')
ylabel('pixels') 
plot(timestamp, xoff, 'LineWidth', 0.25)
plot(timestamp, yoff, 'Linewidth', 0.25)
plot(timestamp, XY, 'Linewidth', 0.25, 'Color', 'k')
xlim([timestamp(1) timestamp(end)])
xlabel('Seconds')
legend({'xoff', 'yoff', 'combined'})

%Perform linear regression between loco and XY trace and compute residuals
p = polyfit(XY_zscored, loco_speed_zscored, 1);
LocoFit = polyval(p, XY_zscored);
LocoResiduals = loco_speed_zscored - LocoFit;
LocoResiduals_zscored = zscore(LocoResiduals,'',2);

subplot(4,1,3); hold on
title('Residuals of Loco + X/Y movement')
ylabel('Z') 
plot(timestamp, XY_zscored, 'LineWidth', 0.25, 'Color', [210/255, 248/255, 210/255])
plot(timestamp, loco_speed_zscored, 'LineWidth', 2, 'Color', [238/255, 144/255, 144/255])
plot(timestamp, LocoResiduals, 'LineWidth', 0.25, 'Color', 'k')
xlim([timestamp(1) timestamp(end)])
xlabel('Seconds')
legend({'Z-scored XY', 'Z-score Loco', 'Residuals'})

subplot(4,1,4); hold on
title('Z movement')
ylabel('Z slices') 
plot(timestamp, zdrift, 'LineWidth', 0.25)
xlim([timestamp(1) timestamp(end)])
xlabel('Seconds')
legend('Z drift')

%% Figure 2
%Understand the gist of xcorr

for f = 1:2
    if f == 1
        loco_speed = SampleData.loco;
        titletag = 'using raw loco data';
    else
        loco_speed = loco_speed_zscored;
        titletag = 'using zscored loco data';
    end
    
    loco_speed_shifted = circshift(loco_speed,60);

    figure; hold on

    subplot(4,3,1); hold on
    plot(timestamp, loco_speed, 'k')
    plot(timestamp, loco_speed, 'k')
    legend({'A', 'B'})
    title('Autocorrelation')
    subplot(4,3,4); hold on
    [cc, lags] = xcorr(loco_speed, loco_speed, maxlag, 'coeff');
    plot(lags, cc)
    [~, peakind] = max(cc);
    ylim([-.2 1])
    vline(0)
    title(['Peak = ' num2str(lags(peakind))])

    subplot(4,3,2); hold on
    plot(timestamp, loco_speed, 'k')
    plot(timestamp, loco_speed_shifted, 'r')
    legend({'A', 'B'})
    title('Trace A precedes B')
    subplot(4,3,5); hold on
    [cc, lags] = xcorr(loco_speed, loco_speed_shifted, maxlag, 'coeff');
    plot(lags, cc)
    [~, peakind] = max(cc);
    ylim([-.2 1])
    vline(0)
    title(['Peak = ' num2str(lags(peakind))])

    subplot(4,3,3); hold on
    plot(timestamp, loco_speed_shifted, 'k')
    plot(timestamp, loco_speed, 'r')
    legend({'A', 'B'})
    title('Trace A follows B')
    subplot(4,3,6); hold on
    [cc, lags] = xcorr(loco_speed_shifted, loco_speed, maxlag, 'coeff');
    plot(lags, cc)
    [~, peakind] = max(cc);
    ylim([-.2 1])
    vline(0)
    title(['Peak = ' num2str(lags(peakind))])

    subplot(4,3,7); hold on
    plot(timestamp, loco_speed, 'k')
    plot(timestamp, loco_speed./2, 'r')
    legend({'A', 'B'})
    title('B is half the magnitude of A')
    subplot(4,3,10); hold on
    [cc, lags] = xcorr(loco_speed, loco_speed./2, maxlag, 'coeff');
    plot(lags, cc)
    [~, peakind] = max(cc);
    ylim([-.2 1])
    vline(0)
    title(['Peak = ' num2str(lags(peakind))])

    subplot(4,3,8); hold on
    plot(timestamp, loco_speed, 'k')
    plot(timestamp, -loco_speed, 'r')
    legend({'A', 'B'})
    title('B is the opposite of A')
    subplot(4,3,11); hold on
    [cc, lags] = xcorr(loco_speed, -loco_speed, maxlag, 'coeff');
    plot(lags, cc)
    [~, peakind] = min(cc);
    ylim([-1 0.2])
    vline(0)
    title(['Min = ' num2str(lags(peakind))])

    sgtitle(['Xcorr results with coeff option and ' titletag])
end

%% Figure 3
%Test zscored and non-zscored loco trace with different normalization methods
% 'none' (default) | 'biased' | 'unbiased' | 'normalized' | 'coeff'

loco_speed = SampleData.loco;

figure; hold on

[cc1,lags] = xcorr(loco_speed_zscored, XY_zscored, maxlag, 'none');
[cc2,lags] = xcorr(loco_speed, XY_zscored, maxlag, 'none');
subplot(6,2,1); plot(lags, cc1); title('Loco zscored - none')
subplot(6,2,2); plot(lags, cc2); title('Loco raw - none')

[cc1,lags] = xcorr(loco_speed_zscored, XY_zscored, maxlag, 'biased');
[cc2,lags] = xcorr(loco_speed, XY_zscored, maxlag, 'biased');
subplot(6,2,3); plot(lags, cc1); title('Loco zscored - biased')
subplot(6,2,4); plot(lags, cc2); title('Loco raw - biased')
subplot(6,2,11); hold on; plot(lags, cc1);

[cc1,lags] = xcorr(loco_speed_zscored, XY_zscored, maxlag, 'unbiased');
[cc2,lags] = xcorr(loco_speed, XY_zscored, maxlag, 'unbiased');
subplot(6,2,5); plot(lags, cc1); title('Loco zscored - unbiased')
subplot(6,2,6); plot(lags, cc2); title('Loco raw - unbiased')
subplot(6,2,11); hold on; plot(lags, cc1);

[cc1,lags] = xcorr(loco_speed_zscored, XY_zscored, maxlag, 'normalized');
[cc2,lags] = xcorr(loco_speed, XY_zscored, maxlag, 'normalized');
subplot(6,2,7); plot(lags, cc1); title('Loco zscored - normalized')
subplot(6,2,8); plot(lags, cc2); title('Loco raw - normalized')
subplot(6,2,11); hold on; plot(lags, cc1);
subplot(6,2,12); hold on; plot(lags, cc2);

[cc1,lags] = xcorr(loco_speed_zscored, XY_zscored, maxlag, 'coeff');
[cc2,lags] = xcorr(loco_speed, XY_zscored, maxlag, 'coeff');
subplot(6,2,9); plot(lags, cc1); title('Loco zscored - coeff')
subplot(6,2,10); plot(lags, cc2); title('Loco raw - coeff')
subplot(6,2,11); hold on; plot(lags, cc1);
legend({'biased', 'unbiased', 'normalized', 'coeff'})
subplot(6,2,12); hold on; plot(lags, cc2);
legend({'normalized', 'coeff'})

sgtitle('Xcorr of loco (zscored or raw) with XY')

figure; hold on
[cc1,lags] = xcorr(loco_speed_zscored, XY_zscored, maxlag, 'coeff');
[cc2,lags] = xcorr(loco_speed, XY_zscored, maxlag, 'coeff');
plot(lags, cc1);
plot(lags, cc2);
legend({'zscored', 'raw'})
sgtitle('Difference between raw and zscored with coeff option')

%% code to understand how XCorr is actually computed

clear all; 

% setup for either sinusoidal wave or cell parameters
z = 0; % if use zscore
data_type = 'sine'; % use simple sine wave or synthetic data (sine or synthetic)

% sine wave
fs = 100; % sampling frequency
sine_time = 2; % length (time) of sine waves

sine1.freq = 2; % freq of sine1 wave
sine1.amp = 2; % amplitude of sine1 wave
sine1.noise = 0.5; % noise of sine1 wave
sine1.dc_shift = 20; % dc_shift of sine2
sine1.phi = 0; % phase shift sine 1

sine2.freq = 2; % freq of sine2 wave
sine2.amp = 2; % amplitude of sine2 wave
sine2.noise = 0.5; % noise of sine2 wave
sine2.dc_shift = 5; % dc_shift of sine2
sine2.phi = 1; % phase shift sine 2

% synthetic calcium data
fs = 30; % sampling rate
syn_time = 120; % total trace time (in seconds)

cell1.transient_freq = 2; % freq of cell 1 transients
cell1.dc_shift = 0; % dc shift of cell 1, default = 0  
cell1.Fs = fs;
cell1.noiseAmount = 3; % noise of cell 1; default = 0.03
cell1.spks = [5 5]; % number of spikes (bursts) per transient

cell2.transient_freq = 2; % freq of cell 1 transients
cell2.dc_shift = 0; % dc shift of cell 1, default = 0  
cell2.Fs = fs;
cell2.noiseAmount = 0.05; % noise of cell 1; default = 0.03
cell2.corrspks = [5 5]; % number of spikes (bursts) per correlated transient
cell2.uncorrspks = [5 5]; % number of spikes per uncorrelated transient
cell2.p_corr = 0.7; % percent of spikes correlated with cell 1
cell2.phase = 0; % phase shift (s) from cell 1

% generate traces
figure; tiledlayout(2,4);

if contains(data_type,'sine')
    % create two time-varying (sinusoidal) signals as test data
    t = 0:(1/fs):sine_time; % time
    x = sine1.amp*sin(2*pi*sine1.freq*t+sine1.phi)+sine1.noise*rand(size(t))+sine1.dc_shift;
    y = sine2.amp*sin(2*pi*sine2.freq*t+sine2.phi)+sine2.noise*rand(size(t))+sine2.dc_shift;
    if z % z_score 
        x=zscore(x);
        y=zscore(y);
    end
end

if contains(data_type,'synthetic')

    t = 0:(1/fs):syn_time; % trace times

    % create a spike train of randomly distributed transients for cell 1
    binned_spike_train = zeros(1,length(t));
    binned_spk_index = randi(length(t),[1 cell1.transient_freq*max(t)]);
    spk_number = randi(cell1.spks,[1 cell1.transient_freq*max(t)]);
    binned_spike_train(binned_spk_index) = spk_number;

    % convolve spikes with calcium transients (GCaMP6f) to create synthetic
    % calcium trace (using Dombeck lab code)
    x = spk2F(binned_spike_train, 'Fs', cell1.Fs, 'noiseAmount', cell1.noiseAmount) + cell1.dc_shift; 
  
    % create a spike train of transients for cell 2, some are correlated
    % with cell 1
    binned_spike_train = zeros(1,length(t));
    corr_transients = round(cell1.transient_freq*max(t)*cell2.p_corr); % number of cell1 and cell 2 correlated transients
    uncorr_transients = round(cell1.transient_freq*max(t))-corr_transients; % number of cell2 transients not correlated to cell 1
    corr_ind = randsample(length(binned_spk_index),corr_transients,false);
    corr_spks = binned_spk_index(corr_ind)+cell2.phase*fs; % spikes correlated with cell #1
    uncorr_spks = randi(length(t),[1 uncorr_transients]); % uncorrelated spks
    corr_spk_number = randi(cell2.corrspks,[1 length(corr_spks)]); % spikes per bin for correlated spikes
    uncorr_spk_number = randi(cell2.uncorrspks,[1 length(uncorr_spks)]); % spikes per bin for uncorrelated spikes
    binned_spike_train(corr_spks) = corr_spk_number;
    binned_spike_train(uncorr_spks) = uncorr_spk_number;

    % convolve spikes with calcium transients (GCaMP6f) to create synthetic
    % calcium trace (using Dombeck lab code)
    y = spk2F(binned_spike_train, 'Fs', cell2.Fs, 'noiseAmount', cell2.noiseAmount) + cell2.dc_shift; 
    
    if z
        x=zscore(x);
        y=zscore(y);
    end
end

% plot both waves
nexttile; plot(t,x,'c','DisplayName', 'X'); 
if contains(data_type,'sine'); xlim([0 5]); else xlim([0 20]); end

hold on; plot(t,y,'m','DisplayName', 'X'); 
xlabel('time(s)'); ylabel('amplitude'); title('test traces'); 

% Perform Matlab's XCorr using different options
% maxlag in seconds
maxlag = 1*fs; % maxlag (seconds * sampling rate)

% perform x corr function
[r,lags] = xcorr(x,y,maxlag,'none');
[r_biased,lags] = xcorr(x,y,maxlag,'biased');
[r_unbiased,lags] = xcorr(x,y,maxlag,'unbiased');
[r_coef,lags] = xcorr(x,y,maxlag,'coef');

% N = length of x and y vectors (same)
N = length(t);

% Compute and plot autocorrelations
[x_auto,lags] = xcorr(x,x,maxlag);
[y_auto,lags] = xcorr(y,y,maxlag);
nexttile; plot(lags,x_auto, 'c','DisplayName', 'X auto corr'); hold on; plot(lags,y_auto, 'm', 'DisplayName', 'Y auto corr'); legend;
xlabel('lags'); ylabel('correlation'); title('X and Y autocorrelations using XCorr'); 

% the following code is how Xcorr is computed (it's simple!) - see the corresponding plots

% pad your x vector with zeros before and after
padded_x = [zeros(maxlag,1); x'; zeros(maxlag,1)];
nexttile(3,[1 2]); plot(padded_x-50, 'c','Linewidth',2,'DisplayName', 'X vector'); hold on;  % plot the x padded vector on the bottom of the figure
for m = 1:maxlag*2+1 % loop across lags
    padded_y = [zeros(m-1,1);y';zeros(2*maxlag-(m-1),1)]; % make shifted vectors of y, padded with zeros
    r_DIY(m) = dot(padded_x,padded_y); % get the dot product for all the lags - sum of the products of all
                                   % timeframes of x with timeframes of each 'lagged y' 
    r_DIY_biased(m) = r_DIY(m)./N; % 'biased' Xcorr is just the correlations divided by the number of points/trace
    r_DIY_unbiased(m) = r_DIY(m)./(N-abs(m-maxlag)); % 'unbiased' XCorr is scaled by how far from the edge of the trace
    r_DIY_coef(m) = r_DIY(m).*(1./sqrt(max(x_auto)*max(y_auto))); % normalized to the sizes of x and y
    plot(padded_y+m*10,'DisplayName', 'lagged Y vectors'); hold on;
end 
xlabel('frames (~time)'); ylabel('amplitude');
title('XCorr is just the sum of padded X times padded Y (dot product) for different lags')

% Plot the matlab xcorr and DIY correlation on the same graph - they are the same! :) 

nexttile; plot(lags./fs,r,'Linewidth', 4, 'DisplayName','Matlab XCorr'); hold on; plot(lags./fs,r_DIY', 'Linewidth', 2, 'DisplayName', 'DYI XCorr'); legend;
xlabel('lags(s)'); ylabel('correlation'); title('XCorr (Matlab vs DIY)');

% plot the matlab xcorr 'biased' and DIY 'biased' correlation on the same graph - they are the same! :)  
nexttile; plot(lags./fs,r_biased,'Linewidth', 4, 'DisplayName','Matlab XCorr Biased'); hold on; plot(lags./fs,r_DIY_biased', 'Linewidth', 2, 'DisplayName', 'DYI XCorr Biased'); legend;
xlabel('lags(s)'); ylabel('correlation'); title('XCorr Biased (Matlab vs DIY)');

% plot the matlab xcorr 'unbiased' and DIY 'unbiased' correlation on the same graph - they are the same! :)  
nexttile; plot(lags./fs,r_unbiased,'Linewidth', 4, 'DisplayName','Matlab XCorr UnBiased'); hold on; plot(lags./fs,r_DIY_unbiased', 'Linewidth', 2, 'DisplayName', 'DYI XCorr UnBiased'); legend;
xlabel('lags(s)'); ylabel('correlation'); title('XCorr UnBiased (Matlab vs DIY)');

% plot the matlab xcorr 'coef' and DIY 'biased' correlation on the same graph - they are the same! :)  
nexttile; plot(lags./fs,r_coef,'Linewidth', 4, 'DisplayName','Matlab XCorr Coef'); hold on; plot(lags./fs,r_DIY_coef', 'Linewidth', 2, 'DisplayName', 'DYI XCorr Coef'); legend;
xlabel('lags(s)'); ylabel('correlation'); title('XCorr Coeff (Matlab vs DIY)');