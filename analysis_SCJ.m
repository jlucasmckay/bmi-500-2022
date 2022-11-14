% reset the workspace
clear
close all

% load spiral drawing data
d = read_trc("lue-spiral.trc");

% set plotting parameters
TL = [0 5];
nr = 2;
nc = 3;

% plot the left hand marker in x-y-z
marker_name = "L.Finger3.M3";
marker_xyz = d{:,find(names(d) == "L.Finger3.M3") + [0:2]};

t = d{:,"Time"};
t_inds = t>min(TL)&t<max(TL);
t_secs = rem(t(t_inds),1)==0;


% plot
figure
subplot(nr,nc,1)
hold on

% TOP LEFT PLOT
plot(t, marker_xyz(:,1), 'LineWidth',4);
plot(t, marker_xyz(:,2), 'LineWidth',4);
plot(t, marker_xyz(:,3), 'LineWidth',4);
xlabel('seconds');
ylabel('mm');
title('Raw Data');
legend('X','Y','Z');
axis([0 5 0 1400])

hold off

% TOP MIDDLE PLOT
subplot(nr,nc,2)
hold on
% Uses only a sixth of the total data to make the plots
plot(marker_xyz(1:600,2),marker_xyz(1:600,3), 'k.', 'MarkerSize', 9);
xlabel('Y')
ylabel('Z')
title('Front View')
hold off

% Filter out large, slow movements with a high-pass butterworth filter at 2
% Hz cutoff and filter out jitter with a low-pass butterworth filter at 20
% Hz cutoff. A 6th order filter is fine.

% sampling freq fs is the reciprocal of the difference between two points
fs = 1/mean(diff(t));

% cutoff frequencies for the filter
fc_hi = 2;
fc_lo = 20;

% You say to use 20Hz for low pass butterworth, but the plot you had uses
% 2 Hz... I am confuse
% Low-pass butterworth filter, 6th-order, 2Hz cutoff
[b_lo, a_lo] = butter(6, 2/(fs/2), "low");

% Use the transfer function coefficients to transform the original data
lowpass_filt = filtfilt(b_lo, a_lo, marker_xyz);

% TOP RIGHT PLOT
subplot(nr,nc,3)
hold on
% Uses only a sixth of the total data to make the plots
plot(lowpass_filt(1:600,2),lowpass_filt(1:600,3), 'k.', 'MarkerSize', 9);
xlabel('Y')
ylabel('Z')
title('Low Frequency Component')
hold off

% [b,a] = butter(n,Wn) returns the transfer function coefficients of an 
% nth-order lowpass digital Butterworth filter with normalized 
% cutoff frequency Wn [https://www.mathworks.com/help/signal/ref/butter.html]

% Band-pass butterworth filter, 6th-order, 20Hz to 2Hz cutoff
[b_band, a_band] = butter(6, [fc_hi/(fs/2) fc_lo/(fs/2)], "bandpass");

% Use the transfer function coefficients to transform the original data
bandpass_filt = filtfilt(b_band, a_band, marker_xyz);

% calculate principal components oh bandpass filter
[pca_coeff, pca_scores, pca_latent] = pca(bandpass_filt);

% Get the first PC and the score's projection
first_princomp = pca_scores(:, 1) * pca_coeff(:, 1)';

% plot projections and line
subplot(nr,nc,4)
hold on
scatter(bandpass_filt(:, 2), bandpass_filt(:, 3), 16, 'filled', 'k')
y_values = -75:75;
z_values = y_values * bandpass_filt(3);
plot(y_values, z_values, 'LineWidth', 3);
axis([-75 75 -75 75])
hold off

% calculate projection onto first PC
proj = bandpass_filt * pca_coeff(:,1);

% smooth with a savitsky-golay smoother
proj_smooth = smoothdata(proj,'sgolay');

% count zero crossings
zcd = dsp.ZeroCrossingDetector();
numZeroCross = cast(zcd(proj_smooth(t_inds)),"double");
tremorFrequency = (numZeroCross/2)/max(TL);

% get envelope from 25 sample moving average
env_width = 25;
env = movmax(proj_smooth(t_inds),env_width);

% use the median of the moving maximum as the estimator of the amplitude
amp = median(env);

ttl = round(tremorFrequency,1) + " Hz, " + round(2*amp,1) + " mm amplitude";

% plot
subplot(nr,nc,[5 6])
hold on
plot(t,proj,'k.')
plot(t,proj_smooth,'r')
h1 = refline(0,amp);
h2 = refline(0,-amp);
h1.Color = 0.5*[1 1 1];
h2.Color = 0.5*[1 1 1];
xlim(TL)
title(ttl)
ylabel("mm")
xlabel("seconds")
