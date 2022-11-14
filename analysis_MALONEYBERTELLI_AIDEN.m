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

%%% YOUR CODE HERE
marker_inds = [find(t == 1), find(t == 2), find(t == 3), find(t == 4)];
plot(t(t_inds), marker_xyz(t_inds, :), '.-')
plot(t(t_inds), marker_xyz(t_inds, :), 'k|', 'MarkerIndices', marker_inds, ...
    'MarkerSize', 10)
xlim([0,max(t(t_inds))])
ylim([0 1400])
title("Raw Data")
xlabel("seconds")
ylabel("mm")
legend("X", "Y", "Z")

subplot(nr,nc,2)
hold on
plot(marker_xyz(t_inds,2), marker_xyz(t_inds,3), 'k.-')
plot(marker_xyz(t_inds,2), marker_xyz(t_inds,3), 'k_', ...
    'MarkerIndices', marker_inds, 'MarkerSize', 10)
title("Front View")
xlabel("Y")
ylabel("Z")

% Filter out large, slow movements with a high-pass butterworth filter at 2
% Hz cutoff and filter out jitter with a low-pass butterworth filter at 20
% Hz cutoff. A 6th order filter is fine.

% sampling freq fs is the reciprocal of the difference between two points
fs = 1/mean(diff(t));

% cutoff frequencies for the filter
fc_hi = 2;
fc_lo = 20;

% [b,a] = butter(n,Wn) returns the transfer function coefficients of an 
% nth-order lowpass digital Butterworth filter with normalized 
% cutoff frequency Wn [https://www.mathworks.com/help/signal/ref/butter.html]

%%% YOUR CODE HERE
[b_fig3, a_fig3] = butter(3, fc_hi/(fs/2));
marker_xyz_fig3 = filtfilt(b_fig3, a_fig3, marker_xyz);

subplot(nr,nc,3)
hold on
plot(marker_xyz_fig3(t_inds,2), marker_xyz_fig3(t_inds,3), 'k.-')
plot(marker_xyz_fig3(t_inds,2), marker_xyz_fig3(t_inds,3), 'k_', ...
    'MarkerIndices', marker_inds, 'MarkerSize', 10)
title("Low frequency component")
xlabel("Y")
ylabel("Z")

[b_hi, a_hi] = butter(3, fc_hi/(fs/2), 'high');
marker_xyz_hi = filtfilt(b_hi, a_hi, marker_xyz);
[b_lo, a_lo] = butter(3, fc_lo/(fs/2));
marker_xyz_hi_lo = filtfilt(b_lo, a_lo, marker_xyz_hi);

% calculate the first PC

%%% YOUR CODE HERE
coeff = pca(marker_xyz_hi_lo);
pc1 = coeff(:,1);

% calculate projection onto first PC

%%% YOUR CODE HERE
proj = marker_xyz_hi_lo * pc1;
pc1_y = -25:1:25;
pc1_z = pc1(3)/pc1(2) * pc1_y;

subplot(nr,nc,4)
hold on
plot(marker_xyz_hi_lo(t_inds,2), marker_xyz_hi_lo(t_inds,3), 'k.-', ...
    'LineWidth', 2, 'MarkerSize', 10)
plot(pc1_y, pc1_z, 'r')
title("High frequency component and 1st PC")
xlim([-75 75])
ylim([-90 90])
xlabel("Y")
ylabel("Z")


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

ylim([-20 20])
