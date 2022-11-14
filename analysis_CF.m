% Chase Fensore: BMI 500, Nov 14, 2022
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

% PLOT 1: ----------------------------------
% Plot the X, Y, Z data of the L.Finger3.M3 marker.
% x-axis: t_secs
% y-axis: X, Y, Z data in marker_xyz
plot(t(t_inds), marker_xyz(t_inds,:), '- .');
xlabel("seconds");
ylabel("mm");
xlim([0 5]);
legend("X", "Y", "Z");
title("Raw Data");


% PLOT 2: ----------------------------------
% 2. Plot the Y-Z front view.
% plot
subplot(nr,nc,2)
hold on

plot(marker_xyz(t_inds,2), marker_xyz(t_inds, 3), '- . black'); % Z vs Y
xlabel("Y");
ylabel("Z");
title("Front View");


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

%%% YOUR CODE HERE:
% but... the above note from Prof says the opposite...

[b, a] = butter(6, fc_hi/(fs/2),'high'); % 6th-order hi butter, cut=2 Hz
[d, c] = butter(6, fc_lo/(fs/2)); % 6th-order lo butter, cut=20 Hz
[f, e] = butter(6, 2/(fs/2)); % low-pass at 2 Hz

% out_lo is used for later analysis data.
% Note that filter produced results not matching Fig 3, so filtfilt
% was used instead.
% For future analaysis, we use all 3600 examples.
out_hi = filtfilt(b, a, marker_xyz); 
out_lo = filtfilt(d, c, out_hi); % *** For use in later analysis

% For plotting in Fig 3., use only 600 examples.
plot_lo_pass = filtfilt(f, e, marker_xyz); 


% PLOT 3 ---------------------------------------
subplot(nr,nc,3)
hold on
% Plot Low frequency component:
plot(plot_lo_pass(1:600,2), plot_lo_pass(1:600,3), 'k', 'LineWidth', 2);
xlabel("Y");
ylabel("Z");
title("Low frequency component");


% calculate the first PC

%%% YOUR CODE HERE

coeff = pca(out_lo); 
% Now, columns in coeff correspond to PCs, each row is all for one input var.

% calculate projection onto first PC

%%% YOUR CODE HERE
% To project out_lo onto 1st PC, we multiply out_lo * coeff[:,1]
% PCA Reference: https://www.mathworks.com/matlabcentral/answers/259957-how-to-apply-pca-correctly

proj = out_lo * coeff; % proj is (3x3)


subplot(nr,nc,4)
hold on
% Plot the first PC
% Plot orig filtered data here.
plot(out_lo(t_inds,2), out_lo(t_inds,3), 'k.', 'MarkerSize', 12);

% Next, plot the PC1 line:
yDomain = -60:0.25:60; % Refinement of 0.25
m = coeff(3,1) / coeff(2,1); % m_PC1 = zCoeff/yCoeff
% Begin at origin (0,0)
y1 = 0;
z1 = 0;
zPCvalues = m*(yDomain-y1) + z1; % Calculate z values for PC1
plot(yDomain,zPCvalues, 'LineWidth', 1); % Plot PC1 line.

xlabel("Y");
ylabel("Z");
xlim([-75 75]);
ylim([-75 75]);
title("High frequency component and 1st PC");
hold off

proj = proj(:,1); % Trim proj to only include X coordinates.

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
