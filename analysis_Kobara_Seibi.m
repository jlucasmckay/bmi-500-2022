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
x = marker_xyz(:,1);
y = marker_xyz(:,2);
z = marker_xyz(:,3);

t = d{:,"Time"};
t_inds = t>min(TL)&t<max(TL);
t_secs = rem(t(t_inds),1)==0;

% plot
figure
subplot(nr,nc,1)
hold on
plot(t,x)
plot(t,y)
plot(t,z)
legend("x","y","z")
xlim(TL)
title("Raw data")
xlabel("seconds")
ylabel("mm")

% plot Y-Z plane view
subplot(nr, nc, 2)
plot(y,z)
title("Front view")
xlabel("Y")
ylabel("Z")

% sampling frequency
fs = 1/mean(diff(t));

% cutoff frequency
fc_hi = 2;
fc_lo = 20;

% bandpass butterworth filter
[B, A] = butter(6, [fc_hi, fc_lo]/(fs/2), "bandpass");
x_band = filter(B, A, x);
y_band = filter(B, A, y);
z_band = filter(B, A, z);

% low pass filter 
[b, a] = butter(6, fc_hi/(fs/2));
x_low = filter(b, a, x);
y_low = filter(b, a, y);
z_low = filter(b, a, z);

% plot low pass filter
subplot(nr, nc, 3)
plot(y_low, z_low)
title("low frequency component")
xlabel("Y")
ylabel("Z")


% calculate PCA
xyz_table = table(x_band, y_band, z_band);
xyz_array = table2array(xyz_table);
[coeff, score, latent] = pca(xyz_array);

% PC1
pca1 = score(:,1);

% projection to pc1
proj = xyz_array*coeff(:,1)*coeff(:,1)';
% X*u*u

% plot high component and projection to PCA1
subplot(nr,nc, 4)
hold on
plot(xyz_array(:,2), xyz_array(:,3), 'b')
plot(proj(:,2), proj(:,3), 'r')
title("High frequency component and 1st PC")
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

