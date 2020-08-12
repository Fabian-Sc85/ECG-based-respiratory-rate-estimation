% Author: Fabian Schrumpf, MSc.
% Laboratory for Biosignal Processing; HTWK Leipzig (Leipzig University of
% Applied Sciences)
% email address: fabian.schrumpf@htwk-leipzig.de
% Website: https://labp.github.io/
% August 2020; Last revision: --

load('example_data/example_data.mat')

% define the ecg features to use for RR estimation
ecg_feature = 'R_amplitude';

% suppress powerline interference using a 2nd order 50 Hz IIR notch filter
wo = 50/(fs/2);  
bw = wo/35;
[b,a] = iirnotch(wo,bw);

ecg = filtfilt(b,a,ecg);

% loop over all time windows and estimate RR
Nwin = length(idx_start);
resp_asdf = zeros(Nwin, 1);
resp_fft = zeros(Nwin, 1);
resp_zc = zeros(Nwin, 1);

for i=1:Nwin
    % determine ECG for the current time window
    ecg_win = ecg(idx_start(i):idx_stop(i));
    
    % determine the QRS fiducial points falling into the current time
    % window
    Q_win = Q(and(Q>idx_start(i),S<idx_stop(i)));
    Q_win = Q_win - idx_start(i) + 1;
    R_win = R(and(Q>idx_start(i),S<idx_stop(i)));
    R_win = R_win - idx_start(i) + 1;
    S_win = S(and(Q>idx_start(i),S<idx_stop(i)));
    S_win = S_win - idx_start(i) + 1;
    
    % estimate the RR
    resp_asdf(i) = 60./ECG_calcRespRate(ecg_win, fs, [Q_win; R_win; S_win], ecg_feature, 'auto_corr');
    resp_fft(i) = 60./ECG_calcRespRate(ecg_win, fs, [Q_win; R_win; S_win], ecg_feature, 'fft');
    resp_zc(i) = 60./ECG_calcRespRate(ecg_win, fs, [Q_win; R_win; S_win], ecg_feature, 'zerocrossing');
    
end

% calculate absolute error between the estimated RR and the true RR
err_asdf = abs(resp_true - resp_asdf);
err_fft = abs(resp_true - resp_fft);
err_zc = abs(resp_true - resp_zc);

% plot the time courses for true and estimated RR
figure(1)
plot(resp_true)
hold on
plot(resp_asdf)
plot(resp_fft)
plot(resp_zc)
legend("ground truth","asdf","fft","zero crossing")
xlabel("time window #");
ylabel("RR / bpm");
grid;

% plot the absolute error for each estimation method
figure(2)
box_labels = {'asdf', 'fft', 'zero crossing'};
boxplot([err_asdf, err_fft, err_zc], 'labels', box_labels)
ylabel("|RR_{true} - RR_{est}| / bpm")
grid;