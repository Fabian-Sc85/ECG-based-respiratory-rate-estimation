function [Q_loc,R_loc,S_loc, varargout] = ECG_QRSdetector(ekg,fs,A)
%ECG_QRSdetector Detect the QRS-Complexes in a given ecg-Signal
%   This function returns the positions (sample index) of QRS-peaks in an
%   electroencephalogram based on a Hamilton-Tompkins QRS-Detector [1,2].
%
% Syntax: [Q_loc,R_loc,S_loc] = ECG_QRSdetector(ekg,fs)
%
% Inputs:
%   ekg     -   Electrocardiogram
%   fs      -   Sampling Frequency of ekg
%   A       -   Initial threshold for peak detection (default A = 0.8)  
%
% Outputs:
%   Q_loc   -   Positions (sample number) of the detected Q-Peaks in ekg
%   R_loc   -   Positions (sample number) of the detected R-Peaks in ekg
%   S_loc   -   Positions (sample number) of the detected R-Peaks in ekg
%
% Literature:
% [1] Hamilton, P., & Tompkins, W. (1986). Quantitative investigation of QRS
%     detection rules using the MIT/BIH arrhythmia database. IEEE Transactions
%     on Biomedical Engineering, 33(12), 1157–1165.
% [2] Sohrt-petersen, L. (2013). Evaluation of Algorithms for ECG Derived
%     Respiration in the Context of Heart Rate Variability Studies. Aalborg University.
%
% Author: Fabian Schrumpf, MSc. 
% Laboratory for Biosignal Processing; HTWK Leipzig (Leipzig University of
% Applied Sciences)
% email address: fabian.schrumpf@htwk-leipzig.de  
% Website: https://labp.github.io/
% September 2015; Last revision: 23.10-2015

if nargin < 3
    A = 0.8;
end;

if (max(abs(ekg)) > 1)
    ekg = ekg ./ max(abs(ekg));
end;

filt_low = designfilt('lowpassfir','PassbandFrequency',15,'StopbandFrequency',30,...
    'PassbandRipple',1,'StopbandAttenuation',80,'SampleRate',fs);
filt_high = designfilt('highpassfir','PassbandFrequency',5,'StopbandFrequency',1,...
    'PassbandRipple',1,'StopbandAttenuation',80,'SampleRate',fs);
delay = unique(grpdelay(filt_high.Coefficients)) + unique(grpdelay(filt_low.Coefficients));

win_len_init = 20*fs;
if win_len_init >= length(ekg)
    N_win = 1;
    idx_start = 1;
    idx_stop = length(ekg);
else
    N_win = ceil(length(ekg)/win_len_init);
    win_len = floor(length(ekg)/N_win);
    idx_start = floor(linspace(1,length(ekg)-win_len,N_win));
    idx_stop = idx_start + win_len - 1;
end;


ekg_avg = zeros(1,length(ekg));
% ROI = [];
for i=1:N_win
    ekg_filt = filtfilt(filt_high,filtfilt(filt_low,ekg(idx_start(i):idx_stop(i))));
    ekg_diff = conv(ekg_filt,[-2 -1 0 1 2],'same')*1/8;
    % ekg_diff = zeros(size(ekg,1),size(ekg,2));
    % for i=4:length(ekg_diff)
    %     ekg_diff(i) = 1/8 * (2*ekg_filt(i) - ekg_filt(i-1) - ekg_filt(i-2) - 2*ekg_filt(i-3));
    % end;
%     ekg_diff = conv(ekg_filt,[-2 -1 0 1 2],'same')*1/8;
    % delay = delay + 5;
    ekg_squared = ekg_diff.^2;
    
    ekg_avg_temp = conv(ekg_squared,ones(1,32),'same')./32;
%     ekg_avg = ekg_avg(round(delay):end);
    ekg_avg_temp = ekg_avg_temp./max(abs(ekg_avg_temp));
    ekg_avg(idx_start(i):idx_stop(i)) = ekg_avg_temp;
end;
ROI = ekg_avg>0.05;


ROI_edge = diff(ROI);
ROI_rising = find(ROI_edge == 1)+1;
ROI_falling = find(ROI_edge == -1);

if ROI_rising(1) > ROI_falling(1)
    ROI_falling(1) = [];
end;

n_QRS = min([length(ROI_rising),length(ROI_falling)]);
R = zeros(1,n_QRS);
R_loc = zeros(1,n_QRS);
Q = zeros(1,n_QRS);
Q_loc = zeros(1,n_QRS);
S = zeros(1,n_QRS);
S_loc = zeros(1,n_QRS);

thresh = A;
for i=1:n_QRS
    [pk_amp, pk_loc] = max(ekg(ROI_rising(i):ROI_falling(i)));
    if pk_amp >= thresh
        R(i) = pk_amp;
        R_loc(i) = pk_loc;
        if i < 7
            thresh = A*median(R(1:i));
        else
            thresh = A*median(R(i-6:i));
        end;
%         thresh = 0.8*R(i);
    else
        thresh = A*thresh;
        continue
    end;
%     [R(i),R_loc(i)] = max(ekg(ROI_rising(i):ROI_falling(i)));
    R_loc(i) = R_loc(i) + ROI_rising(i) - 1;
    [Q(i),Q_loc(i)] = min(ekg(ROI_rising(i):R_loc(i)));
    Q_loc(i) = Q_loc(i) + ROI_rising(i) - 1;
    [S(i), S_loc(i)] = min(ekg(R_loc(i):ROI_falling(i)));
    S_loc(i) = S_loc(i) + R_loc(i) - 1;
end;
Q_loc(Q_loc==0) = [];
R_loc(R_loc==0) = [];
S_loc(S_loc==0) = [];

if nargout > 3
    bpFilt = designfilt('bandpassfir','FilterOrder',80, ...
         'CutoffFrequency1',0.5,'CutoffFrequency2',40, ...
         'SampleRate',fs);
    ekg = filtfilt(bpFilt,ekg);
    t_onset = zeros(1,length(R_loc ));
    t_offset = zeros(1,length(R_loc));
    ecg_e = abs(hilbert(ekg));
    
    r = ceil(fs/250);
    d_ekg_e = zeros(1,length(ekg));
    for i=1+2*r:length(d_ekg_e)-2*r
        d_ekg_e(i) = (2*(ecg_e(i+2*r) - ecg_e(i-2*r)) + ecg_e(i+r) - ecg_e(i-r))/10;
    end;
    AS = 2*d_ekg_e.^2;
    
    for i=1:length(R_loc)
        idx_start = R_loc(i) - round(0.3*fs);
        if idx_start < 1
            idx_start = 1;
        end;
        idx_stop = R_loc(i);
        ASwin=AS(idx_start:idx_stop);
        lambda = max(AS(idx_start:(idx_start+round(0.05*fs))));
        mu0 = mean(AS(idx_start:(idx_start+round(0.05*fs))));
        mu1 = mean(AS(idx_start:(R_loc(i))));
        sigma = var(ASwin);
        g = zeros(1,length(ASwin));
        
        for j=2:length(g)
            p1 = 1/(sigma*sqrt(2*pi)) * exp(-(ASwin(j)-mu1)^2/(2*sigma));
            p0 = 1/(sigma*sqrt(2*pi)) * exp(-(ASwin(j)-mu0)^2/(2*sigma));
            g(j) = g(j-1) + log(p1/p0);
            if g(j) < 0
                g(j) = 0;
            end;
        end;
        idx = find(g>lambda,1,'first');
        t_onset(i) = idx_start + idx;
        
        idx_start = R_loc(i);
        idx_stop = R_loc(i) + round(0.15*fs);
        if idx_stop > length(AS)
            idx_stop = length(AS);
        end;
        ASwin=AS(idx_start:idx_stop);
        lambda = max(AS(R_loc(i)+round(0.1*fs):idx_stop));
        if R_loc(i)+round(0.1*fs) >= idx_stop
            continue;
        end;
        mu0 = mean(AS(R_loc(i)+round(0.1*fs):idx_stop));
        mu1 = mean(AS(R_loc(i):idx_stop));
        ASwin = fliplr(ASwin);
        sigma = var(ASwin);
        g = zeros(1,length(ASwin));
        for j=2:length(g)
            p1 = 1/(sigma*sqrt(2*pi)) * exp(-(ASwin(j)-mu1)^2/(2*sigma));
            p0 = 1/(sigma*sqrt(2*pi)) * exp(-(ASwin(j)-mu0)^2/(2*sigma));
            g(j) = g(j-1) + log(p1/p0);
            if g(j) < 0
                g(j) = 0;
            end;
        end;
        idx = length(g) - find(g>lambda,1,'first');
        t_offset(i) = idx_start + idx;
    end;
    
    varargout{1} = t_onset;
    varargout{2} = t_offset;
end;

% del_idx = [];
% 
% n_QRS = length(Q_loc);
% for i=1:n_QRS
%     if Q_loc(i) == R_loc(i) && R_loc(i) == S_loc(i) && S_loc(i) == Q_loc(i)
%         del_idx = [del_idx, i];
%     end;
% end;
% 
% Q_loc(del_idx) = [];
% R_loc(del_idx) = [];
% S_loc(del_idx) = [];

end

