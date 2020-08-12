function [resp_rate_est, resp_signal_est, auto_corr, resp_depth] = ECG_calcRespRate(ecg,fs,varargin)
%ECG_calcRespRate Derivation of the respiratory rate from a given ecg
%signal
%   This function extracts the respiratory signal from an ecg signal by
%   exploiting features like R-peak amplitude, heart-rate variability and
%   cardiac axis. The respiratory rate is computed from the resulting
%   respiratory signal by either autocorrelation or analysis in frequency
%   domain
%
% Syntax: resp_rate_est = ECG_calcRespRate(ecg,fs,method, criterion)
%
% Inputs:
%   ecg         -       preprocessed ecg signal
%   fs          -       sampling frequency of the ecg signal
%   method      -       Method according to which the respiratory rate should be
%                       calculated. Possible Values:
%                      * QRS_area: area enclosed under the QRS-Complex is
%                        evaluated [1,2]
%                      * R_amplitude: Amplitude of the R-peak is evaluated [1, 3]
%                      * lead_angle: angle of the cardiac axis is evaluated (only
%                        for multi-channel electrocardiograms) [1, 2]
%                      * HRV: Modulation of the length of the RR-intervall (time
%                        between two R-peaks == heart-rate-variability) is
%                        evaluated
%                       default: R_amplitude
%   qrs_complexes   -   precomputed QRS complexes in the N_qrs-by-3 Matrix QRS. 
%                       If they are not provided, QRS-complexes are calculated
%                       as needed
%                      * QRS(:,1): Q-Points
%                      * QRS(:,2): R-Points
%                      * QRS(:,3): S-Points
%   criterion       -   Method according to which the resulting respiration signal
%                       is beeing analyzed. Possible values are:
%                       * auto_corr: Calculation of the base frequency in time
%                         domain using autocorrelation [5,6]
%                       * fft: detection of the dominant frequency in the frequency
%                         domain using fft
%                       default: auto_corr
%   f_range         -   frequency range that is used for EMD-based
%                       repiratory signal estimation
%
% Outputs:
%   resp_rate_est   -   estimated respiratory cylce in seconds
%
% Literature:
% [1] Sohrt-petersen, L. (2013). Evaluation of Algorithms for ECG Derived
%     Respiration in the Context of Heart Rate Variability Studies. Aalborg University.
% [2] Moody, G. B., & Mark, R. G. (1985). Derivation of respiratory signals
%     from multi-lead ECGs. Computers in Cardiology, 12, 113–116. doi:10.1109/FBIE.2008.41
% [3] Mason, C. L., & Tarassenko, L. (2001). Quantitative assessment of respiratory
%     derivation algorithms. 2001 Conference Proceedings of the 23rd Annual
%     International Conference of the IEEE Engineering in Medicine and Biology
%     Society, 2, 1998–2001. doi:10.1109/IEMBS.2001.1020622
% [4] Dingab, S., Zhua, X., Chena, W., & Weia, D. (2004). Derivation of
%     respiratory signal from single-channel ECGs based on Source Statistics.
%     International Journal of Bioelectromagnetism, 6(1).
% [5] McLeod, P., & Wyvill, G. (2005). A smarter way to find pitch. In Proceedings
%     of International Computer Music Conference, ICMC.
% [6] De Cheveigné, A., & Kawahara, H. (2002). YIN, a fundamental frequency
%     estimator for speech and music. The Journal of the Acoustical Society
%     of America, 111(4), 1917 – 1930. doi:10.1121/1.1458024
%
% Author: Fabian Schrumpf, MSc.
% Laboratory for Biosignal Processing; HTWK Leipzig (Leipzig University of
% Applied Sciences)
% email address: fabian.schrumpf@htwk-leipzig.de
% Website: https://labp.github.io/
% July 2015; Last revision: 09.03.2018

expected_method = {'QRS_area','R_amplitude','lead_angle','HRV', 'HRV2'};
expected_criterion = {'auto_corr','fft','zerocrossing'};

p = inputParser;
addOptional(p,'method','R_amplitude',...
    @(x) any(validatestring(x,expected_method)));
addOptional(p,'criterion','auto_corr',...
    @(x) any(validatestring(x,expected_criterion)));
addOptional(p,'qrs_complexes',[],@isnumeric);
addOptional(p,'f_range',[0.1, 0.6],@isnumeric);
parse(p,varargin{:});
% ecg = p.Results.ecg;
% fs = p.Results.fs;
method = p.Results.method;
criterion = p.Results.criterion;
QRS = p.Results.qrs_complexes;
f_range = p.Results.f_range;
t = linspace(0,length(ecg)/fs,length(ecg))';

switch method
    case 'QRS_area'
        T_meas = length(ecg)/fs;
        if T_meas < 3
            warning('Signal duration shorter than 3 s. Accuracy may be insufficient')
        end
        
        if isempty(QRS)
            [Q, R, S] = ECG_QRSdetector(ecg,fs);
        else
            Q = QRS(1,:);
            R = QRS(2,:);
            S = QRS(3,:);
        end
        peak_area = zeros(length(R),1);
        
        for i=1:length(Q)
            if S(i) - Q(i) > 3
                peak_area(i) = trapz(t(Q(i):S(i)),ecg(Q(i):S(i)));
                F = [1 t(Q(i)) ecg(Q(i));...
                    1 t(R(i)) ecg(R(i));...
                    1 t(S(i)) ecg(S(i))];
                peak_area(i) = 0.5 * abs(det(F));
            end
        end
        Q(peak_area==0) = [];
        R(peak_area==0) = [];
        S(peak_area==0) = [];
        peak_area(peak_area==0) = [];
        
        %         peak_area = peak_area - mean(peak_area);
        resp_signal_est = ppval(spline(t(R),peak_area),t(R(1):R(end)));
        resp_signal_est = resp_signal_est - mean(resp_signal_est);
        resp_signal_est = resp_signal_est./max(abs(resp_signal_est));
        %         resp_signal_est = padarray(resp_signal_est,[length(ecg) - length(resp_signal_est),0],'post');
    case 'R_amplitude'
        T_meas = length(ecg)/fs;
        if T_meas < 3
            warning('Signal duration shorter than 3 s. Accuracy may be insufficient')
        end
        
        if isempty(QRS)
            [~, R, S] = ECG_QRSdetector(ecg,fs);
        else
            R = QRS(2,:);
            S = QRS(3,:);
        end

        R_amp = ecg(R);
        S_amp = ecg(S);
        R_amp_interp = ppval(spline(t(R),abs(S_amp - R_amp)),t(R(1):R(end)));
        resp_signal_est = R_amp_interp - mean(R_amp_interp);
        resp_signal_est = resp_signal_est./max(abs(resp_signal_est));
%         [~, resp_signal_est] = EMD_fresp_estimation(resp_signal_est, fs, 'f_range', f_range);
        %         resp_signal_est = padarray(resp_signal_est,[length(ecg) - length(resp_signal_est),0],'post');
    
    case {'HRV2', 'HRV'}
        if isempty(QRS) || size(QRS,2) < 5
            [~, R] = ECG_QRSdetector(ecg,fs);
        else
            R = QRS(2,:);
        end
        
        % there must be at least 5 R-peaks to compute a relieable
        % respiratory rate
        if length(R) < 5
            resp_rate_est = NaN;
            auto_corr = NaN;
            resp_signal_est = NaN;
            resp_depth = NaN;
            return;
        end
        
        t = linspace(0,R(end-1)./fs,R(end-1));
%         smooth_RR = smoothdata(diff(R),'gaussian',1);
        resp_signal_est = ppval(spline(t(R(1:end-1)),diff(R)),t)';
        resp_signal_est = resp_signal_est - mean(resp_signal_est);
        resp_signal_est = resp_signal_est./max(abs(resp_signal_est));
        
        [~, resp_signal_est] = EMD_fresp_estimation(resp_signal_est, fs, 'f_range', f_range);
        
    otherwise
        error(['Method ' method ' not supported']);
end

resp_depth = iqr(resp_signal_est);
switch criterion
    case 'auto_corr'
        
        [resp_rate_est, auto_corr] = funfreq(resp_signal_est,fs);
            
    case 'fft'
        auto_corr = [];
        %             fft_length = 10 * fs;
        Y = abs(fft(resp_signal_est));
        f_res = fs/length(Y);
        f = linspace(0,fs,length(Y));
        Y = Y(f<3);
        [~,idx] = max(Y);
        resp_rate_est = 1/(idx*f_res);

    case 'zerocrossing'
        auto_corr = [];
        zc = zerocrossing(resp_signal_est);
        T_resp = 2*mean(diff(zc));
        resp_rate_est = T_resp/fs;
end
% resp_signal_est = padarray(resp_signal_est,[length(ecg) - length(resp_signal_est),0],'post');
end

