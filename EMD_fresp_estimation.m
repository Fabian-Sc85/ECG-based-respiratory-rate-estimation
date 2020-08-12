function [f_resp, resp_signal, imf] = EMD_fresp_estimation(signal, fs, varargin)
%EMD_fresp_estimation Compute the respiratory frequency from a respiratory
%signal based on the Empirical Mode Decomposition (EMD)
%   F_RESP = EMD_fresp_estimation(SIGNAL, FS) returns the respiratory
%   frequency F_RESP based on the data vector SIGNAL and sampling frequency
%   FS based on the EMD Algorithm [1]
%
%   [... , RESP_SIGNAL] = EMD_fresp_estimation(SIGNAL, FS) returns the IMF,
%   that has been identified as the respiratory signal in the data vector
%   RESP_SIGNAL
%
%   [...] = EMD_fresp_estimation(..., 'f_range', FREQ_RANGE) uses the upper
%   and lower frequency boundary in FREQ_RANGE (FREQ_RANGE(1): lower
%   boundary; FREQ_RANGE(2): uper boundary) during EMD of the respiratory
%   signal
%
% Literature:
% [1] Madhav, K. V., Ram, M. R., Krishna, E. H., Komalla, N. R., & Reddy, 
%     K. A. (2011). Estimation of respiration rate from ECG, BP and PPG 
%     signals using empirical mode decomposition. In Conference Record - IEEE
%     Instrumentation and Measurement Technology Conference (pp. 1661–1664). IEEE.
%
% Author: Fabian Schrumpf, MSc.
% Laboratory for Biosignal Processing; HTWK Leipzig (Leipzig University of
% Applied Sciences)
% email address: fabian.schrumpf@htwk-leipzig.de
% Website: https://labp.github.io/
% March 2018; Last revision:

p = inputParser;
addOptional(p,'f_range',[0.1 0.6],@isnumeric);
parse(p,varargin{:});
f_range = p.Results.f_range;

if f_range(1) == 0
    f_range(1) = 0.01;
end

imf = emd(signal,'MAXMODES',40,'MAXITERATIONS',1000);

bp_peak = zeros(1,size(imf,1));
bp_peakamp = zeros(1,size(imf,1));

for i=1:size(imf,1)
    n_fft = 2^floor(log2(60*fs));                                   % increase spectral resolution by zero padding
    win = hann(size(imf,2))';                                       % use HANN window for FFT
    if n_fft > size(imf,2)
        Y = fft(imf(i,:).*win,n_fft);
    else
        Y = fft(imf(i,:));
    end
    f = linspace(0,fs,length(Y));
    idx_peak = find(abs(Y) == max(abs(Y)),1);                       % find frequency peaks in each IMF-spectrum
    bp_peak(i) = f(idx_peak);
    bp_peakamp(i) = abs(Y(idx_peak));
end

bp_peakamp(~and(bp_peak>=f_range(1),bp_peak<=f_range(2))) = 0;
[~, idxmax] = max(bp_peakamp);
resp_signal = imf(idxmax,:)';

f_resp = funfreq(resp_signal,fs);
end

