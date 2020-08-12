function [per, asdf] = funfreq(x,fs)
%FUNFREQ Estimation of the cycle length of a given periodic function x
%   This function estimates the period of a function by calculating the
%   average squared difference function and computing the distance of its local
%   maxima [1]
%
%   FUNFREQ(x, fs) calculates the period of signal x based on the sampling
%   frequency fs
%
% Literature:
% [1] McLeod, P., & Wyvill, G. (2005). A Smarter Way to Find Pitch. In ICMC
%     Proceedings (S. 138–141).
%
% Author: Fabian Schrumpf, MSc. 
% Laboratory for Biosignal Processing; HTWK Leipzig (Leipzig University of
% Applied Sciences)
% email address: fabian.schrumpf@htwk-leipzig.de  
% Website: https://labp.github.io/
% November 2016; Last revision: April 2017
%
% Edited by Christoph Moench, April 2018: Optimized Code and added the 
% "peak" at (1) and the last peak at length(asdf) before determining the
% frequency. These peaks were not taken into account before.

% compute average squared difference function (asdf)
N = length(x);
if all(x == 0)
    per = 0;
    asdf = zeros(1,length(x));
    return
end;
asdf = (2*ifft(fft(x).*conj(fft(x))) + sum(x.^2) + sum(x.^2))/N;
% asdf = -2*ifft(fft(x).*conj(fft(x)))/(sum(x.^2) + sum(x.^2));
% asdf = asdf - mean(asdf);

% calculate baseline as a 2nd order polynomial and remove it from asdf
f = fit(linspace(1,length(asdf),length(asdf))',asdf,'poly2','normalize','on');
bs = feval(f,linspace(1,length(asdf),length(asdf)));
asdf = asdf - bs;
asdf = asdf - mean(asdf);

% compute zero crossings of asdf
[idx_up, idx_down] = zerocrossing(asdf);

% detect local maxima between zero crossings
locs = zeros(length(idx_up) - 1, 1);
pks = zeros(length(idx_up) - 1, 1);

% we need at least 2 zerocrossings:
if (length(idx_up) < 2 || length(idx_down) < 2)
    per = 0;
    return;
end;

% find peaks between idx_up(j) and idx_down(j+1) (first zerocrossing is alway a "down crossing" in akf)
for j=1:length(idx_up) - 1
    try
        [pks(j),max_idx] = max(asdf(idx_up(j):idx_down(j + 1)));
        locs(j) = idx_up(j) + max_idx;
    catch
    end
end

% keep only peaks greater than 0.2*highestPeakValue:
pks_norm = pks./max(pks);
loc_norm = locs(pks_norm > 0.2);
% add the first "peak" and the one which would occur, if we continue the akf:
loc_norm = [1; loc_norm; length(asdf)];
per = median(diff(loc_norm))/fs;

if nargout == 2
    asdf = asdf';
end

end

function [idx_up, idx_down] = zerocrossing(v)
    sgn(v>=0) = 1;
    sgn(v<0) = -1;
    idx_up = find(diff(sgn) == 2)+1;
    idx_down = find(diff(sgn) == -2)+1;
end