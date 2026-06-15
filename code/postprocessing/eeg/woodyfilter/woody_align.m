function [erpAvg, dg] = woody_align(data, times, alignChans, alignWin, maxShift)
% WOODY_ALIGN  Single-iteration Adaptive Woody filter (Gavin et al., 2019).
% Cluster-driven single shift: one lag per trial is estimated from the
% frontocentral cluster and applied to ALL channels of that trial, so
% topography is preserved. A single iteration is used (Gavin et al. show
% this is sufficient); the template is the condition's own trial average.
%
% INPUTS
%   data       [nChan x nTime x nTrials]  single trials for one condition
%   times      [1 x nTime]                EEG.times, in ms
%   alignChans  vector                    channels driving the shift, e.g. [1 2 33 34]
%   alignWin   [start end]                template-matching window in ms, e.g. [0 300]
%   maxShift    scalar                     max +/- shift in SAMPLES, e.g. 300
%
% OUTPUTS
%   erpAvg     [nChan x nTime]            latency-adjusted average
%   dg         struct                     diagnostics (Gavin's three measures):
%                 .meanR_pre  mean pre-shift correlation across trials
%                 .meanR_post mean max-correlation across trials
%                 .sdShift_ms SD of the shifts across trials, in ms
%                 .shifts/.rPre/.rPost  per-trial values
%
% Requires MATLAB R2016b+ (implicit expansion).

[nChan, nTime, nTrials] = size(data);

% ---- cluster waveform per trial: [nTime x nTrials] ----
clust = squeeze(mean(data(alignChans, :, :), 1));
if size(clust,1) == 1, clust = clust(:); end   % guard single-trial squeeze

% ---- template = average cluster waveform across trials ----
template = mean(clust, 2);

% ---- alignment-window sample indices and shift grid ----
winIdx = find(times >= alignWin(1) & times <= alignWin(2));
winIdx = winIdx(:);                  % [nWin x 1]
shifts = -maxShift:maxShift;         % [1 x nShift]

if min(winIdx)+min(shifts) < 1 || max(winIdx)+max(shifts) > nTime
    error('woody_align:bounds', ...
        'Alignment window +/- maxShift runs off the epoch edge.');
end

% template window, centered once (Pearson r is offset/scale invariant)
tw   = template(winIdx);
tw0  = tw - mean(tw);
twSS = sum(tw0.^2);

idxMat  = winIdx + shifts;           % [nWin x nShift] of sample indices
zeroCol = find(shifts == 0);         % column giving the unshifted (pre) r

bestShift = zeros(nTrials,1);
rPre      = zeros(nTrials,1);
rPost     = zeros(nTrials,1);

% ---- estimate best lag per trial (vectorized across all shifts) ----
for t = 1:nTrials
    seg = clust(:, t);
    S   = seg(idxMat);                                  % [nWin x nShift]
    S0  = S - mean(S, 1);
    r   = (tw0' * S0) ./ sqrt(twSS .* sum(S0.^2, 1));   % [1 x nShift]
    rPre(t)       = r(zeroCol);
    [rPost(t), k] = max(r);
    bestShift(t)  = shifts(k);
end

% ---- realign every channel of each trial, no wrap, then average ----
acc = zeros(nChan, nTime);
cnt = zeros(1, nTime);
for t = 1:nTrials
    s     = bestShift(t);
    src   = (1:nTime) + s;            % realigned(:,tt) = trial(:,tt+s)
    valid = src >= 1 & src <= nTime;
    acc(:, valid) = acc(:, valid) + data(:, src(valid), t);
    cnt(valid)    = cnt(valid) + 1;
end
erpAvg = acc ./ cnt;                  % [nChan x nTime]; edges only lose trials

% ---- diagnostics ----
srate         = 1000 / median(diff(times));     % ms spacing -> Hz
dg.meanR_pre  = mean(rPre);
dg.meanR_post = mean(rPost);
dg.sdShift_ms = std(bestShift) * (1000/srate);
dg.shifts     = bestShift;
dg.rPre       = rPre;
dg.rPost      = rPost;
end
