function [bi, ps, dps] = get_psconf_list(ex)
%%
% INPUT: ex file
% OUTPUT: bi ...
% binary list of trial-by-trial 'confidence' inferred from
% pupil size (column 1) and pupil derivative (column 2)
%                 ps ... trials x pupil size (only during the stimulus
%                 presentation period)
%                 dps ... pupil derivative, same dimension as ps
%

% white noise correction =========================
try
    ex = salvageFrames(ex);
catch
    disp('For this ex-file, "salvageFrames" caused an error')
end

% sort trials chronologically ====================
if isfield(ex.Trials, 'TrialStartDatapixx')
    [~,sorted_index] = sort([ex.Trials.TrialStartDatapixx]);
else
    [~,sorted_index] = sort([ex.Trials.TrialStart]);
end
ex.Trials = ex.Trials(sorted_index);

% extract pupil data =========================
[~, ~, eyep] = getEyedata(ex);
behmat = getBeh(ex);
out = behmat(:,4)==1;
behmat(out, :) = [];
ntr = size(eyep, 1);

% compute pupil metrics
mat = zeros(ntr, 2);
for n = 1:ntr
    % pupil size
    [minima, mini] = min(eyep(n, :));
    maxima = max(eyep(n, mini+1:end));
    mat(n, 1) = maxima - minima;
    
    % pupil derivative
    mat(n, 2) = max(diff(eyep(n, :)));
end

% median split
stm = behmat(:, 5).*sign(behmat(:, 6));
unique_stm = unique(stm);
lenu = length(unique_stm);
bi = zeros(ntr, 2);
for i = 1:lenu
    % trials with the same stimulus 
    cond = find(stm == unique_stm(i));
    
    % median split
    [~, idx2] = median_split(mat(cond, 1));
    bi(cond(idx2), 1) = 1;
    [~, idx2] = median_split(mat(cond, 2));
    bi(cond(idx2), 2) = 1;
end

% outputs
ps = eyep;
dps = [zeros(ntr, 1), diff(eyep, 1, 2)];

% subfunction
function [idx1, idx2] = median_split(v)
[~, idx] = sort(v, 'ascend');
nlen = length(idx);
idx1 = idx(1:floor(nlen/2));
idx2 = idx(floor(nlen/2)+1:end);