function out = single_session_analysis(ex, sesinfo)
%%
% single session analysis in a learning project
%

% input
if nargin < 2; sesinfo = 1; end

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

% extract relevant info =========================
[~, ~, eyep] = getEyedata(ex);
behmat = getBeh(ex);
out = behmat(:,4)==1;
behmat(out, :) = [];
stmmat = getStmseq(ex, 1);
stmmat(out, :) = [];
v = stmmat(:);
disval = unique(v(~isnan(v)));
ntr = size(eyep, 1);

% stimlus types =========================
out = run_fit_routine(behmat, eyep, stmmat, disval);

% store trial matrix =============================
out.mat = [behmat, zeros(ntr, 4)];
ncol = size(behmat, 2);
params = define_params;
for n = 1:ntr
    % pupil size
    [minima, mini] = min(eyep(n, :));
    maxima = max(eyep(n, mini+1:end));
    out.mat(n, ncol+1) = maxima - minima;
    
    % pupil derivative
    dps = diff(eyep(n, :));
    out.mat(n, ncol+2) = max(dps);
    
    % pupil spectrum
    dps = dps - nanmean(dps);
    if n==1
        [out.pupilpow.S, out.pupilpow.f] = mtspectrumc(dps',  params); 
        frange = out.pupilpow.f >= 3 & out.pupilpow.f <= 10;
        [out.mat(n, ncol+3), idx] = max(out.pupilpow.S(frange));
    else
        s = mtspectrumc(dps',  params);
        [out.mat(n, ncol+3), idx] = max(s(frange));
        out.pupilpow.S = out.pupilpow.S + s;
    end    
    ff = out.pupilpow.f(frange);
    out.mat(n, ncol+4) = ff(idx);
end
out.pupilpow.S = out.pupilpow.S/ntr;

% available reward size ==========================
avrew = behmat(:, 14);
uniavrew = unique(avrew);
lena = length(uniavrew);
for a = 1:lena
   out.avrew{a} = run_fit_routine(behmat(avrew==a-1, :), ...
       eyep(avrew==a-1, :), stmmat(avrew==a-1, :), disval);
end

% serial choice bias =============================
out.fitse = serialchoicebias(behmat, [1 0.01]);

% basic session info =============================
if sesinfo == 1
    if isfield(ex,'fileID')
        ID = ex.fileID(1:3);
        filename = ex.fileName;
        if isfield(ex, 'dirName')
            dirname = ex.dirName;
        else
             dirname = [];
        end
    elseif isfield(ex,'Header')
        if isfield(ex.Header,'onlineFileName')
              ID = ex.Header.fileID(1:3);
            filename = ex.Header.onlineFileName;
            dirname = ex.Header.onlineDirName;
        elseif isfield(ex.Header,'Headers')
            ID =  ex.Header.Headers(1).fileID(1:3);
            filename = ex.Header.fileName;
            dirname = ex.Header.Headers(1).onlineDirName;
        else
         ID = ex.Header.fileID(1:3);
         filename = ex.Header.fileName;
         dirname = ex.Header.dirName;
        end
    else
        year = strfind(ex.fileName,'10000');
        dot = strfind(ex.fileName,'.');
        ID = strcat(ex.fileName(year:year+3),' ',ex.fileName(dot(1)+1:dot(2)-1),' ',ex.fileName(dot(2)+1:dot(2)+2));
        filename = ex.fileName;
        if isfield(ex, 'dirname')
                dirname = ex.dirName;
        else
                 dirname = [];
        end
    end
    out.ID = ID; out.filename = filename; out.dirname = dirname;
    out.refreshRate = ex.setup.refreshRate;
    out.stmdur = round(100*size(stmmat, 2)/ex.setup.refreshRate)/100;
    out.ntr = ntr; out.hdx_range = disval'; out.exp = {ex.exp.e1, ex.exp.e2};
    out.stmmat = stmmat;
end

function out = run_fit_routine(behmat, eyep, stmseq, disval)
% analysis as a function of stimulus type
%
stm = behmat(:, 5).*sign(behmat(:, 6));
x = unique(stm)';
lenu = length(x);
y = x;
n = x;
out.pupils = cell(2, 2); % mean and SEM
for i = 1:lenu
    % behavior
    cond = stm == x(i);
    n(i) = sum(cond);
    y(i) = sum(behmat(cond, 11)==1)/n(i);
    
    % pupil size
    ps = eyep(cond, :);
    out.pupils{1,1}(i, :) = nanmean(ps, 1);
    out.pupils{1,2}(i, :) = nanstd(ps, [], 1)/sqrt(n(i));
    
    % pupil derivative
    dps = [zeros(n(i), 1), diff(ps, 1, 2)];
    out.pupils{2,1}(i, :) = nanmean(dps, 1);
    out.pupils{2,2}(i, :) = nanstd(dps, [], 1)/sqrt(n(i));
end
% fit psychometric function
out.fitpm = fitPM(x, y, n, 'Weibull', 'MLE', 0);

% psychophysical kernel
nonsigtr = abs(behmat(:, 5)) < 0.005;
out.pk = getKernel(stmseq(nonsigtr, :), disval, behmat(nonsigtr, 11));

function pk = getKernel(stm, disval, ch)
% trial averaged stimulus distributions split by choice
nd = length(disval);
ntr = length(ch);
svmat = nan(ntr, nd);
for r = 1:ntr
    for d = 1:nd
        svmat(r,d) = sum(stm(r,:)==disval(d));
    end
end
% compute PK for 0% stimulus
pk = nanmean(svmat(ch==1,:)) - nanmean(svmat(ch==0,:));

function params = define_params
params.tapers = [2, 3]; % was [2, 3]. Chalk et al.,2010 used [3,5]
params.err = 0;
params.Fs = 500;
params.fpass = [0.01 100];
params.pad = 0; % nearest 2^N
params.trialave = 1; % if 0, the file is too heavy