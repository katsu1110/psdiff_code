function visualize(data, focus, saveoption)
% visualize results from 'head_free_getData'
% INPUT:
% data ... download data from the followings:
% load('Z:\Katsuhisa\headfree_project\dataset\fixdata.mat')
% load('Z:\Katsuhisa\headfree_project\dataset\disdata.mat')
%
% focus ... either 'params', 'behaviors', 'fixSpan'
%
% saveoption ... 1; save
%

if nargin < 2; focus = 'all'; end
if nargin < 3; saveoption = 0; end

% path
addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_usefulfunc'))
addpath(genpath('Z:\Katsuhisa\headfree_project\code'))
figpath = 'Z:\Katsuhisa\headfree_project\figures\';

%%
% the number of sessions
close all;
lenses = length(data.session);
out = zeros(1, lenses);
for i  = 1:lenses
    if ~isstruct(data.session(i).params) || ~isstruct(data.session(i).psarti)
        out(i) = 1;
        disp(['session ' num2str(i) ' is excluded as parameters are not stored'])
    end
end
data.session(out==1) = [];
% % remove ses 9, 10 (shorter duration) in discrimination
% if length(data.session) < 20
%     data.session(9:10) = [];
% end
% remove DCxDX stimuli
if ~isfield(data.session, 'pm')
    data.session(end-3:end) = [];
end
lenses = length(data.session);
disp([num2str(lenses) ' sessions are included...'])

% task type
if isfield(data.session, 'pm')
    tasktype = '_dis';
else
    tasktype = '_fix';
end

% blue and red
red = [0.7922 0 0.1255];
blue = [0.0196 0.4431 0.6902];

%%
% experimental parameters (fixation duration per trial
% fixation window, fixation dot size)
if strcmp(focus, 'all') || strcmp(focus, 'params')
    params = nan(3, lenses);
    for i = 1:lenses
        params(1,i) = nanmean(data.session(i).params.fixdur);
        params(2,i) = 4*nanmean(data.session(i).params.fixwin(1,:).*data.session(i).params.fixwin(2,:));
        params(3,i) = nanmean(data.session(i).params.fixdotsz);
    end

    paranames = {'fixation duration per trial', 'fixation window (deg^2)', 'fixation dot size'};
    session_plot(params, [], paranames, 'experimental parameters')
    
    % autosave figures
    if saveoption == 1
        savefig(gcf, [figpath 'Figure2_ExperimentalParams\raw_figs\expparams' tasktype '.fig'])
    end
end

%%
% animal's behaviors (working hours, P(fixation breaks), suvival function,
% its fitted parameter)
if strcmp(focus, 'all') || strcmp(focus, 'behaviors')
    params = nan(4, lenses);
    for i = 1:lenses
        params(1,i) = data.session(i).eyedata.eyeveclen/(500*60);
        params(2,i) = 100*sum(data.session(i).eyedata.reward==0)/...
            length(data.session(i).eyedata.reward);
        % median
        params(3,i) = data.session(i).survival.fitparams(1)*log(2)^(1/data.session(i).survival.fitparams(2));
%         params(3,i) = (nanmean(data.session(i).params.fixdur/2)/data.session(i).survival.fitparams(1))...
%             ^data.session(i).survival.fitparams(1);
        params(4,i) = (nanmean(data.session(i).params.fixdur)/data.session(i).survival.fitparams(1))...
            ^data.session(i).survival.fitparams(1);
    end

    % paranames = {'working hours (h)', 'proportion of trials with fixation breaks (%)', 'Weibull \k'};
    paranames = {'working duration (min)', 'proportion of trials with fixation breaks (%)', ...
        'median of weibull dist.', 'cumulative hazard rate'};
    session_plot(params, [], paranames, 'animal behaviors')
    
    % autosave figures
    if saveoption == 1
        savefig(gcf, [figpath 'Figure2_ExperimentalParams\raw_figs\ses_behav' tasktype '.fig'])
    end
    
    % survival functions with a Weibull fit
    figure;
    nr = floor(sqrt(lenses));
    nc = ceil(lenses/nr);
    for i = 1:lenses
        subplot(nr, nc, i)
        stairs(data.session(i).survival.edges, data.session(i).survival.survival, '-k')
        hold on;
        fitx = linspace(data.session(i).survival.edges(1), ...
            data.session(i).survival.edges(end), 100);
        plot(fitx, 1 - cdf('wbl', fitx, ...
            data.session(i).survival.fitparams(1), data.session(i).survival.fitparams(2)), '-r')
        xlim([0 max(data.session(i).params.fixdur)])
        ylim([0 1])
        title(num2str(i))
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    end
    set(gcf, 'Name', 'survival functions', 'NumberTitle', 'off')
    
    % autosave figures
    if saveoption==1
        savefig(gcf, [figpath 'Figure2_ExperimentalParams\raw_figs\survival_all' tasktype '.fig'])
    end
end

% %% 
% % pupil size artifact on gaze position
% if strcmp(focus, 'all') || strcmp(focus, 'psarti')
%     figure;
%     p_thre = 0.05;
%     cmat = nan(lenses, 4);
%     pmat = nan(lenses, 4);
%     for i = 1:lenses
%         cmat(i, :) = [data.session(i).psarti.original.Pearson.r, ...
%             data.session(i).psarti.reg(2).corrected.Pearson.r];
%         pmat(i, :) = [data.session(i).psarti.original.Pearson.p, ...
%             data.session(i).psarti.reg(2).corrected.Pearson.p];
%     end
%     xlabs = {'Pearson r: ps vs x', 'Pearson r: ps vs y'};
%     for k = 1:2
%         for i = 1:lenses
%             subplot(1,2,k)
%             % zero line
%             if i==1
%                 plot([0 0], [0 lenses+1], '-', 'color', 0.5*[1 1 1], 'linewidth', 0.1)
%             end
%             % connect line
%             hold on;
%             plot([cmat(i,k) cmat(i, k+2)], i*[1 1], '-k', 'linewidth', 0.1)
%             % uncorrected
%             hold on;
%             if pmat(i, k) < p_thre
%                 plot(cmat(i,k), i, '^', 'markerfacecolor', red, ...
%                     'markeredgecolor', red, 'markersize', 2)
%             else
%                 plot(cmat(i,k), i, '^', 'markerfacecolor', 'w', ...
%                     'markeredgecolor', red, 'markersize', 2)
%             end
%             % corrected
%             hold on;
%             if pmat(i, k+2) < p_thre
%                 plot(cmat(i,k+2), i, 'o', 'markerfacecolor', blue, ...
%                     'markeredgecolor', blue, 'markersize', 2)
%             else
%                 plot(cmat(i,k+2), i, 'o', 'markerfacecolor', 'w', ...
%                     'markeredgecolor', blue, 'markersize', 2)
%             end
%         end        
%         xlabel(xlabs{k})
%         if k==1
%             ylabel('session')
%         end
%         ylim([0.5 lenses+0.5])
%         set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
%     end
% 
%     % autosave figures
%     savefig(gcf, [figpath 'Figure3_PSartifact\raw_figs\ps_gz_corr' tasktype '.fig'])
% end

%%
% fixation precision (fix span, x var, y var, theta, SI)
if strcmp(focus, 'all') || strcmp(focus, 'fixSpan')
    paranames = {'fix span (log2(deg^2))', 'x variance (deg^2)', 'y variance (deg^2)', ...
            'theta (^o)', 'symmetry index'};
    fignames = {'no ps artifact correction', 'ps artifact corrected (glm)', 'ps artifact corrected (poly 2)'};
    cmap = [colorGradient([1 1 1], [0 0 1], 10); jet(90)]; 
    for c = 1:3
        para.fp(c).params = nan(5, lenses);
    end
    for c = 1:3
        if c==2
            continue
        end
        % fixation span parameters
        for i = 1:lenses
            para.fp(c).params(1, i) = data.session(i).fixPrecision{c}.accuracy;
            para.fp(c).params(2, i) = data.session(i).fixPrecision{c}.variance(1);
            para.fp(c).params(3, i) = data.session(i).fixPrecision{c}.variance(2);
            para.fp(c).params(4, i) = data.session(i).fixPrecision{c}.theta;
            para.fp(c).params(5, i) = data.session(i).fixPrecision{c}.SI;
        end

        % stats
        [~, pval, ~, st] = ttest(para.fp(c).params(2, :), para.fp(c).params(3, :));
        disp([num2str(c) '; x vs y variance: p(ttest) = ' num2str(pval) ', t = ' num2str(st.tstat)])
        [rval, pval] = corrcoef(para.fp(c).params(2, :), para.fp(c).params(3, :));
        disp([num2str(c) '; x vs y Pearson corr=' num2str(rval(1,2)) ...
            ', p(pearson) = ' num2str(pval(1,2))])
        disp(['median fix span = ' num2str(median(para.fp(c).params(1, :)))])
        
        % log2 transform
        para.fp(c).params(1, :) = log2(para.fp(c).params(1, :));
        para.fp(c).params(2, :) = log2(para.fp(c).params(2, :));
        para.fp(c).params(3, :) = log2(para.fp(c).params(3, :));
        
        % plot
        session_plot(para.fp(c).params, [], paranames, 'fixation precision')
        
        % autosave figures
        if saveoption==1
            savefig(gcf, [figpath 'Figure3_FixationPrecision\raw_figs\fixationPrecision' tasktype num2str(c) '.fig'])
        end
        
        % 2D hist of eye positions
        figure;
        nr = floor(sqrt(lenses));
        nc = ceil(lenses/nr);
        smoothr = 3;
        for i = 1:lenses
            subplot(nr, nc, i)
            m = smooth2a(data.session(i).fixPrecision{c}.Prob', smoothr, smoothr);
            imagesc(data.session(i).fixPrecision{c}.edge_x, data.session(i).fixPrecision{c}.edge_y, m)
            colormap(cmap)
%             caxis([0 0.01])
    %         xlim(max(data.session(i).params.fixwin(1,:))*[-1 1])
    %         ylim(max(data.session(i).params.fixwin(2,:))*[-1 1])
            axis(0.7*[-1 1 -1 1])
            title(num2str(i))
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        end
        set(gcf, 'Name', ['2D hist of eye positions: ' fignames{c}], 'NumberTitle', 'off')
        
        % autosave figures
        if saveoption==1
            savefig(gcf, [figpath 'Figure3_FixationPrecision\raw_figs\fixSpan_all' tasktype num2str(c) '.fig'])
        end
    end
    % before and after comparison
    figure;
    for i = 1:5
        subplot(1,5,i)
        unity_scatter(para.fp(1).params(i,:)', para.fp(end).params(i,:)')
        title(paranames{i})
        if i==1
            xlabel('before')
            ylabel('after')
        end
        axis square
        hold on;
    end    
    % autosave figures
    if saveoption==1
        savefig(gcf, [figpath 'Figure3_FixationPrecision\raw_figs\artifactcorr_beforeafter.fig'])
    end
end

%%
% microsaccade 
if strcmp(focus, 'all') || strcmp(focus, 'microsaccade')
    % load data
    load(['Z:\Katsuhisa\headfree_project\dataset\ms_uncor' tasktype '.mat'])
    ms = ms_uncor;
    lenses = length(ms.session); 
    
    foldername = 'Figure4_DiscriminationTask';
    
    % microsaccade parameters
    params = nan(5, lenses);
    errors = zeros(5, lenses);
    for i = 1:lenses
        params(1, i) = ms.session(i).results.counts/...
            (0.002*length(ms.session(i).results.event));
        params(2, i) = mean(ms.session(i).results.amp);
        params(3, i) = mean(ms.session(i).results.peakv);
        params(4, i) = mean(ms.session(i).results.duration);
        params(5, i) = mean(ms.session(i).results.angle);
        errors(2, i) = nanstd(ms.session(i).results.amp)/sqrt(ms.session(i).results.counts);
        errors(3, i) = nanstd(ms.session(i).results.peakv)/sqrt(ms.session(i).results.counts);
        errors(4, i) = nanstd(ms.session(i).results.duration)/sqrt(ms.session(i).results.counts);
        errors(5, i) = nanstd(ms.session(i).results.angle)/sqrt(ms.session(i).results.counts);
    end
    paranames = {'rate (/sec)', 'amplitude (deg)', 'peak velocity (deg/sec)', ...
        'duration (sec)', 'angle (^o)'};
    session_plot(params, errors, paranames, 'microsaccade parameters')
    
    % autosave figures
    if saveoption==1
        savefig(gcf, [figpath foldername '\raw_figs\ms_params' tasktype '.fig'])
    end
    
    % all sessions (amplitude vs peak velocity)
    figure;
    nr = floor(sqrt(lenses));
    nc = ceil(lenses/nr);    
    amp = []; peakv = [];
    for i = 1:lenses
        subplot(nr, nc, i)
        try
            plot(ms.session(i).results.amp, ms.session(i).results.peakv, 'o', 'color', 'k', 'markersize', 1)
            hold on;
            title(num2str(i))
%             set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            
            amp = [amp, ms.session(i).results.amp];
            peakv = [peakv, ms.session(i).results.peakv];
        catch
            continue
        end
    end
    set(gcf, 'Name', 'amp vs peakv', 'NumberTitle', 'off')

    % autosave figures
    if saveoption==1
        savefig(gcf, [figpath foldername '\raw_figs\ampVSpeakv_all' tasktype '.fig'])
    end
    
    % all collapsed
    figure; 
    [edgex, edgey, N] = ndhist(log(amp'), log(peakv'), 'bins', 5);
    cols = flipud(bone(100));
    colormap(cols)
    xlabel('amplitude (deg)', 'fontsize', 6)
    ylabel('peak velocity (deg/sec)', 'fontsize', 6)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    set(gcf, 'Name', 'amplitude vs peak velocity', 'NumberTitle', 'off')
    
    % autosave figures
    if saveoption==1
        savefig(gcf, [figpath foldername '\raw_figs\ampVSpeakv' tasktype '.fig'])
    end
end

%%
% discrimination task
if strcmp(focus, 'all') || strcmp(focus, 'afc')
    if isfield(data.session, 'pm')        
        foldername = 'Figure4_DiscriminationTask';
        
        % pschometric functions
        figure;
        nr = floor(sqrt(lenses));
        nc = ceil(lenses/nr);
        for i = 1:lenses
            subplot(nr, nc, i)
            try
                plot(data.session(i).pm.fitx, data.session(i).pm.fity, '-k', 'linewidth', 0.1)
                hold on;
                plot(data.session(i).pm.raw(1,:), data.session(i).pm.raw(2,:), 'o', ...
                    'markersize', 2.5, 'markerfacecolor', 'w', 'markeredgecolor', 'k')
                title(num2str(i))
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            catch
                continue
            end
        end
        set(gcf, 'Name', 'PMs', 'NumberTitle', 'off')
        
        % autosave figures
        if saveoption == 1
            savefig(gcf, [figpath foldername '\raw_figs\PM_all' tasktype '.fig'])
        end
        
        % reaction times
        figure;
        for i = 1:lenses
            try
                idx = abs(data.session(i).eyedata.reward)>0;
                signs = ones(1, length(data.session(i).eyedata.or));
                signs(data.session(i).eyedata.or<45) = -1;
                dcs = data.session(i).eyedata.dc.*signs;
                unidc = unique(dcs);
                lenuni = length(unidc);
                rts = nan(4, lenuni);
                for u = 1:lenuni
                    for c = 1:2
                        try
                            if c==1 % error
                                acc = data.session(i).eyedata.reward < 0;
                            else % correct
                                acc = data.session(i).eyedata.reward > 0;
                            end
                            rts(2*c-1, u) = nanmean(data.session(i).eyedata.rt(idx & acc & dcs==unidc(u)));
                            rts(2*c, u) = nanstd(data.session(i).eyedata.rt(idx & acc & dcs==unidc(u)))...
                                /sqrt(sum(idx & acc & dcs==unidc(u)));
                        catch
                            continue
                        end
                    end
                end
                subplot(nr, nc, i)
                errorbar(unidc, rts(1,:), rts(2,:), 'color', red, 'capsize', 0)
                hold on;
                errorbar(unidc, rts(3,:), rts(4,:), 'color', blue, 'capsize', 0)
                title(num2str(i))
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            catch
                continue
            end
        end
        set(gcf, 'Name', 'RTs', 'NumberTitle', 'off')
        
        % autosave figures
        if saveoption==1
            savefig(gcf, [figpath foldername '\raw_figs\RT_all' tasktype '.fig'])
        end
        
        % choice saccade
        figure;
        rng(19891220);
        nsample = 30;
        dpp = 0.0117;
        m = 100;
        for i = 1:lenses
            subplot(nr, nc, i)
            try
                idx = datasample(find(abs([data.session(i).eyedata.reward])>0), nsample);

                % windows
                fixwin = [max(data.session(i).params.fixwin(1,:)), max(data.session(i).params.fixwin(2,:))];
                plot(fixwin(1)*[-1 1 1 -1 -1], fixwin(2)*[-1 -1 1 1 -1], '-r', 'linewidth', 0.1)
                hold on;
                plot(m*dpp*[-1 1 1 -1 -1], m*dpp*[-1 -1 1 1 -1] + dpp*180, '-k', 'linewidth', 0.1)
                hold on;
                plot(m*dpp*[-1 1 1 -1 -1], m*dpp*[-1 -1 1 1 -1] - dpp*180, '-k', 'linewidth', 0.1)
                hold on;

                % choice saccades
                for k = 1:nsample
                    xv = data.session(i).eyedata.choicevec(idx(k)).x;
                    yv = data.session(i).eyedata.choicevec(idx(k)).y;
                    if abs(xv(1)) < fixwin(1) && abs(yv(1)) < fixwin(2) && ...
                            abs(xv(end)) < m*dpp && abs(yv(end)) < (m+180)*dpp
                        p = plot(xv(2:end-1), yv(2:end-1), '-', 'linewidth', 0.5, 'color', 0.5*[1 1 1]);
                        p.Color(4) = 0.1;
                        hold on;
                        plot(xv(1), yv(1), 'or', 'markersize', 2)
                        hold on;
                        plot(xv(end), yv(end), 'ok', 'markersize', 2)
                        hold on;
                    end
                end
                title(num2str(i))
                axis(3.5*[-1 1 -1 1])
                set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            catch
                continue
            end
        end
        set(gcf, 'Name', 'choice saccade', 'NumberTitle', 'off')
        
        % autosave figures
        if saveoption == 1
            savefig(gcf, [figpath foldername '\raw_figs\chvec_all' tasktype '.fig'])
        end
        
        % pupil size
        nt = 1000;
        out = zeros(1, lenses);
        for i = 1:lenses
            cand = length(data.session(i).ps.auc);
            disp(['session ' num2str(i) '; the number of frames: ' num2str(cand)])
            if cand < 900
                out(i) = 1;
                continue
            end
            if nt > cand
                nt = cand;
            end
        end
        for k = 1:2
            pstc.avrew(k).tc = nan(lenses, nt);
            pstc.avrew(k).ntr = zeros(lenses, 1);
        end
        for k = 1:5
            pstc.timebin(k).tc = nan(lenses, nt);
            pstc.timebin(k).ntr = zeros(lenses, 1);
        end
%         pstc.auc.tc = nan(lenses, nt);
        psych = nan(lenses, 2);
        psavrew = nan(lenses, 2);
        pstime = nan(lenses, 5);
        for i = 1:lenses
            for k = 1:2
                % time-course
                if out(i)==0
                    pstc.avrew(k).tc(i,:) = data.session(i).ps.avrew(k).tc.mean(end-nt+1:end);
                    pstc.avrew(k).ntr(i) = data.session(i).ps.avrew(k).tc.ntr;
                end
                % response
                lenps = length(data.session(i).ps.avrew(k).tc.mean);
                [minv, minidx] = min(data.session(i).ps.avrew(k).tc.mean(1:round(lenps/2)));
                psavrew(i, k) = max(data.session(i).ps.avrew(k).tc.mean(minidx+1:end)) ...
                    - minv;
            end
            for k = 1:5
                % time-course
                if out(i)==0
                    pstc.timebin(k).tc(i,:) = data.session(i).ps.timebin(k).tc.mean(end-nt+1:end);
                    pstc.timebin(k).ntr(i) = data.session(i).ps.timebin(k).tc.ntr;
                end
%                 
                % response
                lenps = length(data.session(i).ps.timebin(k).tc.mean);
%                 [minv, minidx] = min(data.session(i).ps.timebin(k).tc.mean(1:round(lenps/2)));
%                 pstime(i, k) = max(data.session(i).ps.timebin(k).tc.mean(minidx+1:end)) ...
%                     - minv;
                pstime(i, k) = mean(data.session(i).ps.timebin(k).tc.mean(end-round(lenps/2)+1:end));
            end
%             pstc.auc.tc(i,:) = data.session(i).ps.auc(1:nt);
            psych(i, :) = data.session(i).pm.params;
        end
        
        figure;
%         subplot(2,3,[1 2])
%         imagesc([1:nt]/500, 1:lenses, pstc.auc.tc)
%         colorbar('eastoutside')
%         xlabel('time after stimulus onset (sec)')
%         ylabel('sessions')
%         title('aROC (easy vs hard) in pupil size')
%         set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
                
        subplot(2,3,1)
        plot(1:lenses, abs(100*psych(:,1)), '-ok', 'linewidth', 0.1, ...
            'markersize', 3, 'markerfacecolor', 'k', 'markeredgecolor', 'k')
        xlim([1 lenses])
        xlabel('sessions')
        ylabel('|bias| (% signal)')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        
        subplot(2,3,2)
        plot(1:lenses, 100*psych(:,2), '-ok',  'linewidth', 0.1, ...
            'markersize', 3, 'markerfacecolor', 'k', 'markeredgecolor', 'k')
        pr = ranksum(100*psych(1:8,2), 100*psych(9:16,2));
        [r, p] = corr([1:lenses]', 100*psych(:,2));
        title({['p(ranksum) = ' num2str(pr)],['r(spe) = ' num2str(r) ', p(spe) = ' num2str(p)]})
        xlim([1 lenses])
        xlabel('sessions')
        ylabel('threshold (% signal)')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        
        subplot(2,3,3)
        idx = 8;
        t = [1:length(data.session(idx).ps.avrew(1).tc.mean)]/500;
        me = (data.session(idx).ps.avrew(1).tc.mean.*data.session(idx).ps.avrew(1).tc.ntr ...
            + data.session(idx).ps.avrew(2).tc.mean.*data.session(idx).ps.avrew(2).tc.ntr)/ ...
            (data.session(idx).ps.avrew(1).tc.ntr + data.session(idx).ps.avrew(2).tc.ntr);
        sem = (data.session(idx).ps.avrew(1).tc.sem.*data.session(idx).ps.avrew(1).tc.ntr ...
            + data.session(idx).ps.avrew(2).tc.sem.*data.session(idx).ps.avrew(2).tc.ntr)/ ...
            (data.session(idx).ps.avrew(1).tc.ntr + data.session(idx).ps.avrew(2).tc.ntr);
        fill_between(t, me - sem, me + sem, [0 0 0], 0.4)
        hold on;
        plot(t, me, '-', 'color', [0 0 0], 'linewidth', 0.5)
        hold on;
        [minv, minidx] = min(me);
        [maxv, maxidx] = max(me(minidx+1:end));
        plot([t(minidx)-0.1 t(maxidx+minidx+1)+0.1], minv*[1 1], '--k', 'linewidth', 0.1)
        hold on;
        plot(t(maxidx+minidx+1)*[1 1], [minv maxv], '-k', 'linewidth', 0.1)
        text(t(maxidx+minidx+1)+0.1, mean([minv maxv]), 'pupil response', 'fontsize', 6)
        xlim([1 nt]/500)
        xlabel('time after stimulus onset (sec)')
        ylabel('z-scored pupil size (a.u.)')
        title(['session ' num2str(idx)], 'fontsize', 6)
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        
        subplot(2,3,4)
        lent = length(data.session(end).ps.avrew(1).tc.mean);
        cols = {blue, red};
        for k = 1:2
            me = data.session(end).ps.avrew(k).tc.mean;
            sem = data.session(end).ps.avrew(k).tc.sem;
            fill_between([1:lent]/500, me - sem, me + sem, cols{k}, 0.4)
            hold on;
            plot([1:lent]/500, me, '-', 'color', cols{k}, 'linewidth', 0.5)
            hold on;
        end
        xlim([1 lent]/500)
        xlabel('time after stimulus onset (sec)')
        ylabel('z-scored pupil size (a.u.)')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        
        subplot(2,3,5)
%         minima = min(psavrew(:));
%         maxima = max(psavrew(:));
%         axrange = [minima - 0.05*(maxima - minima), ...
%             maxima + 0.05*(maxima - minima)];
        axrange = [0 0.6];
        plot(axrange, axrange, '-k', 'linewidth', 0.1)
        hold on;
        plot(psavrew(:,1), psavrew(:,2), 'o', 'markersize', 3, 'color', 'k')
        pval = signrank(psavrew(:,1), psavrew(:,2));
        text(axrange(1)+0.45*(axrange(2)-axrange(1)), axrange(1)+0.3*(axrange(2)-axrange(1)), ...
            ['p (signrank) '], 'fontsize', 6)
        text(axrange(1)+0.45*(axrange(2)-axrange(1)), axrange(1)+0.1*(axrange(2)-axrange(1)), ...
            ['= ' num2str(pval)], 'fontsize', 6)
        axis([axrange axrange])
        title('pupil size response', 'fontsize', 6)
        xlabel('small avrew.')
        ylabel('large avrew.')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        
        subplot(2,3,6)
        errorbar(1:5, nanmean(pstime, 1), nanstd(pstime, [], 1)/sqrt(lenses), ...
            '-k', 'capsize', 0)
        xlim([0.8 5.2])
        xlabel('time bin')
        ylabel('normalized pupil size')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
        
        % autosave figures
        set(gcf, 'Name', 'behavior_pupil', 'NumberTitle', 'off')
        if saveoption==1
            savefig(gcf, [figpath foldername '\raw_figs\beh_ps' tasktype '.fig'])
        end
%         subplot(2,3,3)
%         col = copper(5);
%         colors =  cell(1,5);
%         for t = 1:5
%                 colors{t} = col(t,:);
%         end
%         pstc_plot(pstc.timebin, colors, {'-','-','-','-','-'})   
%         ylabel('z-scored pupil size')
%         
%         subplot(2,3,6)
%         pstc_plot(pstc.avrew, {blue, red}, {'-','-'})    
%         xlabel('time after stimulus onset (sec)')
%         ylabel('z-scored pupil size')
    else
        disp('')
    end    
end

% subfunctions
function session_plot(params, errors, paranames, figtitle)
if isempty(errors)
    errors = zeros(size(params));
end
figure;
lenses = size(params, 2);
lenp = length(paranames);
firstn = 14;
for i = 1:lenp
    subplot(1, lenp, i)
    errorbar(1:lenses, params(i, :), errors(i,:), '-ok', 'linewidth', 0.1, ...
        'markersize', 1, 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'capsize', 0)
    xlim([1 lenses])
    xlabel('session')
    ylabel(paranames{i})
%     [rr,pp] = corrcoef([1:lenses]', params(i,:)');
%     [rr, pp] = corr([1:firstn]', params(i,1:firstn)', 'type', 'Spearman');
    pr = ranksum(params(i,1:23), params(i,24:end));
    [r, p] = corr([1:lenses]', params(i,:)', 'type', 'Spearman');
    title({['p(ranksum) = ' num2str(pr)], ...
        ['r(spe) = ' num2str(r) ' , p(spe) = ' num2str(p)]})
%     title({['r(1:' num2str(firstn) ') = ' num2str(rr) ' , p(1:' num2str(firstn) ') = ' num2str(pp)], ...
%         ['r(spe) = ' num2str(r) ' , p(spe) = ' num2str(p)]})
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end
set(gcf, 'Name', figtitle, 'NumberTitle', 'off')



function pstc_plot(tcs, colors, linetypes)
% data extraction
lenc = length(tcs);
nses = size(tcs(1).tc, 1);
nframe = size(tcs(1).tc, 2);
xvec = [1:nframe]/500;
l = zeros(1,lenc);
for c = 1:lenc
    hold on;
    [me, sem] = weighted(tcs(c).tc, tcs(c).ntr);
    fill_between(xvec, me-sem, me+sem, colors{c}, 0.4)
    hold on;
    l(c) = plot(xvec,me, 'color',colors{c}, 'linestyle', linetypes{c}, 'linewidth', 0.5);
end
% stats
p_vals = nan(2, nframe);
xlim([xvec(1) xvec(end)])
yy = get(gca, 'YLim');
if lenc > 1
    for t = 1:nframe
        statsmat = nan(nses, lenc);
        for c = 1:lenc
            statsmat(:,c) = tcs(c).tc(:, t);
        end
        if lenc > 2
            p_vals(1,t) = anova1(statsmat, [], 'off');
            p_vals(2,t) = kruskalwallis(statsmat, [], 'off');
        elseif lenc==2
            [~,p_vals(1,t)] = ttest(statsmat(:,1), statsmat(:,2));
            p_vals(2,t) = signrank(statsmat(:,1), statsmat(:,2));
        end
    end
    pvals = nan(size(p_vals));
    pvals(p_vals < 0.05/nframe) = 1;
    hold on;
    plot(xvec, (yy(1)+0.1*(yy(2)-yy(1)))*pvals(1,:), 'color', 0.8*ones(1,3), 'linewidth', 2);
end
hold on;
yy = get(gca, 'YLim');
yy(1) = floor(100*yy(1))/100;
yy(2) = ceil(100*yy(2))/100;
ylim(yy)
set(gca, 'XTick', [xvec(1) xvec(end)])
set(gca,'YTick', yy)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')