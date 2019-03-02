function format(fignum)
% format figures 

if nargin < 1; fignum = 'all'; end
close all;
gr = (1 + sqrt(5))/2; % golden ratio to 1
fz = 6; % font size
figpath = 'Z:\Katsuhisa\headfree_project\figures\';
addpath(genpath('Z:\Katsuhisa\code\integrated\matlab_usefulfunc'))
addpath(genpath('Z:\Katsuhisa\code\integrated\cbrewer'))

%%
% Fig. 1: head-free setup
% a: cup schematics
% b: movie of training
% c: example eye traces
% 
if strcmp(fignum, 'all') || fignum==1
    % figure format convention
    [figPars, axPars] = setPlotPars;
    figPos = [5 5 21 20]; % this is in cm
    figure(figPars);
    set(gcf, 'position', figPos, 'paperposition', figPos);

    % spacing parameters
    xbegin = 2.5;
    ybegin = 5;
    sq = 1.5;
    offset_figlab = [1.8, -0.5];
    figspace_x = 4;
    figspace_y = 3;
    
    % foldername
    foldername = 'Figure1_HeadFreeSetup';
    
    % a
    ax_new = axes(axPars, 'position', [xbegin ybegin+2*figspace_y 2*gr*sq sq]);
    fig = openfig([figpath foldername '\raw_figs\advantage.fig'], 'invisible');
    % axis object
    axesObjs = get(fig, 'Children');
    xlabs = {'head-post surgery', 'waiting period (a few months) to ensure logevity of the implant', 'start training', 'with head-fixation'};
    copyobj(axesObjs.Children, ax_new); delete(fig);
    set(gca, 'XTick', [1 3 5], 'XTickLabel', xlabs, 'YTick', [])
    ylim([0.05 0.3])
    offset_axis([0.05 0.05], axPars)
    
    % a panel number
    axes(axPars,'position', [xbegin-offset_figlab(1) ybegin+2*figspace_y-offset_figlab(2)...
        2*gr*sq sq])
    title('a','fontsize',8)
    axis off
    
    % c
    idx = {'2015.09.04', 5};
    load(['Z:\Katsuhisa\headfree_project\dataset\eyes\' idx{1} '_eyemat.mat']);
    vnames = {'x position (deg)', 'y position (deg)', 'pupil size (a.u.)'};
    ylims = {[-0.2 0.2], [-0.4 0.4], [-0.3 0.3]};
    cols = cbrewer('qual', 'Accent', 3);
    for l = 1:3
        % eye traces
        itr = eyemat(end, :) == idx{2};
        t = [1:sum(itr)]/500;
        ax_new = axes(axPars, 'position', [xbegin+figspace_x/2 ybegin+figspace_y*(3-l)/3 2*gr*sq sq/3]);
        plot(ax_new, t, eyemat(l, itr), '-', 'color', cols(l,:), 'linewidth', 0.75)
        text(ax_new, -1, 0, vnames{l}, 'fontsize', fz)
        if l < 3
            set(gca, 'XTick', [])
        else
            set(gca, 'XTick', [0 1 2])
        end        
        ylim(ylims{l})
        set(gca, 'YTick', [ylims{l}(1) 0 ylims{l}(2)])
        set(ax_new, 'box', 'off'); set(ax_new, 'TickDir', 'out')
        offset_axis([0.05 0.05], axPars)
    end    
    xlabel({'time after fixation start', '(sec)'}, 'fontsize', fz)
    
    % c panel number
    axes(axPars,'position', [xbegin-offset_figlab(1) ybegin+figspace_y/3-offset_figlab(2) gr*sq sq])
    title('c','fontsize',8)
    axis off
    
    % autosave
    savefig([figpath foldername '\formatted_figs\fig_paper.fig'])
    print(gcf,'-dpdf',[figpath foldername '\formatted_figs\fig_paper.pdf'], sprintf('-r%d',300))
    disp([foldername ' saved!'])
end

%%
% Fig. 2: experimental parameters
% a: fixation task schematics
% b: fixation duration per trial
% c: fixation window area
% d: working duration
% e: example survivor function and 
%     scale parameters (all)
if strcmp(fignum, 'all') || fignum==2
    % figure format convention
    [figPars, axPars] = setPlotPars;
    figPos = [5 5 21 20]; % this is in cm
    figure(figPars);
    set(gcf, 'position', figPos, 'paperposition', figPos);

    % spacing parameters
    xbegin = 1.5;
    ybegin = 5;
    sq = 1.5;
    offset_figlab = [1.8, -0.5];
    figspace_x = 4;
    figspace_y = 3;
    xrange = [1 46];
    
    % foldername
    foldername = 'Figure2_ExperimentalParams';
    pnames = {'b', 'c', 'f', 'd', '', 'e'};
        
    % a
    axes(axPars,'position', [xbegin-offset_figlab(1) ybegin+figspace_y-offset_figlab(2) gr*sq sq])
    title('a','fontsize',8)
    axis off
    
    % b, c
    labs = {{'fixation duration per trial', 'required to get a reward', '(sec)'}, ...
        {'area of fixation window', '(deg^2)'}};
    for k = 1:length(labs)
        place = [xbegin+figspace_x*k ybegin+figspace_y gr*sq sq];
        ax_new = axes(axPars, 'position', place);
        fig = openfig([figpath foldername '\raw_figs\expparams_fix.fig'], 'invisible');
        % axis object
        axesObjs = get(fig, 'Children');
        copyobj(axesObjs(4-k).Children, ax_new); delete(fig);
        xlabel('session', 'fontsize', fz)
        ylabel(labs{k}, 'fontsize', fz)
        xlim(xrange)
        yy = get(gca, 'YLim');
        yy(2) = ceil(yy(2));
        set(gca, 'XTick', xrange, 'YTick', [yy(1) mean(yy) yy(2)])
        ylim(yy)
        set(gca,'DefaultTextFontSize',fz)
        offset_axis([0.05 0.05], axPars)
        
        % panel legend
        place(1:2) = place(1:2) - offset_figlab;
        axes(axPars,'position', place)
        title(pnames{k},'fontsize',8)
        axis off
    end   
        
    % f (left)
    ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin gr*sq sq]);
    fig = openfig([figpath foldername '\raw_figs\survival_all_fix.fig'], 'invisible');
    % axis object
    axesObjs = get(fig, 'Children');
    lena = length(axesObjs);
    ind = [16 27 36];     
    cols = cbrewer('qual', 'Set2', length(ind));
    for l = 1:length(ind)
        % survivor function
        stairs(ax_new, axesObjs(lena - ind(l) + 1).Children(2).XData, ...
            axesObjs(lena - ind(l) + 1).Children(2).YData, '-', 'color', ...
            cols(l,:), 'linewidth', 0.1)
        hold(ax_new, 'on');
        % fit
        p = plot(ax_new, axesObjs(lena - ind(l) + 1).Children(1).XData, ...
            axesObjs(lena - ind(l) + 1).Children(1).YData, '-', 'color', ...
            cols(l,:), 'linewidth', 1);
        p.Color(4) = 0.5;
        hold(ax_new, 'on');
    end
    delete(fig);
    set(ax_new, 'box', 'off'); set(ax_new, 'TickDir', 'out')
    xlabel({'time after fixation start', '(sec)'}, 'fontsize', fz)
    ylabel('P(keep fixation)', 'fontsize', fz)    
    text(1.3, 0.95, 'session:', 'fontsize', 6)
    for l = 1:length(ind)
        text(1.5, 0.95-0.1*l, num2str(ind(l)), 'color', cols(l,:), 'fontsize', 6)
    end
    set(ax_new, 'XTick', [0 1 2], 'YTick', [0.5 1])
    set(ax_new,'DefaultTextFontSize',fz)
    offset_axis([0.05 0.05], axPars)
    
    % panel legend
    axes(axPars,'position', [xbegin+figspace_x-offset_figlab(1) ybegin-offset_figlab(2) gr*sq sq])
    title(pnames{3},'fontsize',8)
    axis off
    
    % e, f (right)
    labs = {'fixation break (%)', {'median of', 'survivor function','(sec)'}, ...
        'working duration (min)'};
    idx = [3, 2, 4];
    for k = 1:length(labs)
        if k==1
            place = [xbegin+figspace_x*3 ybegin+figspace_y gr*sq sq];
        elseif k==2            
            place = [xbegin+figspace_x*2 ybegin gr*sq sq];
        elseif k==3
            place = [xbegin ybegin gr*sq sq];
        end
        ax_new = axes(axPars, 'position', place);
        fig = openfig([figpath foldername '\raw_figs\ses_behav_fix.fig'], 'invisible');
        % axis object
        axesObjs = get(fig, 'Children');
        copyobj(axesObjs(idx(k)).Children, ax_new); delete(fig);
        if k==2
            for l = 1:length(ind)                
                hold(ax_new, 'on')
                plot(ind(l)*[1 1], [3 4], '-', 'color', cols(l,:), 'linewidth', 0.1)
            end
        end
        xlabel('session', 'fontsize', fz)
        ylabel(labs{k}, 'fontsize', fz)        
        xlim(xrange)
        yy = round(get(gca, 'YLim'));
        if k==2
            yy(2) = 4;
        end
        set(gca, 'XTick', xrange, 'YTick', [yy(1) mean(yy) yy(2)])
        ylim(yy)
        set(gca,'DefaultTextFontSize',fz)
        offset_axis([0.05 0.05], axPars)
        
        % panel legend
        if ismember(k, [1 3])
            place(1:2) = place(1:2) - offset_figlab;
            axes(axPars,'position', place)
            title(pnames{3+k},'fontsize',8)
            axis off
        end
    end

    % autosave
    savefig([figpath foldername '\formatted_figs\fig_paper.fig'])
    print(gcf,'-dpdf',[figpath foldername '\formatted_figs\fig_paper.pdf'], sprintf('-r%d',300))
    disp([foldername ' saved!'])
end

% %%
% % Fig. 3: eye signal confounded by pupil size
% % a: example scatter (eye vs ps; uncorrected)
% % b: example scatter (eye vs ps; corrected)
% % c, d: all the pearson r (uncorrected - corrected)
% if strcmp(fignum, 'all') || fignum==3
%     % figure format convention
%     [figPars, axPars] = setPlotPars;
%     figPos = [5 5 21 20]; % this is in cm
%     figure(figPars);
%     set(gcf, 'position', figPos, 'paperposition', figPos);
% 
%     % spacing parameters
%     xbegin = 2.5;
%     ybegin = 5;
%     sq = 1.5;
%     offset_figlab = [1.8, -0.5];
%     figspace_x = 3.5;
%     figspace_y = 2.75;
%     
%     % foldername
%     foldername = 'Figure3_PSartifact';
%     panels = {'a','', 'b', 'c', 'd'};
%     
%     % a, b
%     idx = '2015.09.03';
%     load(['Z:\Katsuhisa\headfree_project\dataset\eyes\' idx '_eyemat.mat']);
%     eyemat = trav_eyemat(eyemat);
%     fixwin = 2;
%     x = linspace(min(eyemat(3,:)), max(eyemat(3,:)), 100);
%     exs = [1,2,6,7];
%     xx = 2*[-1 1];
%     yy = fixwin*[-1 1];
%     for l = 1:4
%         % eye traces
%         if l < 3
%             place = [xbegin+figspace_x*(l-1) ybegin+figspace_y sq sq];
%         else
%             place = [xbegin+figspace_x*(l-3) ybegin sq sq];
%         end
%         ax_new = axes(axPars, 'position', place);
%         s = exs(l);
%         eyetrace = eyemat(s,:);
%         pupil = eyemat(3,:);
% %         out = abs(eyetrace) > fixwin;
% %         eyetrace(out) = []; 
% %         pupil(out) = [];
% %         plot(ax_new, pupil, eyetrace, 'ok', 'markersize', 1)
%         scatter(ax_new, pupil, eyetrace, 1, 'o', 'filled', 'markerfacecolor', 'k', ...
%             'markerfacealpha', 0.1)
%         p = polyfit(eyemat(3, :), eyemat(s, :), 1);
%         hold(ax_new, 'on');
%         plot(x, polyval(p, x), '-r', 'linewidth', 0.5)
%         r = corrcoef(eyemat(3, :), eyemat(s, :));
%         r = round(1000*r(1,2))/1000;
%         text(ax_new, xx(1)+0.95*(xx(2)-xx(1)), ...
%             yy(1)+0.9*(yy(2)-yy(1)), ['r = ' num2str(r)], ...
%             'color', 'r', 'fontsize', fz)
%         axis([xx yy])
%         if l==1
%             ylabel({'x position', '(deg)'}, 'fontsize', fz)
%             title('uncorrected', 'fontsize', fz)
%         elseif l==2
%             ylabel({'y position', '(deg)'}, 'fontsize', fz)
%         elseif l==3
%             title('corrected', 'fontsize', fz)
%             xlabel({'pupil size', '(a.u.)'}, 'fontsize', fz)
%         end
%         set(ax_new, 'box', 'off'); set(ax_new, 'TickDir', 'out')
%         offset_axis([0.05 0.05], axPars)
%         
%         % panel legend
%         if mod(l, 2)==1
%             place(1:2) = place(1:2) - offset_figlab;
%             axes(axPars,'position', place)
%             title(panels{l},'fontsize',8)
%             axis off
%         end
%     end    
%     
%     % c, d
%     for l = 1:2
%         place = [xbegin+figspace_x*(l+1) ybegin 1.2*sq 1.2*gr*sq];
%         ax_new = axes(axPars, 'position', place);
%         fig = openfig([figpath foldername '\raw_figs\ps_gz_corr_fix.fig'], 'invisible');
%         % axis object
%         axesObjs = get(fig, 'Children');
%         copyobj(axesObjs(3-l).Children, ax_new); delete(fig);
%         if l==1
%             xlabel('pupil size vs horizontal position', 'fontsize', fz)
%             ylabel('session', 'fontsize', fz)
%         elseif l==2
%             xlabel('pupil size vs vertical position', 'fontsize', fz)
%         end
%         xlim([-0.6 0.6])
%         ylim([1 50])
%         set(gca, 'XTick', [-0.6 0 0.6], 'YTick', [1 50])
%         set(ax_new, 'box', 'off'); set(ax_new, 'TickDir', 'out')
%         offset_axis([0.05 0.05], axPars)
% 
%         % panel legend
%         place(1:2) = place(1:2) - offset_figlab;
%         axes(axPars,'position', place)
%         title(panels{3+l},'fontsize',8)
%         axis off
%     end
%     
%     % autosave
%     savefig([figpath foldername '\formatted_figs\fig_paper.fig'])
%     print(gcf,'-dpdf',[figpath foldername '\formatted_figs\fig_paper.pdf'], sprintf('-r%d',300))
%     disp([foldername ' saved!'])
% end

%%
% Fig. 3: fixation precision
% a: fixation presion with example 2d hist (uncorrected)
% b: fixation presion with example 2d hist (corrected)
if strcmp(fignum, 'all') || fignum==4
   % figure format convention
    [figPars, axPars] = setPlotPars;
    figPos = [5 5 21 20]; % this is in cm
    figure(figPars);
    set(gcf, 'position', figPos, 'paperposition', figPos);

    % spacing parameters
    xbegin = 2.5;
    ybegin = 5;
    sq = 1.5;
    offset_figlab = [1.8, -0.5];
    figspace_x = 4;
    figspace_y = 3;
    xrange = [1 46];
    
    % color
    cmap = [colorGradient([1 1 1], [0 0 1], 10); jet(90)]; 
    
    % foldername
    foldername = 'Figure3_FixationPrecision';
    pnames = {'a', 'b'};
    
    % examples
    ind = [6, 18, 29, 38];
    
    % a
    place = [xbegin ybegin gr*sq sq];
    ax_new = axes(axPars, 'position', place);
    fig = openfig([figpath foldername '\raw_figs\fixationPrecision_fix1.fig'], 'invisible');
    % axis object
    axesObjs = get(fig, 'Children');
    x = axesObjs(4).Children.YData;
    y = axesObjs(4).Children.YData;
    delete(fig);
    yy = [-6 1];
    offset = 0.5;
    cols = cbrewer('qual', 'Accent', 3);
    plot(ax_new, [1:length(x)]-offset, x, '-o','linewidth', 0.1, 'color', cols(1,:), ...
        'markersize', 1, 'markerfacecolor', cols(1,:), 'markeredgecolor', cols(1,:))
    hold(ax_new, 'on');
    plot(ax_new, [1:length(y)]+offset, y, '-o','linewidth', 0.1, 'color', cols(2,:), ...
        'markersize', 1, 'markerfacecolor', cols(2,:), 'markeredgecolor', cols(2,:))
    xlabel('session', 'fontsize', fz)
    ylabel('variance (deg^2)', 'fontsize', fz)
    text(20, 0.5,'x position', 'color', cols(1,:), 'fontsize', fz)
    text(20, -1,'y position', 'color', cols(2,:), 'fontsize', fz)
    xlim([xrange(1)-offset xrange(2)+offset] )
    ylim(yy)
    set(gca, 'XTick', xrange, 'YTick', [yy(1) -3 yy(2)], 'YTickLabel', {'2^{-6}','2^{-3}','2^1'})
    offset_axis([0.05 0.05], axPars)

    % panel legend
    place(1:2) = place(1:2) - offset_figlab;
    axes(axPars,'position', place)
    title(pnames{1},'fontsize',8)
    axis off
    
    % b
    yy = [-3 4];
    for k = 1:2
        if k==2
            break
        end
        % sessions
        place = [xbegin+figspace_x*k ybegin gr*sq sq];
        ax_new = axes(axPars, 'position', place);
        fig = openfig([figpath foldername '\raw_figs\fixationPrecision_fix' num2str(2*(k-1)+1) '.fig'], 'invisible');
        % axis object
        axesObjs = get(fig, 'Children');
        copyobj(axesObjs(end).Children, ax_new); delete(fig);
        for l = 1:length(ind)
            hold(ax_new, 'on');
            plot(ind(l)*[1 1], [2.5 5], '-k', 'linewidth', 0.1)
        end       
        xlabel('session', 'fontsize', fz)
        if k==1
%             title('uncorrected', 'fontsize', fz)
            ylabel('fixation span (deg^2)', 'fontsize', fz)
        elseif k==2
            title('corrected', 'fontsize', fz)
        end
        xlim(xrange)
        ylim(yy)
        set(gca, 'XTick', xrange, 'YTick', [yy(1) 0 yy(2)], 'YTickLabel', {'2^{-3}','2^0','2^4'})
        offset_axis([0.05 0.05], axPars)
        
        % panel legend
        place(1:2) = place(1:2) - offset_figlab;
        axes(axPars,'position', place)
        title(pnames{k+1},'fontsize',8)
        axis off
        
        % examples
        for l = 1:length(ind)
            place = [xbegin+figspace_x*k+(l-1)*figspace_x/3 ybegin+figspace_y sq/3 sq/3];
            ax_new = axes(axPars, 'position', place);
            fig = openfig([figpath foldername '\raw_figs\fixSpan_all_fix' num2str(2*(k-1)+1) '.fig'], 'invisible');
            % axis object
            axesObjs = get(fig, 'Children');
            copyobj(axesObjs(end-ind(l)+1).Children, ax_new); delete(fig);
            colormap(cmap)
%             caxis
            caxis([0 0.0008])
            title(num2str(ind(l)), 'fontsize', fz)
            axis(0.8*[-1 1 -1 1])
            if k > 1 || l > 1
                set(gca, 'XTick', [], 'YTick', [])
            end
            set(gca,'DefaultTextFontSize',fz)
%             offset_axis(0, axPars)
        end
        
    end   
       
    % autosave
    savefig([figpath foldername '\formatted_figs\fig_paper.fig'])
    print(gcf,'-dpdf',[figpath foldername '\formatted_figs\fig_paper.pdf'], sprintf('-r%d',300))
    disp([foldername ' saved!'])
end

%%
% Fig. 4: discrimination task (part 1)
% a: discrimination task schematics
% b: choice saccades examples
% c: psychophysical threshold with example PMs
% d: pupil size (mean, sem)
% e: pupil size (available reward size)
if strcmp(fignum, 'all') || fignum==5
    % figure format convention
    [figPars, axPars] = setPlotPars;
    figPos = [5 5 21 20]; % this is in cm
    figure(figPars);
    set(gcf, 'position', figPos, 'paperposition', figPos);

    % spacing parameters
    xbegin = 2;
    ybegin = 4;
    sq = 1.5;
    offset_figlab = [1.8, -0.5];
    figspace_x = 4;
    figspace_y = 3.5;
    
    % foldername
    foldername = 'Figure4_DiscriminationTask';
    
    % a
    axes(axPars,'position', [xbegin ybegin+4*figspace_y-offset_figlab(2) gr*sq sq])
    title('a','fontsize',8)
    axis off
    
    % b
    ind_ch = [3, 6, 12, 15];
    leni = length(ind_ch);
    mag = 1;
    for i = 1:leni
        place = [xbegin+figspace_x+(figspace_x/2)*(i-1) ybegin+3.5*figspace_y mag*sq mag*sq];
        ax_new = axes(axPars, 'position', place);
        fig = openfig([figpath foldername '\raw_figs\chvec_all_dis.fig'], 'invisible');
        % axis object
        axesObjs = get(fig, 'Children');
        copyobj(axesObjs(end-ind_ch(i)+1).Children, ax_new); delete(fig);
        axis(3.5*[-1 1 -1 1])
        title(num2str(ind_ch(i)), 'fontsize', fz)
        axis off
    end    
    % panel legend
    axes(axPars,'position', [xbegin+figspace_x-offset_figlab(1) ...
        ybegin+4*figspace_y-offset_figlab(2) gr*sq sq])
    title('b','fontsize',8)
    axis off
    
    % c (example pms)
    ind_pm = [7, 11, 13, 16];
    leni = length(ind_pm);
    mag = 0.5;
    for i = 1:leni
        place = [xbegin+(figspace_x/2)*(i-1) ybegin+3*figspace_y mag*sq mag*sq];
        ax_new = axes(axPars, 'position', place);
        fig = openfig([figpath foldername '\raw_figs\PM_all_dis.fig'], 'invisible');
        % axis object
        axesObjs = get(fig, 'Children');
        copyobj(axesObjs(end-ind_pm(i)+1).Children, ax_new); delete(fig);
        title(num2str(ind_pm(i)), 'fontsize', fz)
        set(gca, 'YTick', [0 1])
        offset_axis([0.05 0.05], axPars)
    end   
        
    % c (threshold)
    ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin+2*figspace_y gr*sq sq]);
    fig = openfig([figpath foldername '\raw_figs\beh_ps_dis.fig'], 'invisible');
    % axis object
    axesObjs = get(fig, 'Children');
    copyobj(axesObjs(5).Children, ax_new); delete(fig);
    xlabel('sessions', 'fontsize', fz)
    ylabel({'psychophysical threshold', '(% signal)'}, 'fontsize', fz)
    xlim([1 16])
    ylim([10 125])
    set(gca, 'XTick', [1 16], 'YTick', [10 100])
    offset_axis([0.05 0.05], axPars)
    
    % panel legend
    axes(axPars,'position', [xbegin+figspace_x-offset_figlab(1) ...
        ybegin+2*figspace_y-offset_figlab(2) gr*sq sq])
    title('c','fontsize',8)
    axis off
    
    % autosave
    savefig([figpath foldername '\formatted_figs\fig_paper.fig'])
    print(gcf,'-dpdf',[figpath foldername '\formatted_figs\fig_paper.pdf'], sprintf('-r%d',300))
    disp([foldername ' saved!'])   
end

%% 
% Fig. 5: pupil size
if strcmp(fignum, 'all') || fignum==6
    % figure format convention
    [figPars, axPars] = setPlotPars;
    figPos = [5 5 21 20]; % this is in cm
    figure(figPars);
    set(gcf, 'position', figPos, 'paperposition', figPos);

    % spacing parameters
    xbegin = 2;
    ybegin = 4;
    sq = 1.5;
    offset_figlab = [1.8, -0.5];
    figspace_x = 4;
    figspace_y = 3.5;
    
    % foldername
    foldername = 'Figure5_PupilSize';
    
    ax_new = axes(axPars, 'position', [xbegin ybegin+figspace_y gr*sq sq]);
    fig = openfig([figpath foldername '\raw_figs\beh_ps_dis.fig'], 'invisible');
    % axis object
    axesObjs = get(fig, 'Children');
    copyobj(axesObjs(4).Children, ax_new); delete(fig);
    xlabel({'time after stimulus onset', '(sec)'}, 'fontsize', fz)
    ylabel({'z-scored pupil size', '(a.u.)'}, 'fontsize', fz)
    xx = get(gca, 'XLim');
    ylim(0.4*[-1 1])
    set(gca, 'XTick', xx, 'XTickLabel', [0 2], ...
        'YTick', [-0.4 0 0.4])
    offset_axis([0.05 0.05], axPars)
    
    % panel legend
    axes(axPars,'position', [xbegin-offset_figlab(1) ybegin+figspace_y-offset_figlab(2)...
        gr*sq sq])
    title('a','fontsize',8)
    axis off

    % e
    ax_new = axes(axPars, 'position', [xbegin+1.5*figspace_x ybegin+figspace_y sq sq]);
    fig = openfig([figpath foldername '\raw_figs\beh_ps_dis.fig'], 'invisible');
    % axis object
    axesObjs = get(fig, 'Children');
    copyobj(axesObjs(2).Children, ax_new); delete(fig);
    xlabel({'pupil response (a.u.)', '(small available reward)'}, 'fontsize', fz)
    ylabel({'pupil response (a.u.)', '(large available reward)'}, 'fontsize', fz)
    set(gca, 'XTick', [0 0.6], 'YTick', [0 0.6])
    axis([0 0.6 0 0.6])
    offset_axis([0.05 0.05], axPars)
    
    % panel legend
    axes(axPars,'position', [xbegin+1.5*figspace_x-offset_figlab(1) ybegin+figspace_y-offset_figlab(2)...
        gr*sq sq])
    title('b','fontsize',8)
    axis off
    
    % autosave
    savefig([figpath foldername '\formatted_figs\fig_paper.fig'])
    print(gcf,'-dpdf',[figpath foldername '\formatted_figs\fig_paper.pdf'], sprintf('-r%d',300))
    disp([foldername ' saved!'])   
end

%%
% Fig. 6: Microsaccade
% a: microsaccade example traces 
% b: position space, velocity space
% c: amplitude vs peak velocity
if strcmp(fignum, 'all') || fignum==7
    % figure format convention
    [figPars, axPars] = setPlotPars;
    figPos = [5 5 21 20]; % this is in cm
    figure(figPars);
    set(gcf, 'position', figPos, 'paperposition', figPos);

    % spacing parameters
    xbegin = 2;
    ybegin = 5;
    sq = 1.5;
    offset_figlab = [1.8, -0.5];
    figspace_x = 3;
    figspace_y = 3;
    
    % foldername
    foldername = 'Figure6_Microsaccade';
    
    % a, b
%     idx = '2016.01.07';
%     e = load(['Z:\Katsuhisa\headfree_project\dataset\eyes\' idx '_eyemat.mat']);
%     j = 50; stop = 0;
%     while stop==0
%         itr = e.eyemat(end, :) == j;
%         ms = microsaccade_detection(e.eyemat(6,itr), e.eyemat(7,itr), 500, 'engbert', 1);
%         disp(['trial ' num2str(j) '; ' num2str(ms.counts)  ' microsaccades'])
%         if ms.counts < 3
%             j = j + 1;
%             delete(gcf);
%         else
%             stop = 1;
%         end
%     end
%     fig = gcf;
%     exfigpath =  [figpath foldername '\raw_figs\ms_example_' idx '_' num2str(j) '.fig'];
%     savefig(fig, exfigpath)
%     delete(fig);
%     yy= {0.4*[-1 1], 0.7*[-1 1]};
%     for l = 1:4
%         if l < 3
%             place = [xbegin ybegin+figspace_y*(1+0.5*(2-l)) gr*sq sq/2];
%         else
%             place = [xbegin+figspace_x*(l-3) ybegin sq sq];
%         end          
%         ax_new = axes(axPars, 'position', place);
%         fig = openfig(exfigpath, 'invisible');
%         % axis object
%         axesObjs = get(fig, 'Children');
%         copyobj(axesObjs(8-l).Children, ax_new); delete(fig);
%         switch l
%             case 1
%                 ylabel('x (deg)', 'fontsize', fz)
%                 set(gca, 'XTick', [], 'YTick', [yy{l}(1) 0 yy{l}(2)])
%                 ylim(yy{l})
%             case 2
%                 xlabel('time after stimlus onset (sec)', 'fontsize', fz)
%                 ylabel('y (deg)', 'fontsize', fz)
%                 set(gca, 'XTick', [0 2], 'YTick', [yy{l}(1) 0 yy{l}(2)])
%                 ylim(yy{l})
%             case 3
%                 xlim(0.4*[-1 1])
%                 ylim(0.7*[-1 1])
%                 set(gca, 'XTick', 0.4*[-1 0 1], 'YTick', 0.7*[-1 0 1])
%                 xlabel('x position (deg)', 'fontsize', fz)
%                 ylabel('y position (deg)', 'fontsize', fz)
%             case 4
%                 xlim(35*[-1 1])
%                 ylim(60*[-1 1])
%                 set(gca, 'XTick', 35*[-1 0 1], 'YTick', 60*[-1 0 1])
%                 xlabel('x velocity (deg/sec)', 'fontsize', fz)
%                 ylabel('y velocity (deg/sec)', 'fontsize', fz)
%         end
%         offset_axis([0.05 0.05], axPars)
%         
%         % panel legend
%         if l == 1
%             place(1:2) = place(1:2) - offset_figlab;
%             axes(axPars,'position', place)
%             title('a','fontsize',8)
%             axis off
%         end
%     end

    % a
    listses = {'2016.01.07_48', '2016.01.02_19'};
    for k = 1:length(listses)
        place = [xbegin ybegin+figspace_y*(2-k) gr*sq sq/2];
        ax_new = axes(axPars, 'position', place);
        fig = openfig([figpath foldername '\raw_figs\' listses{k} '.fig'], 'invisible');
        % axis object
        axesObjs = get(fig, 'Children');
        copyobj(axesObjs(3).Children, ax_new); delete(fig);        
        xx = [0, 2];
        yy = 0.5*[-1 1];
        axis([xx yy])
        set(gca, 'XTick', [xx(1) 1 xx(2)], 'YTick', [yy(1) 0 yy(2)])
        offset_axis([0.05 0.05], axPars)
        
        % panel legend
        if k==1
            place(1:2) = place(1:2) - offset_figlab;
            axes(axPars,'position', place)
            title('a','fontsize',8)
            axis off
        else
            xlabel('time after stimulus onset (sec)', 'fontsize', fz)
            ylabel('eye position (deg)', 'fontsize', fz)
        end
    end
    
    % b
    mag = 1.2;
    place = [xbegin+1.5*figspace_x ybegin mag*sq mag*sq];
    ax_new = axes(axPars, 'position', place);
    fig = openfig([figpath foldername '\raw_figs\ampVSpeakv_uneye.fig'], 'invisible');
%     fig = openfig([figpath foldername '\raw_figs\ampVSpeakv_all_dis.fig'], 'invisible');
    % axis object
    axesObjs = get(fig, 'Children');
    cols = axesObjs(2).Colormap;
    copyobj(axesObjs(2).Children, ax_new); delete(fig);
    xlabel('amplitude (deg)', 'fontsize', fz)
    ylabel('peak velocity (deg/sec)', 'fontsize', fz)
    colormap(ax_new, cols)
    xx = [0, 1.1];
    yy = [0, 90];
    axis([xx yy])
    set(gca, 'XTick', [0 1], 'YTick', yy)    
    
    c = colorbar('eastoutside');
    cpos = c.Position;
    cpos(1) = 1.2*cpos(1);
    cpos(3) = 1.1*cpos(3);
    c.Position = cpos;
    c.AxisLocation = 'out';
    c.Box = 'off';
    c.FontSize = fz;
    c.TickDirection = 'out';
    c.Ticks = [1 80 160];
    c.Limits = [1 160];
    c.Label.String = 'the number of microsaccades';
    
    offset_axis([0.05 0.05], axPars)
    % panel legend
    place(1:2) = place(1:2) - offset_figlab;
    axes(axPars,'position', place)
    title('b','fontsize',8)
    axis off
       
    % autosave
    savefig([figpath foldername '\formatted_figs\fig_paper.fig'])
    print(gcf,'-dpdf',[figpath foldername '\formatted_figs\fig_paper.pdf'], sprintf('-r%d',300))
    disp([foldername ' saved!'])   
end