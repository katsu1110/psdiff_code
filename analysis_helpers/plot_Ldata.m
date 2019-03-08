function plot_Ldata(Ldata)
%%
% visualization function for single Ldata from single_session_analysis
%

close all;
n_split = length(Ldata) - 1;
cols = [[0 0 0]; cool(n_split)];

%%
% PK ========================
figure;
if n_split==0
    % psychophysical kernel
    plot(Ldata.hdx_range, Ldata.pk, '-', 'color', 'k', 'linewidth', 3)
    hold on;
    plot([Ldata.hdx_range(1) Ldata.hdx_range(end)], [0 0], '--k')
    hold on;
    yy = get(gca, 'YLim');
    plot([0 0], yy, '--k')
    xlabel('hdx')
    ylabel('# frames')
    set(gca, 'box', 'off', 'tickdir', 'out')
else
    for n = 1:n_split+1
        % psychophysical kernel
        subplot(1, 2, 1)
        plot(Ldata{1}.hdx_range, Ldata{n}.pk, '-', 'color', cols(n, :))
        hold on;

        % amplitude
        subplot(1, 2, 2)
        if n == 1
            plot([1 n_split], sum(Ldata{1}.pk.*Ldata{n}.pk)*[1 1], '-', 'color', cols(n, :)) 
        else
            scatter(n-1, sum(Ldata{1}.pk.*Ldata{n}.pk), 30, 'o', 'markerfacecolor', cols(n, :), ...
                'markeredgecolor', [1 1 1 ])
        end
        hold on;
    end
    % format
    subplot(1, 2, 1)
    xlabel('hdx')
    ylabel('# frames')
    set(gca, 'box', 'off', 'tickdir', 'out')

    subplot(1, 2, 2)
    xlabel('time bin within session')
    ylabel('PKA')
    set(gca, 'box', 'off', 'tickdir', 'out')
end
    
set(gcf, 'Name', 'PK', 'NumberTitle', 'off')

%%
% PM ========================
paramnames = {'bias', 'threshold', 'lapse rate'};
lenp = length(paramnames);
if n_split > 0
    figure;
    for i = 1:lenp
        for n = 1:n_split+1
            subplot(1, lenp, i)
            if n == 1
                plot([1 n_split], Ldata{1}.fitpm.params(i)*[1 1], '-', 'color', cols(n, :)) 
            else
                scatter(n-1, Ldata{n}.fitpm.params(i), 60, 'o', 'markerfacecolor', cols(n, :), ...
                    'markeredgecolor', [1 1 1 ])
            end
            hold on;
        end
        % format
        subplot(1, lenp, i)
        xlabel('time bin within session')
        ylabel(paramnames{i})
        set(gca, 'box', 'off', 'tickdir', 'out')
    end
end

set(gcf, 'Name', 'PM', 'NumberTitle', 'off')

%%
% serial choice bias =============
figure;
paramnames = {'Ch_{n-1}', 'StmSign_{n-1}'};
lenp = length(paramnames);
if n_split == 0
    bar(1:lenp, Ldata.fitse.beta(1:lenp))
    set(gca, 'XTick', 1:lenp, 'XTickLabel', paramnames)
    xtickangle(45)
    set(gca, 'box', 'off', 'tickdir', 'out')
else
    for i = 1:lenp
        for n = 1:n_split+1
            subplot(1, lenp, i)
            if n == 1
                plot([1 n_split], Ldata{1}.fitse.beta(i)*[1 1], '-', 'color', cols(n, :)) 
            else
                scatter(n-1, Ldata{n}.fitse.beta(i), 60, 'o', 'markerfacecolor', cols(n, :), ...
                    'markeredgecolor', [1 1 1 ])
            end
            hold on;
        end
        % format
        subplot(1, lenp, i)
        xlabel('time bin within session')
        ylabel(paramnames{i})
        set(gca, 'box', 'off', 'tickdir', 'out')
    end
end

set(gcf, 'Name', 'choice history bias', 'NumberTitle', 'off')

%%
% PM & PK split ========================
if n_split == 0
    figure;
    paramnames = {'pupil size', 'pupil derivative', 'pupil pow'};
    lenp = length(paramnames);
    x = Ldata.fitpm.raw(1, :);
    lenx = length(x);
    v = [1:size(Ldata.mat, 1)]';
    for i = 1:lenp
        % median split
        idx1 = []; idx2 = [];
        for k = 1:lenx
            idx = find(Ldata.mat(:, 5).*sign(Ldata.mat(:,6)) == x(k));
            [idx1_temp, idx2_temp] = median_split(Ldata.mat(idx, 14+i));
            idx1 = [idx1; idx(idx1_temp)];
            idx2 = [idx2; idx(idx2_temp)];
        end
        idx1 = ismember(v, idx1); idx2 = ismember(v, idx2);
        
        % PM
        subplot(2, lenp, i)
        [x, y, n] = getPM(Ldata.mat(idx1, :));
%         data = fitPM(x,y,n,'Weibull','MLE',0);
        data = fitPM(x,y,n,'Gaussian','MLE',0);
        plot(data.fitx, data.fity, '-', 'color', 'c')
        hold on;
        plot(x, y, 'o', 'color', 'c')
        hold on;
        text(0, 0.4, 'bias, sensitivity, lapse rate', 'color', 'k', 'fontsize', 8)
%         text(0, 0.3, ['low: ' num2str(data.params(1)) ', ' num2str(data.params(2))...
%             ', ' num2str(data.params(3))], 'color', 'c', 'fontsize', 8)
        text(0, 0.3, ['low: ' num2str(data.params(1)) ', ' num2str(data.params(2))]...
            , 'color', 'c', 'fontsize', 8)
        [x, y, n] = getPM(Ldata.mat(idx2, :));
%         data = fitPM(x,y,n,'Weibull','MLE',0);
        data = fitPM(x,y,n,'Gaussian','MLE',0);
        plot(data.fitx, data.fity, '-', 'color', 'm')
        hold on;
        plot(x, y, 'o', 'color', 'm')    
%         text(0, 0.2, ['high: ' num2str(data.params(1)) ', ' num2str(data.params(2))...
%             ', ' num2str(data.params(3))], 'color', 'm', 'fontsize', 8)
        text(0, 0.2, ['high: ' num2str(data.params(1)) ', ' num2str(data.params(2))...
            ], 'color', 'm', 'fontsize', 8)
        xlabel('stimulus')
        ylabel('% far choice')
        title(paramnames{i})
        set(gca, 'box', 'off', 'tickdir', 'out')

        % PK
        subplot(2, lenp, i+lenp)
        zeroidx = abs(Ldata.mat(:, 5)) < 0.15;
        pk = getKernel(Ldata.stmmat(zeroidx & idx1, :), Ldata.hdx_range, ...
            Ldata.mat(zeroidx & idx1, 11));
        plot(Ldata.hdx_range, pk, '-', 'color', 'c')
        hold on;
        pk = getKernel(Ldata.stmmat(zeroidx & idx2, :), Ldata.hdx_range, ...
            Ldata.mat(zeroidx & idx2, 11));
        plot(Ldata.hdx_range, pk, '-', 'color', 'm')
        xlabel('hdx')
        ylabel('# frame')
        title(paramnames{i})
        set(gca, 'box', 'off', 'tickdir', 'out')
    end
else
    for i = 1:lenp
        for n = 1:n_split+1
            subplot(1, lenp, i)
            if n == 1
                plot([1 n_split], Ldata{1}.fitpm.params(i)*[1 1], '-', 'color', cols(n, :)) 
            else
                scatter(n-1, Ldata{n}.fitpm.params(i), 60, 'o', 'markerfacecolor', cols(n, :), ...
                    'markeredgecolor', [1 1 1 ])
            end
            hold on;
        end
        % format
        subplot(1, lenp, i)
        xlabel('time bin within session')
        ylabel(paramnames{i})
        set(gca, 'box', 'off', 'tickdir', 'out')
    end
end

set(gcf, 'Name', 'PM', 'NumberTitle', 'off')

%%
% PS ========================
figure;
if n_split == 0
    Ldata = {Ldata};
end
paramnames = {'pupil size', 'pupil derivative', 'pupil pow'};
lenp = length(paramnames);
cols = [1 0.7 0; ...
    0    0.7    0];
pmcols = [0 0 0; 0 0 1; 1 0 0];
t = linspace(0, Ldata{1}.stmdur, size(Ldata{1}.pupils{1, 1}, 2));
stm = Ldata{1}.mat(:,5).*sign(Ldata{1}.mat(:,6));
cor = Ldata{1}.mat(:,9);
for i = 1:lenp
    % stimulus types
    s = Ldata{1}.fitpm.raw(1, :);
    
   % time-course =================
   if i < 3
        subplot(2, lenp, i)
        % easy
        idx = abs(s) > 0.4;
        if sum(idx)==0
            idx = [1, length(s)];
        end
        plot(t, mean(Ldata{1}.pupils{i, 1}(idx, :), 1), '-', 'color', cols(1, :))
        hold on;

        % hard
        idx = abs(s) > 0 & abs(s) < 0.15;
        if sum(idx)==0
            idx = [floor(length(s)/2), floor(length(s)/2)+2];
        end
        plot(t, mean(Ldata{1}.pupils{i, 1}(idx, :), 1), '-', 'color', cols(2, :))
        hold on;
        
        % format
        subplot(2, lenp, i)
        xlabel('time from stimulus onset')
        ylabel(paramnames{i})
        set(gca, 'box', 'off', 'tickdir', 'out')
        legend('easy', 'hard', 'location', 'southeast')
        legend('boxoff')
   end
    
    % PM-like ==========================
    subplot(2, lenp, i+lenp)
    pmv = nan(3, length(s));
    for k = 1:length(s)
        pmv(1, k) = nanmean(Ldata{1}.mat(stm==s(k), 14+i));
        pmv(2, k) = nanmean(Ldata{1}.mat(stm==s(k) & cor==1, 14+i));
        pmv(3, k) = nanmean(Ldata{1}.mat(stm==s(k) & cor==0, 14+i));
    end    
    for k = 1:3
        plot(s, pmv(k,:), '-o', 'color', pmcols(k, :))
        hold on;
    end
    
%     for n = 1:n_split + 1
%         
%     end
    % format
    subplot(2, lenp, i+lenp)
    xlabel('stimulus')
    ylabel('a.u.')
    set(gca, 'box', 'off', 'tickdir', 'out')
    legend('all', 'correct', 'error', 'location', 'northwest')
    legend('boxoff')
end

set(gcf, 'Name', 'pupils', 'NumberTitle', 'off')

% subfuncton
function [x, y, n] = getPM(behmat)
stm = behmat(:, 5).*sign(behmat(:, 6));
x = unique(stm)';
lenu = length(x);
y = x;
n = x;
for i = 1:lenu
    % behavior
    cond = stm == x(i);
    n(i) = sum(cond);
    y(i) = sum(behmat(cond, 11)==1)/n(i);
end

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

function [idx1, idx2] = median_split(v)
[~, idx] = sort(v, 'ascend');
nlen = length(idx);
idx1 = idx(1:floor(nlen/2));
idx2 = idx(floor(nlen/2)+1:end);