function plot_Ldata(Ldata)
%%
% visualization function for single Ldata from single_session_analysis
%

close all;
n_split = length(Ldata) - 1;
cols = [[0 0 0]; cool(n_split)];

% PK ========================
figure;
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
    
set(gcf, 'Name', 'PK', 'NumberTitle', 'off')

% PM ========================
figure;
paramnames = {'bias', 'threshold', 'lapse rate'};
lenp = length(paramnames);
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

set(gcf, 'Name', 'PM', 'NumberTitle', 'off')

% serial choice bias =============
figure;
paramnames = {'Ch_{n-1}', 'StmSign_{n-1}'};
lenp = length(paramnames);
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

set(gcf, 'Name', 'choice history bias', 'NumberTitle', 'off')

% PS ========================
figure;
paramnames = {'pupil size', 'pupil derivative'};
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
    
    % PM-like ==========================
    subplot(2, lenp, i+2)
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
    subplot(2, lenp, i)
    xlabel('time from stimulus onset')
    ylabel(paramnames{i})
    set(gca, 'box', 'off', 'tickdir', 'out')
    subplot(2, lenp, i+2)
    xlabel('stimulus')
    ylabel('a.u.')
    set(gca, 'box', 'off', 'tickdir', 'out')
end

set(gcf, 'Name', 'pupils', 'NumberTitle', 'off')