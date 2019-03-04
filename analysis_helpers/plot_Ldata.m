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
            scatter(n-1, Ldata{1}.fitpm.params(i), 30, 'o', 'markerfacecolor', cols(n, :), ...
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
            scatter(n-1, Ldata{1}.fitse.beta(i), 30, 'o', 'markerfacecolor', cols(n, :), ...
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

% PS ========================
figure;
paramnames = {'pupil size', 'pupil derivative'};
lenp = length(paramnames);
cols = [0.9922    0.6824    0.3804; ...
    0.6706    0.8510    0.9137];
t = linspace(0, Ldata{1}.stmdur, size(Ldata{1}.pupils{1, 1}, 2));
for i = 1:lenp
    subplot(1, lenp, i)
    % easy
    idx = abs(Ldata{1}.fitpm.raw(1, :)) > 0.4;
    plot(t, mean(Ldata{1}.pupils{i, 1}(idx, :), 1), '-', 'color', cols(1, :))
    hold on;
    
    % hard
    idx = abs(Ldata{1}.fitpm.raw(1, :)) > 0 & abs(Ldata{1}.fitpm.raw(1, :)) < 0.15;
    plot(t, mean(Ldata{1}.pupils{i, 1}(idx, :), 1), '-', 'color', cols(2, :))
    hold on;
    
%     for n = 1:n_split + 1
%         
%     end
    % format
    subplot(1, lenp, i)
    xlabel('time from stimulus onset')
    ylabel(paramnames{i})
    set(gca, 'box', 'off', 'tickdir', 'out')
end