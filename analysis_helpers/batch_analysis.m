function batch_analysis
%%
% main analysis script
%
% - history bias ... learned or congenital?
% - relationship between history bias and pupil or derivative of pupil 
% - within trial and session (split into 3 in each session)
% - how bias and threshold changes 
%

% path dependent on environment 
if ispc
    mypath = 'Z:/';
else % cluster
    mypath = '/gpfs01/nienborg/group/';
end

% path of analysis scripts
addpath(genpath([mypath 'Katsuhisa/HNlab_meta']))
addpath(genpath([mypath 'Katsuhisa/code/integrated/matlab_usefulfunc']))

% load list of sessions
% list = listMaker({'kiwi', 'kaki', 'mango'});
ses_list = load([mypath 'data/ses_list.mat']);
list = ses_list.list;
num_animal = length(list);

% batch processing
n_split = 4;
for a = 1:num_animal
    % number of sessions
    n_ses = length(list{a});
    list_a = list{a};
    parfor n = 1:n_ses
        try
            % load ex-file
            stridx = strfind(list_a{n}, 'data');
            ex = load([mypath list_a{n}(stridx:end)]);
            ex = ex.ex;        

            % split ex
            exs = ex_split(ex, n_split);
            Ldata = cell(1, n_split + 1);

            % single session analysis
            Ldata{1} = single_session_analysis(ex);

            % within session split
            for i = 1:n_split
                % single session analysis
                Ldata{i+1} = single_session_analysis(exs{i}, 0);
            end        

            % autosave
            ses_saver(Ldata, mypath)
        catch
            disp(['ERR: session ' num2str(n)])
        end
    end
end

% subfunction
function exs = ex_split(ex, n_split)
ex.Trials = ex.Trials(abs([ex.Trials.Reward]) > 0);
ntr = length(ex.Trials);
q = floor(quantile(1:ntr, 1/n_split));
exs = cell(1, n_split);
begin = 1;
for n = 1:n_split
    ex_temp = ex;
    if n == n_split
        ex_temp.Trials = ex.Trials(begin:end);
    else
        ex_temp.Trials = ex.Trials(begin:begin + q - 1);
    end
    exs{n} = ex_temp;
    begin = begin + q;
end

function ses_saver(Ldata, mypath)
slash = strfind(Ldata{1}.dirname, '/');
fname = [Ldata{1}.filename(1:2) '/' Ldata{1}.dirname(slash(end)+1:end), '.', ...
    Ldata{1}.filename(1:end-4), '_Ldata.mat'];
save([mypath 'Katsuhisa/learning_project/data/' fname], 'Ldata')
disp([fname ' saved!'])