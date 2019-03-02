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
    mypath = '/gpfs01/nienborg/users/';
end

% path of analysis scripts
addpath(genpath([mypath 'Katsuhisa/code/HNlab_meta']))
addpath(genpath([mypath 'Katsuhisa/code/integrated/matlab_usefulfunc']))

% load list of sessions
ses_list = load([mypath 'data/ses_list.mat']);
list = ses_list.list;
num_animal = length(list);

% batch processing
n_split = 4;
animal = {'kaki', 'mango', 'kiwi'};
for a = 1:num_animal
    % number of sessions
    n_ses = length(list{a});
    for n = 1:n_ses
        % load ex-file
        stridx = strfind(list{a}{n}, 'data');
        ex = load([mypath list{a}{n}(stridx:end)]);
        ex = ex.ex;
        
        
        % split ex
        delta = quantile(ntr, 1/n_split);
        c = 1;
        for i = 1:n_split
            
        end        
    end
    % autosave 
    save([mypath 'Katsuhisa/learning_project/data/Ldata_' animal{a} '.mat'], 'Ldata', '-v7.3')
    disp(['Ldata for ' animal{a} ' is saved!'])
end