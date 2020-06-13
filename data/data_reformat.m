%% Description

%{

Loads all source localised data and saves into a single .mat file

Storage structure is a cell array
    (participants x conditions x pairs)
    (2 x 3 x 10)
Each cell contains a cell array
    (songs x 1)
    (5-6 x 1)
Each cell contains a matrix
    (sources x samples x epoch-trials)
    (148 x 768 x 20-25)

%}

%% Setup

data_dir = '../data_raw/';

% tm = interpersonal flow
% tmrv = interpersonal non-flow
% tmoc = occluded individual flow
cond_names = {'tm', 'tmrv', 'tmoc'};

out_file = 'data_sources.mat';

%% Get info about data-files

str = fileread('../data_raw/data_dates.txt');
data_dates = regexp(str, '\r\n|\r|\n', 'split');

nPairs = length(data_dates);

%% Initialise full data structure

% Constant values: 3 conditions; 2 participants per pair; nPairs pairs
data = cell(2, length(cond_names), nPairs);

%% Load each pair and store

for pair = 1 : nPairs
    disp(pair);
    
    pair_string = data_dates{pair};
    
    % Split string based on spaces into (in order)
    %   date, pIDs, tm-trials, rv-trials, oc-trials
    pair_info = regexp(pair_string, '\s', 'split');
    pIDs = regexp(pair_info{2}, ',', 'split');
    
    for p = 1 : length(pIDs) % for each participant
        
        pID = pIDs{p};
        
        for c = 1 : 3 % for each condition
            
            songs = str2num(pair_info{2+c}); % Working on a vector, not scalar, so DO NOT use str2double as suggested
            
            data{p, c, pair} = cell(length(songs), 1);
            
            for s = 1 : length(songs) % for each song
                
                filename = [pair_info{1} '_' pID '_T' num2str(songs(s)) '_' cond_names{c} '_PLAY_gp.mat'];
                
                % Load file and store
                tmp = load([data_dir filename]);
                data{p, c, pair}{s} = tmp.Value;
                
            end
            
        end
        
    end
    
end

%% Save altogether
% ~ 16 seconds

disp('saving');
tic;
save('data_all.mat', 'data', 'cond_names');
toc
disp('saved');
