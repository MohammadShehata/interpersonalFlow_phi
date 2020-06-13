%% Description

%{

Joins participant pairs, conditions, and timescales into a single file

%}

%% Constants

nChannels = 4; nStates = nChannels^2;
ps = (1:10);
taus = [20 40];

% TODO - load from file
sets = cell(4, 2); % networks x participants
sets(:, 1) = {[6 5 12], [13], [1 2 6], [11]};
sets(:, 2) = {[13], [6 5 12], [11], [1 2 6]};

cond_names = {'tm', 'tmrv', 'tmoc'};

tpm_dir = ['../../tpms/' num2str(nChannels) 'ch_diff_perSong_2t2_GP7Left/'];
phi_dir = ['../' num2str(nChannels) 'ch_diff_perSong_2t2_GP7Left/'];

out_file = [num2str(nChannels) 'ch_diff_perSong_phi3_2t2_GP7Left.mat'];

%% Load data file (for determining number of songs)

song_data = load('../../data/data_all.mat');

%% Load network orders

load([tpm_dir 'network_orders.mat']); % gives variable 'part_networks'

%% Load and join

song_max = 6;

pair_songs = zeros(length(ps), length(cond_names));
phis = zeros(length(sets), length(ps), length(cond_names), length(taus), song_max);
state_counters = zeros(nStates, length(sets), length(ps), length(cond_names), length(taus), song_max);
mips = cell(nStates, length(sets), length(ps), length(cond_names), length(taus), song_max);

for tau = 1 : length(taus)
    for cond = 1 : size(cond_names, 2)
        for pair = 1 : length(ps)
            tic;
            % Determine number of songs
            nSongs = length(song_data.data{1, cond, pair});
            pair_songs(pair, cond) = nSongs;
            
            parfor net_c = 1 : size(part_networks, 1)
                network = net_c;
                
                song_mips = cell(nStates, song_max);
                song_phis = zeros(song_max, 1);
                song_stateCounters = zeros(nStates, song_max);
                
                for song = 1 : nSongs
                    
                    filename = ['interPart_' cond_names{cond} '_p' num2str(pair) 'n' num2str(network) 's' num2str(song) 'tau' num2str(taus(tau)) '.mat'];
                    
                    data_tpm = load([tpm_dir filename]);
                    data_tpm = data_tpm.data;
                    data_phi = load([phi_dir filename]);
                    
                    % Weight by occurrence of states
                    phi_sum = data_phi.phis .* data_tpm.state_counters;
                    phi_weighted = sum(phi_sum, 1) ./ sum(data_tpm.state_counters, 1);
                    
                    song_phis(song) = phi_weighted;
                    song_stateCounters(:, song) = data_tpm.state_counters;
                    %phis(net_c, pair, cond, tau, song) = phi_weighted;
                    %state_counters(:, net_c, pair, cond, tau, song) = data_tpm.state_counters;
                    mips_tmp = [data_phi.sias{:}];
                    cuts = [mips_tmp.cut_subsystem];
                    song_mips(:, song) = {cuts.cut};
                    %mips(:, net_c, pair, cond, tau, song) = {cuts.cut};
                end
                
                phis(net_c, pair, cond, tau, :) = song_phis;
                state_counters(:, net_c, pair, cond, tau, :) = song_stateCounters;
                mips(:, net_c, pair, cond, tau, :) = song_mips;

            end
            
            if nSongs < song_max
                phis(:, pair, cond, tau, end-(song_max-nSongs)+1:end) = nan;
            end
            
            % Average across songs
            %phis(:, pair, cond, tau) = mean(song_phis, 2);
            toc
        end
    end
end

%% Save
sets = part_networks;
save(out_file, 'phis', 'state_counters', 'mips', 'pair_songs', 'sets', 'nChannels', 'cond_names', 'taus');