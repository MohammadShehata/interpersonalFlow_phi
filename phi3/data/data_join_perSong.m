%% Description

%{

Joins participant pairs, conditions, and timescales into a single file

%}

%% Constants

max_songs = 6;
nChannels = 2; nStates = nChannels^2;
ps = (1:10);
taus = [20 40 60 80];
networks = nchoosek(28, 2);

cond_names = {'tm', 'tmrv', 'tmoc'};

tpm_dir = ['../../tpms/' num2str(nChannels) 'ch_diff_perSong/'];
phi_dir = ['../' num2str(nChannels) 'ch_diff_perSong/'];

out_file = [num2str(nChannels) 'ch_diff_perSong_phi3.mat'];

%% Load data file (for determining number of songs)

song_data = load('../../data/data_all.mat');

%% Load and join

phis = zeros(networks, length(ps), length(cond_names), length(taus));
phis_perSong = zeros(networks, length(ps), length(cond_names), length(taus), max_songs);
phis_perSong(:) = nan;

for tau = 1 : length(taus)
    for cond = 1 : length(cond_names)
        for pair = 1 : length(ps)
            
            % Determine number of songs
            nSongs = length(song_data.data{1, cond, pair});
            
            song_phis = zeros(networks, nSongs);
            
            for song = 1 : nSongs
                
                filename = [cond_names{cond} '_p' num2str(pair) 's' num2str(song) 'tau' num2str(taus(tau)) '.mat'];
                
                data_tpm = load([tpm_dir filename]);
                data_tpm = data_tpm.data;
                data_phi = load([phi_dir filename]);
                
                % Weight by occurrence of states
                phi_sum = data_phi.phis .* data_tpm.state_counters;
                phi_weighted = sum(phi_sum, 1) ./ sum(data_tpm.state_counters, 1);
                
                song_phis(:, song) = phi_weighted;
                
                phis_perSong(:, pair, cond, tau, song) = phi_weighted;
                
            end
            
            % Average across songs
            phis(:, pair, cond, tau) = mean(song_phis, 2);
            
        end
    end
end

%% Save

networks = nchoosek((1:28), 2);
save(out_file, 'phis', 'phis_perSong', 'networks', 'nChannels', 'cond_names', 'taus');