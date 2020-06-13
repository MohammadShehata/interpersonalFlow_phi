%% Description

%{

Joins participant pairs, conditions, and timescales into a single file

%}

%% Constants

nChannels = 4; nStates = nChannels^2;
ps = (1:10);
taus = [20 40];
%sets = [373 733 2698 5308];
%sets = [2299 2496 2513 2872 3940 6575 6614 7215 7801 8189];
sets = [373 733 2299 2496 2513 2698 2872 3940 5308 6575 6614 7215 7801 8189];

cond_names = {'tm', 'tmrv', 'tmoc'};

tpm_dir = ['../../tpms/' num2str(nChannels) 'ch_diff_perSong/'];
phi_dir = ['../' num2str(nChannels) 'ch_diff_perSong/'];

out_file = [num2str(nChannels) 'ch_diff_perSong_phi3.mat'];

%% Load data file (for determining number of songs)

song_data = load('../../data/data_all.mat');

%% Load and join

song_max = 6;

phis = zeros(length(sets), length(ps), length(cond_names), length(taus), song_max);

for tau = 1 : length(taus)
    for cond = 1 : length(cond_names)
        for pair = 1 : length(ps)
            
            % Determine number of songs
            nSongs = length(song_data.data{1, cond, pair});
            
            song_phis = zeros(length(sets), nSongs);
            
            for net_c = 1 : length(sets)
                network = sets(net_c);
                
                for song = 1 : nSongs
                    
                    filename = ['interPart_' cond_names{cond} '_p' num2str(pair) 'n' num2str(network) 's' num2str(song) 'tau' num2str(taus(tau)) '.mat'];
                    
                    data_tpm = load([tpm_dir filename]);
                    data_tpm = data_tpm.data;
                    data_phi = load([phi_dir filename]);
                    
                    % Weight by occurrence of states
                    phi_sum = data_phi.phis .* data_tpm.state_counters;
                    phi_weighted = sum(phi_sum, 1) ./ sum(data_tpm.state_counters, 1);
                    
                    phis(net_c, pair, cond, tau, song) = phi_weighted;
                    
                end
                
            end
            
            if nSongs < song_max
                phis(:, pair, cond, tau, end-(song_max-nSongs)+1:end) = nan;
            end
            
            % Average across songs
            %phis(:, pair, cond, tau) = mean(song_phis, 2);
            
        end
    end
end

%% Save

save(out_file, 'phis', 'sets', 'nChannels', 'cond_names', 'taus');