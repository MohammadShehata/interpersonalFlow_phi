%% Description

%{

Builds a single TPM for each condition, per trial

Builds a TPM for every combination of channels within a participant-pair

Saves TPMs into individual files, one per participant-pair and condition

%}

%% Constants

nChannels = 2;
nValues = 2;

% sampling rate multipliers = 5, 8, 10, 13, 15, 18, 20
% base sampling rate = 4ms; taus = 20, 32, 40, 52, 60, 72, 80
% base sampling rate = 3.9ms; taus = 19.5, 31.2, 39, 50.7, 58.5, 70.2, 78
tau_base = 4;
tau_multipliers = [5 10 15 20];

% Output directory
out_dir = '2ch_diff_perSong/';

cond_names = {'tm', 'tmrv', 'tmoc'};

%% Load

load('../data/data_all.mat');

[~, nConditions, nPairs] = size(data);

%% Build

for pair = 1 : nPairs
    for condition = 1 : nConditions
        pair_data = data(:, condition, pair);
        
        % For storing how many channels per participant
        % so we can identify later which sets of channels are within or
        % across participants
        channels_per_participant = zeros(length(pair_data), 1);
        for part = 1 : length(pair_data)
            channels_per_participant(part) = size(pair_data{part}{1}, 1);
        end
        
        % All possible combinations of channels
        networks = nchoosek((1:sum(channels_per_participant)), nChannels);
        
        for song = 1 : length(pair_data{1})
            
            song_data = cell(size(pair_data));
            
            % Join participants together
            for part = 1 : length(pair_data)
                song_data{part} = pair_data{part}{song};
            end
            song_data = cat(1, song_data{:});

            song_data = permute(song_data, [2 1 3]); % (time x channels x trials)
            
            % Downsample
            parfor multiplier = 1 : length(tau_multipliers) % ~90 seconds per tau for pair 1
                tic;
                tau_multiplier = tau_multipliers(multiplier);
                tau = tau_base * tau_multiplier;
                
                % Downsample
                pair_data_down = downsample_meanShift(song_data, tau_multiplier);
                
                % Binarise
                pair_data_binary = binarise_diff(pair_data_down);
                
                % Build TPM per combination of channels
                tpms = zeros(nChannels^nValues, nChannels, size(networks, 1));
                state_counters = zeros(nChannels^nValues, size(networks, 1));
                for network = 1 : size(networks, 1)
                    
                    [sbs_tpm, state_counters(:, network)] = build_tpm(pair_data_binary(:, networks(network, :), :), 1, nValues)
                    tpms(:, :, network) = tpm_sbs2sbn(sbs_tpm);
                    
                end
                
                % Save
                tpm_struct = struct();
                tpm_struct.tpms = tpms;
                tpm_struct.state_counters = state_counters;
                out_file = [cond_names{condition} '_p' num2str(pair) 's' num2str(song) 'tau' num2str(tau) '.mat'];
                parsave([out_dir out_file], tpm_struct);
                
                toc
                
            end
            
        end
    end
end