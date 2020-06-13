%% Description

%{

Builds a single TPM for each condition, per trial

Builds a TPM for every combination of channels within a participant-pair

Saves TPMs into individual files, one per participant-pair and condition

ONLY for inter-person sets, 2 channels from each person

%}

%% Constants for 4ch phi

nChannels = 4;
nValues = 2;
channels_per_participant = 14;

% sampling rate multipliers = 5, 8, 10, 13, 15, 18, 20
% base sampling rate = 4ms; taus = 20, 32, 40, 52, 60, 72, 80
% base sampling rate = 3.9ms; taus = 19.5, 31.2, 39, 50.7, 58.5, 70.2, 78
tau_base = 4;
tau_multipliers = [5 10]; %[5 10 15 20];

% Output directory
out_dir = '4ch_diff_perSong_3t1_GP7Left/';

file_prefix = 'interPart_';

cond_names = {'tm', 'tmrv', 'tmoc'};

%% Load and re-order channels according to labelling scheme

% Time series data
tic;
load('../data/data_all.mat');
[~, nConditions, nPairs] = size(data);
toc

% Grouping order
scout_regions_order = [3 8 5 12 10 6 1 4 9 11 13 14 7 2];
scout_labels = {'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
[~, scout_regions] = sort(scout_regions_order);

% Re-order
tic;
for x = numel(data)
    for song = numel(data{x})
        data{x}{song} = data{x}{song}(scout_regions, :, :);
    end
end
toc

%% Create networks (3-1)

channels = [1 2 3 4 5 6 7];
fixed_channels = [7];

part_networks = cell(1, 2); % channel sets x participants
net_c = 1;
for fixed_channel = 1 : length(fixed_channels)
    
    % List of channels excluding the fixed one
    valid_channels = channels(channels ~= fixed_channels(fixed_channel));
    
    % Combinations of valid channels
    channel_selections = nchoosek(valid_channels, 2);
    
    % Combine with fixed channels
    for c = 1 : size(channel_selections, 1)
        
        % 3-1
        part_networks{net_c, 1} = sort([fixed_channels(fixed_channel) channel_selections(c, :)]);
        part_networks{net_c, 2} = fixed_channels(fixed_channel);
        net_c = net_c + 1;
        
        % 1-3
        part_networks{net_c, 2} = sort([fixed_channels(fixed_channel) channel_selections(c, :)]);
        part_networks{net_c, 1} = fixed_channels(fixed_channel);
        net_c = net_c + 1;
        
    end
    
end

save([out_dir 'network_orders.mat'], 'part_networks');

%% Build TPMs

for pair = 1 : nPairs
    for condition = 1 : nConditions
        pair_data = data(:, condition, pair);
        parfor network_c = 1 : size(part_networks, 1)
            network = network_c;
            tic;
            
            for song = 1 : length(pair_data{1})
                
                song_data = cell(size(pair_data));
                
                % Get relevant channels for each participant
                for part = 1 : length(pair_data)
                    song_data{part} = pair_data{part}{song}(part_networks{network, part}, :, :);
                end
                
                % Join participant channels together
                song_data = cat(1, song_data{:});
                song_data = permute(song_data, [2 1 3]); % (time x channels x trials)
                
                % Downsample
                for multiplier = 1 : length(tau_multipliers) % ~90 seconds per tau for pair 1
                    %tic;
                    tau_multiplier = tau_multipliers(multiplier);
                    tau = tau_base * tau_multiplier;
                    
                    % Downsample
                    pair_data_down = downsample_meanShift(song_data, tau_multiplier);
                    
                    % Binarise
                    pair_data_binary = binarise_diff(pair_data_down);
                    
                    % Build TPM per combination of channels
                    tpms = zeros(nChannels^nValues, nChannels, 1);
                    state_counters = zeros(nChannels^nValues, 1);
                    
                    [sbs_tpm, state_counters(:)] = build_tpm(pair_data_binary, 1, nValues);
                    tpms(:, :) = tpm_sbs2sbn(sbs_tpm);
                    
                    % Save
                    tpm_struct = struct();
                    tpm_struct.tpms = tpms;
                    tpm_struct.state_counters = state_counters;
                    tpm_struct.network = part_networks(network_c, :);
                    tpm_struct.network_type = '3t1';
                    out_file = [file_prefix cond_names{condition} '_p' num2str(pair) 'n' num2str(network) 's' num2str(song) 'tau' num2str(tau) '.mat'];
                    parsave([out_dir out_file], tpm_struct);
                    
                    %toc
                    
                end
                toc
            end
        end
    end
end