%% Description

%{

Builds a single TPM for each condition, per trial

Builds a TPM for every combination of channels within a participant-pair

Saves TPMs into individual files, one per participant-pair and condition

ONLY for inter-person sets, 2 channels from each person

%}

%% Constants

nChannels = 4;
nValues = 2;
channels_per_participant = 14;

% sampling rate multipliers = 5, 8, 10, 13, 15, 18, 20
% base sampling rate = 4ms; taus = 20, 32, 40, 52, 60, 72, 80
% base sampling rate = 3.9ms; taus = 19.5, 31.2, 39, 50.7, 58.5, 70.2, 78
tau_base = 4;
tau_multipliers = [5 10]; %[5 10 15 20];

% Output directory
out_dir = '4ch_diff_perSong_3t1/';
file_prefix = 'interPart_';

cond_names = {'tm', 'tmrv', 'tmoc'};

%% Load

% Time series data
load('../data/data_all.mat');
[~, nConditions, nPairs] = size(data);

%% Specification of which sets to compute for

% Sets - {6,5,12}+{13}
% Sets - {1,2,6}+{11}
part_networks = cell(4, 2); % networks x participants
part_networks(:, 1) = {[6 5 12], [13], [1 2 6], [11]};
part_networks(:, 2) = {[13], [6 5 12], [11], [1 2 6]};

% TODO - save part_networks in a file

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