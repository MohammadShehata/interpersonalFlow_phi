%% Description

%{

Builds a single TPM for each condition, per trial

Builds a TPM for every combination of channels within a participant-pair

Saves TPMs into individual files, one per participant-pair and condition

ONLY for inter-person sets, 2 channels from each person

%}

%%

out_dir = '4ch_diff_perSong_3t1_2chOrder/';

%% Constants for 2ch phi

nChannels = 2;

channels_per_brain = 14;

%% Load

source_dir = '../phi3/data/';

source_file = [num2str(nChannels) 'ch_diff_perSong_phi3.mat'];

phi_data = load([source_dir source_file]); % phi_data.phis should be (networks x pairs x conditions x taus)
%phi_data.phis_perSong(isnan(phi_data.phis_perSong)) = 0;
%phi_data.phis = mean(phi_data.phis_perSong, 5, 'omitnan');

%% Identify which channel sets are intra-brain and which are inter-brain

networks = nchoosek((1:28), nChannels);

brain_id = networks > channels_per_brain;

intra_inds = logical(diff(brain_id, 1, 2) == 0);
inter_inds = ~intra_inds;

intra_nets = networks(intra_inds, :);
inter_nets = networks(inter_inds, :);

% Separate intra from inter
phis_intra = phi_data.phis(intra_inds, :, :, :);
phis_inter = phi_data.phis(inter_inds, :, :, :);

%% Create causation matrices for intra-person, inter-person separately

% Intra-person - average across participants
intra_mat = zeros(channels_per_brain, channels_per_brain, size(phi_data.phis, 2), size(phi_data.phis, 3), size(phi_data.phis, 4));
intra_mat_counter = zeros(size(intra_mat));
for tau = 1 : size(phi_data.phis, 4)
    for condition = 1 : size(phi_data.phis, 3)
        for pair = 1 : size(phi_data.phis, 2)
            for n = 1 : size(intra_nets, 1)
                network = sort(intra_nets(n, :));
                network = mod(network, channels_per_brain);
                network(network == 0) = channels_per_brain;
                intra_mat(network(1), network(2), pair, condition, tau) = intra_mat(network(1), network(2), pair, condition, tau) + phis_intra(n, pair, condition, tau);
                intra_mat_counter(network(1), network(2), pair, condition, tau) = intra_mat_counter(network(1), network(2), pair, condition, tau) + 1;
            end
        end
    end
end
intra_mat = intra_mat ./ intra_mat_counter;
intra_mat(isnan(intra_mat)) = 0;

% Inter-person - average upper and lower triangular parts
inter_mat = zeros(channels_per_brain, channels_per_brain, size(phi_data.phis, 2), size(phi_data.phis, 3), size(phi_data.phis, 4));
inter_mat_counter = zeros(size(inter_mat));
for tau = 1 : size(phi_data.phis, 4)
    for condition = 1 : size(phi_data.phis, 3)
        for pair = 1 : size(phi_data.phis, 2)
            for n = 1 : size(inter_nets, 1)
                network = inter_nets(n, :);
                network = mod(network, channels_per_brain);
                network(network == 0) = channels_per_brain;
                network = sort(network);
                inter_mat(network(1), network(2), pair, condition, tau) = inter_mat(network(1), network(2), pair, condition, tau) + phis_inter(n, pair, condition, tau);
                inter_mat_counter(network(1), network(2), pair, condition, tau) = inter_mat_counter(network(1), network(2), pair, condition, tau) + 1;
            end
        end
    end
end
inter_mat = inter_mat ./ inter_mat_counter;
inter_mat(isnan(inter_mat)) = 0;

%% Create reference (average across conditions)

% Intra-person
intra_ref = mean(intra_mat, 4);

% Inter-person
inter_ref = mean(inter_mat, 4);

%% Selection of 4ch sets - for each channel, find best 3 pairs

taus = (1:2);
pairs = (1:10);
select_cond = 1; % which condition to use when sorting channel pairs

% Use inter results
phi_mat_condDiff = zeros(size(inter_mat));
for condition = 1 : size(inter_mat, 4)
    phi_mat_condDiff(:, :, :, condition, :) = inter_mat(:, :, :, condition, :) - inter_ref;
end

% Average across participant pairs and taus
phi_mat_pMean = permute(mean(mean(phi_mat_condDiff(:, :, pairs, :, taus), 5), 3), [1 2 4 3 5]);

% Convert phi matrix to vector, with channel labels
phi_vec = zeros(nchoosek(channels_per_brain, nChannels)+channels_per_brain, 1); % triangle and diagonal of matrix
net_map = zeros(length(phi_vec), nChannels);
c = 1;
for row = 1 : size(phi_mat_pMean, 1)
    for col = row : size(phi_mat_pMean, 2)
        phi_vec(c) = phi_mat_pMean(row, col, select_cond);
        net_map(c, :) = [row col];
        c = c + 1;
    end
end

% Sort by magnitude - largest first
% note that these are the original channel labels
% Get top 3 "best" channels for each channel
part_networks = cell(4, 2); % networks x participants
part_networks(:, 1) = {[6 5 12], [13], [1 2 6], [11]};
part_networks(:, 2) = {[13], [6 5 12], [11], [1 2 6]};
net_counter = 1;
for channel = 1 : channels_per_brain
    
    % Get networks associated with channel
    channel_net_ids = net_map(:, 1) == channel | net_map(:, 2) == channel;
    channel_phis = phi_vec(channel_net_ids, :);
    channel_nets = net_map(channel_net_ids, :);
    
    % Exclude the channel of interest
    other_channel = channel_nets ~= channel;
    other_channels = zeros(size(channel_nets, 1), 1);
    for net = 1 : size(channel_nets, 1)
        tmp = channel_nets(net, other_channel(net, :));
        if isempty(tmp)
            other_channels(net) = channel;
        else
            other_channels(net) = tmp;
        end
    end
    channel_nets = other_channels;
    
    % Get top three linked channels
    [channel_phis_sorted, order] = sort(channel_phis, 'desc');
    channel_nets_sorted = channel_nets(order, :);
    
    % 1-3
    part_networks{net_counter, 1} = channel;
    part_networks{net_counter, 2} = channel_nets_sorted(1:3);
    
    net_counter = net_counter + 1;
    
    % 3-1
    part_networks{net_counter, 1} = channel_nets_sorted(1:3);
    part_networks{net_counter, 2} = channel;
    
    net_counter = net_counter + 1;
    
end

save([out_dir 'network_orders.mat'], 'part_networks');

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
file_prefix = 'interPart_';

cond_names = {'tm', 'tmrv', 'tmoc'};

%% Load

% Time series data
load('../data/data_all.mat');
[~, nConditions, nPairs] = size(data);

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