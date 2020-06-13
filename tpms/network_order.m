%% Description

%{

Generates random order of networks

%}

%% Constants

nChannels = 4;
channels_per_participant = 14;

% Output directory
out_dir = '4ch_diff_perSong/';

cond_names = {'tm', 'tmrv', 'tmoc'};

%% Load

load('../data/data_all.mat');

[~, nConditions, nPairs] = size(data);

%% All possible within participants

intra_nets = nchoosek((1:channels_per_participant), nChannels);

%% All possible pairings of channel pairs

a = nchoosek((1:channels_per_participant), 2);

inter_nets = zeros(size(a, 1)^2, nChannels);
net_c = 1;
for net_a = 1 : size(a, 1)
    for net_b = 1 : size(a, 1)
        inter_nets(net_c, :) = [a(net_a, :) a(net_b, :)];
        net_c = net_c + 1;
    end
end

%% Generate random order of channel sets per participant pair

network_labels = {'intra', 'inter'};
networks = cell(2, 1);
networks{1} = intra_nets;
networks{2} = inter_nets;

network_orders = cell(size(networks));

for net_type = 1 : length(networks)
    network_orders{net_type} = zeros(size(networks{net_type}, 1), nPairs);
    for pair = 1 : nPairs
        network_orders{net_type}(:, pair) = randperm(size(networks{net_type}, 1));
    end
end

%% Save ordering

save([out_dir 'network_orders.mat'], 'networks', 'network_labels', 'network_orders');