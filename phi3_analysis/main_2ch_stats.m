%% Description

%{

LME example

%}

%% Constants

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

%% Relabel networks with actual labels

% Grouping order
scout_regions_order = [3 8 5 12 10 6 1 4 9 11 13 14 7 2];
scout_labels = {'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
[~, scout_regions] = sort(scout_regions_order);

networks = nchoosek([scout_regions scout_regions+max(scout_regions)], nChannels);
intra_nets = networks(intra_inds, :);
inter_nets = networks(inter_inds, :);

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

%% LME table

% Difference from reference
phi_ref = repmat(mean(phi_data.phis_perSong, 3, 'omitnan'), [1 1 size(phi_data.phis_perSong, 3) 1 1]);
phi_vals = phi_data.phis_perSong - phi_ref;

% Convert matrix to table
dims = size(phi_vals);
idx = cell(length(dims), 1);
table_mat = zeros(numel(phi_vals), length(dims) + 1);
for row = 1 : numel(phi_vals)
    [idx{:}] = ind2sub(dims, row);
    table_mat(row, :) = [[idx{:}] phi_vals(row)];
end

% Remove rows with NaN
table_mat(isnan(table_mat(:, end)), :) = [];

% Convert to table
phi_table = array2table(table_mat, 'VariableNames', {'set', 'pair', 'cond', 'tau', 'song', 'phi'});
phi_table.set = nominal(phi_table.set);
phi_table.pair = nominal(phi_table.pair);
phi_table.cond = nominal(phi_table.cond);
phi_table.song = nominal(phi_table.song);

%% LME

model_spec = 'phi ~ cond + tau + cond:tau + (1|pair) + (1|pair:set) + (1|song)';
model_full = fitlme(phi_table, model_spec);

model_null_specs{1} = 'phi ~ cond + tau + (1|pair) + (1|pair:set) + (1|song)';
model_null_specs{2} = 'phi ~ cond + cond:tau + (1|pair) + (1|pair:set) + (1|song)';
model_null_specs{3} = 'phi ~ tau + cond:tau + (1|pair) + (1|pair:set) + (1|song)';
model_nulls = cell(size(model_null_specs));
for null = 1 : length(model_null_specs)
    model_nulls{null} = fitlme(phi_table, model_null_specs{null});
    compare(model_nulls{null}, model_full)
end
