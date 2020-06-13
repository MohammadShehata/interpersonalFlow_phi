%% Description

%{

2-channel region-based analysis

Note - assumes that the atlas is identical across all participants and
trials

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

%% Sort networks based on some meaningful ordering of regions

% % Load one participant's atlas 
% map = load('../data_raw/170719_EM1_T2_tm_PLAY_gp.mat');
% 
% [regions, region_ids] = sort({map.Atlas.Scouts.Region});
% networks_new = nchoosek([region_ids region_ids+length(region_ids)], nChannels);
% 
% % Reorder networks (all networks together)
% phis = phi_data.phis;
% for network = 1 : size(phis, 1)
%     phis(network, :, :, :) = phi_data.phis(all(networks == sort(networks_new(network, :)), 2), :, :, :);
% end
% 
% % Reorder networks (intra and inter separately)
% net_new_intra = networks_new(intra_inds, :);
% net_new_inter = networks_new(inter_inds, :);
% for network = 1 : size(phis_intra, 1)
%     phis_intra(network, :, :, :) = phis_intra(all(intra_nets == sort(net_new_intra(network, :)), 2), :, :, :);
% end
% for network = 1 : size(phis_inter, 1)
%     phis_inter(network, :, :, :) = phis_inter(all(inter_nets == sort(net_new_inter(network, :)), 2), :, :, :);
% end

%% Sort networks based on re-ordering of regions

% % Grouping order
% scout_regions = [3 8 5 12 10 6 1 4 9 11 13 14 7 2];
% scout_labels = {'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
% 
% networks_new = nchoosek([scout_regions scout_regions+max(scout_regions)], nChannels);
% 
% % Reorder networks (all networks together)
% phis = phi_data.phis;
% for network = 1 : size(phis, 1)
%     phis(network, :, :, :) = phi_data.phis(all(networks == sort(networks_new(network, :)), 2), :, :, :);
% end
% 
% % Reorder networks (intra and inter separately)
% net_new_intra = networks_new(intra_inds, :);
% net_new_inter = networks_new(inter_inds, :);
% for network = 1 : size(phis_intra, 1)
%     phis_intra(network, :, :, :) = phis_intra(all(intra_nets == sort(net_new_intra(network, :)), 2), :, :, :);
% end
% for network = 1 : size(phis_inter, 1)
%     phis_inter(network, :, :, :) = phis_inter(all(inter_nets == sort(net_new_inter(network, :)), 2), :, :, :);
% end

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

%% Plot intra for each condition, tau

figure;

phi_mat_pMean = log(permute(mean(intra_mat, 3), [1 2 4 5 3]));

clim = [min(phi_mat_pMean(:)) max(phi_mat_pMean(:))];

subplot_counter = 1;
for condition = 1 : size(intra_mat, 4)
    for tau = 1 : size(intra_mat, 5)
        subplot(size(intra_mat, 4), size(intra_mat, 5), subplot_counter);
        imagesc(phi_mat_pMean(:, :, condition, tau)', clim); colorbar;
        title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        axis square
        subplot_counter = subplot_counter + 1;
    end
end

% Average across taus
figure;
phi_mat_pMean = log(permute(mean(mean(intra_mat, 5), 3), [1 2 4 3 5]));
clim = [min(phi_mat_pMean(:)) max(phi_mat_pMean(:))];
for cond = 1 : size(intra_mat, 4)
    subplot(1, size(intra_mat, 4), cond);
    imagesc(phi_mat_pMean(:, :, cond)', clim); colorbar;
    title([phi_data.cond_names{cond}]);
end

%% Plot in reference to one condition for each tau

figure;

cpos = colormap('autumn');
cneg = colormap('winter');
cmap = cat(1, cneg(1:32, :), cpos(33:64, :));
colormap(cmap);

dims = size(intra_mat);
phi_mat_condRef = intra_mat(:, :, :, 1, :);
phi_mat_condDiff = zeros(dims(1), dims(2), dims(3), dims(4)-1, dims(5));
for condition = 2 : size(intra_mat, 4)
    phi_mat_condDiff(:, :, :, condition-1, :) = phi_mat_condRef - intra_mat(:, :, :, condition, :);
end

pairs = (1:10);
phi_mat_pMean = permute(mean(phi_mat_condDiff(:, :, pairs, :, :), 3), [1 2 4 5 3]);

biggest = max(abs([min(phi_mat_pMean(:)) max(phi_mat_pMean(:))]));
clim = [-biggest biggest];

subplot_counter = 1;
for condition = 2 : size(intra_mat, 4)
    for tau = 1 : size(intra_mat, 5)
        subplot(size(intra_mat, 4)-1, size(intra_mat, 5), subplot_counter);
        imagesc(phi_mat_pMean(:, :, condition-1, tau)', clim); colorbar;
        title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        axis square
        subplot_counter = subplot_counter + 1;
    end
end

%% Plot in reference to average across conditions

figure;

cpos = colormap('autumn');
cneg = colormap('winter');
cmap = cat(1, cneg(1:32, :), cpos(33:64, :));
%cmap = cat(1, flipud(cneg(1:2:64, :)), cpos(1:2:64, :));

colormap(cmap);

phi_mat_condDiff = zeros(size(intra_mat));
for condition = 1 : size(intra_mat, 4)
    phi_mat_condDiff(:, :, :, condition, :) = intra_mat(:, :, :, condition, :) - intra_ref;
end

pairs = (1:10);
phi_mat_pMean = permute(mean(phi_mat_condDiff(:, :, pairs, :, :), 3), [1 2 4 5 3]);

biggest = max(abs([min(phi_mat_pMean(:)) max(phi_mat_pMean(:))]));
clim = [-biggest biggest];

subplot_counter = 1;
for condition = 1 : size(intra_mat, 4)
    for tau = 1 : size(intra_mat, 5)
        subplot(size(intra_mat, 4), size(intra_mat, 5), subplot_counter);
        imagesc(phi_mat_pMean(:, :, condition, tau)', clim); colorbar;
        title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        axis square
        subplot_counter = subplot_counter + 1;
    end
end

% Average across taus
figure;
colormap(cmap)
phi_mat_pMean = permute(mean(mean(phi_mat_condDiff(:, :, :, :, 1:2), 5), 3), [1 2 4 3 5]);
biggest = max(abs([min(phi_mat_pMean(:)) max(phi_mat_pMean(:))]));
clim = [-biggest biggest];
for cond = 1 : size(intra_mat, 4)
    subplot(1, size(intra_mat, 4), cond);
    imagesc(phi_mat_pMean(:, :, cond)', clim); colorbar;
    title([phi_data.cond_names{cond}]);
    axis square
end

%% Plot inter for each condition, tau

figure;

phi_mat_pMean = log(permute(mean(inter_mat, 3), [1 2 4 5 3]));

clim = [min(phi_mat_pMean(:)) max(phi_mat_pMean(:))];

subplot_counter = 1;
for condition = 1 : size(inter_mat, 4)
    for tau = 1 : size(inter_mat, 5)
        subplot(size(inter_mat, 4), size(inter_mat, 5), subplot_counter);
        imagesc(phi_mat_pMean(:, :, condition, tau)', clim); colorbar;
        title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        axis square
        subplot_counter = subplot_counter + 1;
    end
end

% Average across taus
figure;
phi_mat_pMean = log(permute(mean(mean(inter_mat, 5), 3), [1 2 4 3 5]));
clim = [min(phi_mat_pMean(:)) max(phi_mat_pMean(:))];
for cond = 1 : size(inter_mat, 4)
    subplot(1, size(inter_mat, 4), cond);
    imagesc(phi_mat_pMean(:, :, cond)', clim); colorbar;
    title([phi_data.cond_names{cond}]);
end

%% Plot in reference to one condition for each tau

cpos = colormap('autumn');
cneg = colormap('winter');
cmap = cat(1, cneg(1:32, :), cpos(33:64, :));

figure; colormap(cmap);

dims = size(inter_mat);
phi_mat_condRef = inter_mat(:, :, :, 1, :);
phi_mat_condDiff = zeros(dims(1), dims(2), dims(3), dims(4)-1, dims(5));
for condition = 2 : size(inter_mat, 4)
    phi_mat_condDiff(:, :, :, condition-1, :) = phi_mat_condRef - inter_mat(:, :, :, condition, :);
end

pairs = (1:10);
phi_mat_pMean = permute(mean(phi_mat_condDiff(:, :, pairs, :, :), 3), [1 2 4 5 3]);

biggest = max(abs([min(phi_mat_pMean(:)) max(phi_mat_pMean(:))]));
clim = [-biggest biggest];

subplot_counter = 1;
for condition = 2 : size(inter_mat, 4)
    for tau = 1 : size(inter_mat, 5)
        subplot(size(inter_mat, 4)-1, size(inter_mat, 5), subplot_counter);
        imagesc(phi_mat_pMean(:, :, condition-1, tau)', clim); colorbar;
        title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        axis square
        subplot_counter = subplot_counter + 1;
    end
end

%% Plot in reference to average across conditions

figure;

cpos = colormap('autumn');
cneg = colormap('winter');
cmap = cat(1, cneg(1:32, :), cpos(33:64, :));
%cmap = cat(1, flipud(cneg(1:2:64, :)), cpos(1:2:64, :));

colormap(cmap);

phi_mat_condDiff = zeros(size(inter_mat));

for condition = 1 : size(inter_mat, 4)
    phi_mat_condDiff(:, :, :, condition, :) = inter_mat(:, :, :, condition, :) - inter_ref;
end

pairs = (1:10);
phi_mat_pMean = permute(mean(phi_mat_condDiff(:, :, pairs, :, :), 3), [1 2 4 5 3]);

biggest = max(abs([min(phi_mat_pMean(:)) max(phi_mat_pMean(:))]));
clim = [-biggest biggest];

subplot_counter = 1;
for condition = 1 : size(inter_mat, 4)
    for tau = 1 : size(inter_mat, 5)
        subplot(size(inter_mat, 4), size(inter_mat, 5), subplot_counter);
        imagesc(phi_mat_pMean(:, :, condition, tau)', clim); colorbar;
        title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        axis square
        subplot_counter = subplot_counter + 1;
    end
end

% Average across taus
figure;
%colormap(cmap)
phi_mat_pMean = permute(mean(mean(phi_mat_condDiff(:, :, pairs, :, 1:2), 5), 3), [1 2 4 3 5]);
biggest = max(abs([min(phi_mat_pMean(:)) max(phi_mat_pMean(:))]));
clim = [-biggest biggest];
for cond = 1 : size(intra_mat, 4)
    subplot(1, size(intra_mat, 4), cond);
    imagesc(phi_mat_pMean(:, :, cond)', clim); colorbar;
    title([phi_data.cond_names{cond}]);
    axis square
end

%% Sort channel-pairs by phi, for selection of 4ch sets

taus = (1:2);
pairs = (1:10);
select_cond = 1; % which condition to use when sorting channel pairs

% Use inter results
phi_mat_condDiff = zeros(size(inter_mat));
for condition = 1 : size(inter_mat, 4)
    phi_mat_condDiff(:, :, :, condition, :) = inter_mat(:, :, :, condition, :) - inter_ref;
end

phi_mat_pMean = permute(mean(mean(phi_mat_condDiff(:, :, pairs, :, taus), 5), 3), [1 2 4 3 5]);

% Convert to vector with channel labels
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

% Convert channel labels to old channel labels
net_map_orig = relabel_channels(scout_regions, 0, net_map);

% Sort by magnitude, largest first
[phi_vec_sorted, order] = sort(phi_vec, 'desc');
net_map_sorted = net_map(order, :);
net_map_orig_sorted = net_map_orig(order, :);

% Save

%% LME stats

% Difference from reference
phi_ref = repmat(mean(phi_data.phis_perSong, 3, 'omitnan'), [1 1 size(phi_data.phis_perSong, 3) 1 1]);
phi_data = phi_data.phis_perSong - phi_ref;

% Convert matrix to table
dims = size(phi_data);
idx = cell(length(dims), 1);
table_mat = zeros(numel(phi_data), length(dims) + 1);
for row = 1 : numel(phi_data)
    [idx{:}] = ind2sub(dims, row);
    table_mat(row, :) = [[idx{:}] phi_data(row)];
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

%% Create causation matrix

phi_mat = zeros(max(networks(:)), max(networks(:)), size(phis, 2), size(phis, 3), size(phis, 4));
for tau = 1 : size(phis, 4)
    for condition = 1 : size(phis, 3)
        for pair = 1 : size(phis, 2)
            
            for network = 1 : size(networks, 1)
                phi_mat(networks(network, 1), networks(network, 2), pair, condition, tau) = phis(network, pair, condition, tau);
            end
            
        end
    end
end

%% Plot for each condition, tau

figure;

phi_mat_pMean = permute(mean(phi_mat, 3), [1 2 4 5 3]);

clim = [min(phi_mat_pMean(:)) max(phi_mat_pMean(:))];

subplot_counter = 1;
for condition = 1 : size(phi_mat, 4)
    for tau = 1 : size(phi_mat, 5)
        subplot(size(phi_mat, 4), size(phi_mat, 5), subplot_counter);
        imagesc(phi_mat_pMean(:, :, condition, tau), clim); colorbar;
        title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        axis square
        subplot_counter = subplot_counter + 1;
    end
end

%% Plot in reference to one condition for each tau

figure;

dims = size(phi_mat);
phi_mat_condRef = phi_mat(:, :, :, 1, :);
phi_mat_condDiff = zeros(dims(1), dims(2), dims(3), dims(4)-1, dims(5));
for condition = 2 : size(phi_mat, 4)
    phi_mat_condDiff(:, :, :, condition-1, :) = phi_mat_condRef - phi_mat(:, :, :, condition, :);
end

pairs = (1:10);
phi_mat_pMean = permute(mean(phi_mat_condDiff(:, :, pairs, :, :), 3), [1 2 4 5 3]);

clim = [min(phi_mat_pMean(:)) max(phi_mat_pMean(:))];

subplot_counter = 1;
for condition = 2 : size(phi_mat, 4)
    for tau = 1 : size(phi_mat, 5)
        subplot(size(phi_mat, 4)-1, size(phi_mat, 5), subplot_counter);
        imagesc(phi_mat_pMean(:, :, condition-1, tau), clim); colorbar;
        title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        axis square
        subplot_counter = subplot_counter + 1;
    end
end

%% Correlate between each pair of participants

tmp = reshape(permute(phi_data.phis, [1 3 4 2]), [378*3*4 10]);
