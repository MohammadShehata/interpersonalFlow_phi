%% Description

%{

Intra-brain vs inter-brain at each condition, at each timescale

%}

%% Constants

nChannels = 4;

channels_per_brain = 14;

%% Load

source_dir = '../phi3/data/';

source_file = [num2str(nChannels) 'ch_diff_perSong_phi3_2t2_GP7Left.mat'];

phi_data = load([source_dir source_file]); % phi_data.phis should be (networks x pairs x conditions x taus)

%% Count MIP types
% Assumes that first two channels are from P1, second two are from P2

% Find proportion of MIPs separating the two participants

% Identify MIP type - separates two participants or not
mip_interCut = zeros(size(phi_data.mips)); % 1 if the mip cut is between participants
mip_interCut(:) = nan;
mip_equalCut = zeros(size(mip_interCut)); mip_equalCut(:) = nan; % track which cuts divide the system equally (note - this is a superset of inter-participant cuts)
for state = 1 : size(phi_data.mips, 1)
    for net_c = 1 : size(phi_data.mips, 2)
        for pair = 1 : size(phi_data.mips, 3)
            for cond = 1 : size(phi_data.mips, 4)
                for tau = 1 : size(phi_data.mips, 5)
                    for song = 1 : size(phi_data.mips, 6)
                        
                        if ~isempty(phi_data.mips{state, net_c, pair, cond, tau, song})
                            
                            % Get mip
                            mip = phi_data.mips{state, net_c, pair, cond, tau, song};
                            
                            % 0-indexed to 1-indexed
                            mip.from_nodes = mip.from_nodes + 1;
                            mip.to_nodes = mip.to_nodes + 1;
                            
                            % Identify which channel labels belong to which
                            % participants
                            channel_owners = [ones(size(phi_data.sets{net_c, 1}))'; ones(size(phi_data.sets{net_c, 2}))'+1];
                            
                            from_owners = channel_owners(mip.from_nodes);
                            to_owners = channel_owners(mip.to_nodes);
                            
                            mip_interCut(state, net_c, pair, cond, tau, song) = all(diff(from_owners) == 0) & all(diff(to_owners) == 0);
                            mip_equalCut(state, net_c, pair, cond, tau, song) = length(mip.from_nodes) == length(mip.to_nodes);
                            
                        end
                        
                    end
                end
            end
        end
    end
end

% Count and weight by occurrence of states (all cuts)
interCut_counter = mip_interCut .* phi_data.state_counters;
interCut_prop_perSong = permute(sum(interCut_counter, 1) ./ sum(phi_data.state_counters, 1), [2 3 4 5 6 1]);

% Count and weight by occurrence of states (only for equal cuts)
% interCut_counter = mip_interCut .* phi_data.state_counters;
% equalCut_counter = mip_equalCut .* phi_data.state_counters;
% equalCut_counter(equalCut_counter == 0) = 1; % to deal with division by 0
% interCut_prop_perSong = permute(sum(interCut_counter, 1) ./ sum(equalCut_counter, 1), [2 3 4 5 6 1]);

% Only consider the most frequent state
% [~, idx] = max(phi_data.state_counters, [], 1);
% interCut_counter = permute(mip_interCut(idx), [2 3 4 5 6 1]);
% interCut_prop_perSong = sum(interCut_counter, 5) ./ size(interCut_counter, 5);

%% Average across songs and create reference 'condition'

interCut_prop = mean(interCut_prop_perSong, 5, 'omitnan');

interCut_prop_ref = repmat(mean(interCut_prop, 3), [1 1 size(interCut_prop, 3), 1]);

%% Plot proportions

phi_condDiff = interCut_prop;
clim = [min(phi_condDiff(:)) max(phi_condDiff(:))];

% Average across taus
%phi_condDiff = mean(phi_condDiff, 4);

% Colorplot
figure;
cpos = colormap('autumn');
cneg = colormap('winter');
cmap = cat(1, cneg(1:32, :), cpos(33:64, :));
colormap(cmap);
subplot_c = 1;
for tau = 1 : size(phi_condDiff, 4)
    for cond = 1 : size(phi_condDiff, 3)
        subplot(size(phi_condDiff, 4), size(phi_condDiff, 3), subplot_c);
        
        imagesc(phi_condDiff(:, :, cond, tau), clim);
        colorbar;
        
        ylabel('set');
        xlabel('pair');
        
        title([phi_data.cond_names{cond} ' ' num2str(phi_data.taus(tau)) 'ms']);
        
        subplot_c = subplot_c + 1;
    end
end

% Errorbar plot
figure;
tau_offsets = [-0.1 0.1];
for tau = 1 : size(phi_condDiff, 4)
    tmp = reshape(phi_condDiff(:, :, :, tau), [size(phi_condDiff, 1)*size(phi_condDiff, 2) size(phi_condDiff, 3)]); % collapse sets and participant pairs
    errorbar(...
        (1:size(phi_data.phis, 3))+tau_offsets(tau),...
        mean(tmp, 1),...
        std(tmp, [], 1));
    hold on;
end
set(gca, 'XTick', (1:size(phi_data.phis, 3)), 'XTickLabel', phi_data.cond_names);
xlabel('condition');
ylabel('inter-cut / all-cuts');
legend('20ms', '40ms');
title('Inter-participant MIP proportion');

%% Plot difference from average

%phi_condDiff = phi_data.phis - phi_ref;
phi_condDiff = interCut_prop - interCut_prop_ref;
biggest = max(abs([min(phi_condDiff(:)) max(phi_condDiff(:))]));
clim = [-biggest biggest];

% Average across taus
%phi_condDiff = mean(phi_condDiff, 4);

% Colorplot
figure;
cpos = colormap('autumn');
cneg = colormap('winter');
cmap = cat(1, cneg(1:32, :), cpos(33:64, :));
colormap(cmap);
subplot_c = 1;
for tau = 1 : size(phi_condDiff, 4)
    for cond = 1 : size(phi_condDiff, 3)
        subplot(size(phi_condDiff, 4), size(phi_condDiff, 3), subplot_c);
        
        imagesc(phi_condDiff(:, :, cond, tau), clim);
        colorbar;
        
        ylabel('set');
        xlabel('pair');
        
        title([phi_data.cond_names{cond} ' ' num2str(phi_data.taus(tau)) 'ms']);
        
        subplot_c = subplot_c + 1;
    end
end

% Errorbar plot
figure;
tau_offsets = [-0.1 0.1];
for tau = 1 : size(phi_condDiff, 4)
    tmp = reshape(phi_condDiff(:, :, :, tau), [size(phi_condDiff, 1)*size(phi_condDiff, 2) size(phi_condDiff, 3)]); % collapse sets and participant pairs
    errorbar(...
        (1:size(phi_data.phis, 3))+tau_offsets(tau),...
        mean(tmp, 1),...
        std(tmp, [], 1));
    hold on;
end
set(gca, 'XTick', (1:size(phi_data.phis, 3)), 'XTickLabel', phi_data.cond_names);
xlabel('condition');
ylabel('\Delta inter-cut / all-cuts');
legend('20ms', '40ms');
title('cond-mean(tm,tmrv,tmoc)');

%% Plot individual observations

phi_condDiff = interCut_prop - interCut_prop_ref;

figure;

for part = 1 : size(phi_condDiff, 2)
    for net = 1 : size(phi_condDiff, 1)
        for tau = 1 : size(phi_condDiff, 4)
            plot((1:size(phi_condDiff, 3)), squeeze(phi_condDiff(net, part, :, tau)));
            hold on;
        end
    end
end

%% Regress proportion onto average 2ch phi values

taus = (1:2);

% Load 2ch phi values
phi_data_2ch = load([source_dir '2ch_diff_perSong_phi3.mat']);

% Get relevant pair and average
pair_phis = zeros(size(interCut_prop));
for net_c = 1 : size(phi_data.sets, 1)
    net = sort([phi_data.sets{net_c, 1}' phi_data.sets{net_c, 2}']);
    
    sets_2ch = nchoosek(net, 2);
    
    match = all(ismember(phi_data_2ch.networks, sets_2ch), 2);
    
    pair_phis(net_c, :, :, :) = mean(mean(phi_data_2ch.phis_perSong(match, :, :, taus, :), 5, 'omitnan'), 1);
end

% Plot for all participants
figure;
cond_colours = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
for part = 1 : size(interCut_prop, 2)
    subplot(2, 5, part);
    for cond = 1 : size(interCut_prop, 3)
        a = (pair_phis(:, part, cond, :) - mean(pair_phis(:, part, :, :), 3));
        a = log(mean(pair_phis(:, part, cond, :), 4));
        b = mean(interCut_prop(:, part, cond, :), 4);
        scatter(a(:), b(:), [], cond_colours(cond, :), '.'); hold on;
        
        % Line of best fit
        fit = polyfit(a(:), b(:), 1);
        plot(a(:), polyval(fit, a(:)), 'Color', cond_colours(cond, :), 'LineWidth', 2);
        
    end
    
    %xlim([-7.5  -4]);
    ylim([0 1]);
    
    if part == 1
        xlabel('log(mean(2ch-phi))'); ylabel('inter-cut prop.');
        legend('tm', '', 'tmrv', '', 'tmoc', '', 'Location', 'northoutside');
    end
    
    axis square
    
end

% Plot average across participants
figure;
cond_colours = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
for cond = 1 : size(interCut_prop, 3)
    a = (pair_phis(:, part, cond, :) - mean(pair_phis(:, part, :, :), 3));
    a = log(mean(mean(pair_phis(:, :, cond, :), 4), 2));
    b = mean(mean(interCut_prop(:, :, cond, :), 4), 2);
    scatter(a(:), b(:), [], cond_colours(cond, :), '.'); hold on;
    
    % Line of best fit
    fit = polyfit(a(:), b(:), 1);
    plot(a(:), polyval(fit, a(:)), 'Color', cond_colours(cond, :), 'LineWidth', 2);
    
end

xlabel('log(mean(2ch-phi))'); ylabel('inter-cut prop.');
legend('tm', '', 'tmrv', '', 'tmoc', '', 'Location', 'northwest');

axis square

%% Build phi-table for stats

phi_condDiff_perSong = phis_perSong - phi_ref_perSong;

% Convert matrix to table
dims = size(phi_condDiff_perSong);
idx = cell(length(dims), 1);
table_mat = zeros(numel(phi_condDiff_perSong), length(dims) + 1);
for row = 1 : numel(phi_condDiff_perSong)
    [idx{:}] = ind2sub(dims, row);
    table_mat(row, :) = [[idx{:}] phi_condDiff_perSong(row)];
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