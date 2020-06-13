%% Description

%{

Intra-brain vs inter-brain at each condition, at each timescale

%}

%% Constants

nChannels = 2;

channels_per_brain = 14;

%% Load

source_dir = '../phi3/data/';

source_file = [num2str(nChannels) 'ch_diff_perSong_phi3.mat'];

phi_data = load([source_dir source_file]); % phi_data.phis should be (networks x pairs x conditions x taus)

%% Identify which channel sets are intra-brain and which are inter-brain

networks = nchoosek((1:28), nChannels);

brain_id = networks > channels_per_brain;

intra_inds = logical(diff(brain_id, 1, 2) == 0);
inter_inds = ~intra_inds;

intra_nets = networks(inter_inds, :);
inter_nets = networks(intra_inds, :);

%% Separate intra/inter into different conditions
% But the dimensions are different...

figure;
subplot_counter = 1;
for condition = 1 : size(phi_data.phis, 3)
    for tau = 1 : size(phi_data.phis, 4)
        subplot(size(phi_data.phis, 3), size(phi_data.phis, 4), subplot_counter);
        
        inter = phi_data.phis(inter_inds, :, condition, tau);
        intra = phi_data.phis(intra_inds, :, condition, tau);
                
        histogram(log(intra), 'Normalization', 'pdf', 'LineStyle', 'none'); hold on;
        histogram(log(inter), 'Normalization', 'pdf', 'LineStyle', 'none');
        
        legend('intra', 'inter');
        
        %title(['tau' num2str(phi_data.taus(tau)) ' ' phi_data.cond_names{condition}]);
        title([num2str(mean(log(intra(:)+1)), 2) ' ' num2str(mean(log(inter(:)+1)), 2)]);
        %xlim([0 0.01]);
        
        subplot_counter = subplot_counter + 1;
    end
end

%% Separate intra/inter into different conditions
% cdfplot from 0 to 95% (don't need to show extreme tails

cond_colours = {'r', 'g', 'b'};

phi_means = zeros(size(phi_data.phis, 3), size(phi_data.phis, 4), 2);
phi_stds = zeros(size(phi_means));

figure;
subplot_counter = 1;
for tau = 1 : size(phi_data.phis, 4)
    subplot(1, size(phi_data.phis, 4), subplot_counter);
    hold on;
    for condition = 1 : size(phi_data.phis, 3)
        
        inter = phi_data.phis(inter_inds, :, condition, tau);
        intra = phi_data.phis(intra_inds, :, condition, tau);
        
        h = cdfplot(log(intra(:)+1)); h.LineStyle = '-'; h.Color = cond_colours{condition};
        h = cdfplot(log(inter(:)+1)); h.LineStyle = '--'; h.Color = cond_colours{condition};
        
        % Store for summary plot
        phi_means(condition, tau, 1) = mean(log(intra(:)+1));
        phi_stds(condition, tau, 1) = std(log(intra(:)+1));
        phi_means(condition, tau, 2) = mean(log(inter(:)+1));
        phi_stds(condition, tau, 2) = std(log(inter(:)+1));
        
    end
    
    %ylim([0 0.95]);
    %xlim([0 0.005]);
    
    title(['tau' num2str(phi_data.taus(tau))]);
    
    subplot_counter = subplot_counter + 1;
end

%% Bar plot

tau_offsets = [-0.15 -0.05 0.05 0.15];

figure;
for ch_cond = 1 : size(phi_means, 3)
    subplot(1, size(phi_means, 3), ch_cond);
    hold on;
    for tau = 1 : size(phi_means, 2)
        errorbar((1:length(phi_data.cond_names))+tau_offsets(tau), phi_means(:, tau, ch_cond), phi_stds(:, tau, ch_cond));
    end
end

%% Errorbar plot

ch_cond_names = {'intra', 'inter'};
tau_offsets = [0 5 10 15];

figure;
for ch_cond = 1 : size(phi_means, 3)
    subplot(1, size(phi_means, 3), ch_cond);
    hold on;
    for tau = 1 : size(phi_means, 2)
        errorbar((1:length(phi_data.cond_names))+tau_offsets(tau), phi_means(:, tau, ch_cond), phi_stds(:, tau, ch_cond));
    end
    
    set(gca, 'XTick', (1:length(phi_data.cond_names)), 'XTickLabel', phi_data.cond_names);
    xtickangle(45);
    
    ylabel('log(\Phi+1)');
    
    title(['2ch ' ch_cond_names{ch_cond} ' 20, 40, 60, 80ms']);
end


%% Boxplot

%% Compare ratio intra/inter across conditions

% Average across participants, at each channel pair



%% Hemispheric effects?

