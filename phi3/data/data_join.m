%% Description

%{

Joins participant pairs, conditions, and timescales into a single file

%}

%% Constants

nChannels = 2; nStates = nChannels^2;
ps = (1:10);
taus = [20 40 60 80];
networks = nchoosek(28, 2);

cond_names = {'tm', 'tmrv', 'tmoc'};

tpm_dir = ['../../tpms/' num2str(nChannels) 'ch_diff/'];
phi_dir = ['../' num2str(nChannels) 'ch/'];

out_file = [num2str(nChannels) 'ch_diff_phi3.mat'];

%% Load and join

phis = zeros(networks, length(ps), length(cond_names), length(taus));

for tau = 1 : length(taus)
    for cond = 1 : length(cond_names)
        for pair = 1 : length(ps)
            
            filename = [cond_names{cond} '_p' num2str(pair) 'tau' num2str(taus(tau)) '.mat'];
            
            data_tpm = load([tpm_dir filename]);
            data_tpm = data_tpm.data;
            data_phi = load([phi_dir filename]);
            
            % Weight by occurrence of states
            phi_sum = data_phi.phis .* data_tpm.state_counters;
            phi_weighted = sum(phi_sum, 1) ./ sum(data_tpm.state_counters, 1);
            
            phis(:, pair, cond, tau) = phi_weighted;
            
        end
    end
end

%% Save

save(out_file, 'phis', 'nChannels', 'cond_names', 'taus');