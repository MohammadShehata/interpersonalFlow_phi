function [ tpm, state_counter ] = build_tpm(data, tau, n_values)
% Builds a state-by-state tpm from time-series data
%
% Inputs:
%   data = matrix of discretised data (samples x channels x trials)
%   tau = integer; the lag between current and future states
%   n_values = integer; number of states each channel can take at each time
%       (e.g. 2 for ON and OFF)
%
% Outputs:
%   tpm = matrix (n_values^channels x n_values^channels)

% Determine number of possible system states
n_states = n_values ^ size(data, 2);

tpm = zeros(n_states, n_states);

state_counter = zeros(n_states, 1); % counter for each state

for trial = 1 : size(data, 3)
    
    for sample = 1 : size(data, 1) - tau
        sample_current = data(sample, :, trial);
        sample_future = data(sample+tau, :, trial);
        
        % Identify current state
        state_current = state2loli_index(sample_current);
        
        % Identify future state
        state_future = state2loli_index(sample_future);
        
        % Increment TPM transition by 1
        tpm(state_current, state_future) = tpm(state_current, state_future) + 1;
        
        % Increment transition counter
        state_counter(state_current) = state_counter(state_current) + 1;
    end
end

% Divide elements in TPM by transition counter
% If counter is 0, then transition never occurred - to avoid dividing by 0
for counter = 1 : length(state_counter)
    if state_counter(counter) == 0
        state_counter(counter) = 1;
    end
end
for future_state = 1 : size(tpm, 2)
    tpm(:, future_state) = tpm(:, future_state) ./ state_counter;
end

    function [index] = state2loli_index(state)
        % Function to mimic pyphi.convert.state2loli_index()
        % Converts LOLI bit-index (fortran order) into decimal index
        % Inputs:
        %   state = vector of 1s and 0s
        % Outputs:
        %   index = corresponding loli index (1-indexed, not 0-indexed as in Python)
        
        % TODO: flip bit-order if using C order (as opposed to Fortran
        % order)
        
        % LOLI indexing is fortran order (as opposed to C order)
        index = 0;
        for bit = 1 : length(state)
            index = index + state(bit) * 2^(bit-1);
        end
        index = index + 1;
    end

end
