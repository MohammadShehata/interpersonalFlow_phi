function [med_data] = binarise_diff(data)
% Binarises data based on gradient
%
% Inputs:
%   data = matrix of continuous data (samples x channels x trials)
%
% Outputs:
%   med_data = matrix of binary data (samples-1 x channels x trials)
%       Note - same number of samples as original

[nTimes nChannels nTrials] = size(data);

% Reshape all trials into one super trial
data_allt = permute(data, [2 1 3]); % channels x time x trials
data_allt = reshape(data_allt, [nChannels nTimes*nTrials]);

% Find median across all trials
data_med = median(data_allt, 2);
data_med = repmat(data_med, [1 nTimes*nTrials]);

% 1 if above median
med_data = data_allt > data_med;

% Reshape back individual trials
med_data = reshape(med_data, [nChannels nTimes nTrials]);
med_data = permute(med_data, [2 1 3]);

end
