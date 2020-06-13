function [data_down] = downsample_meanShift(data, multiplier)
% Downsamples data by calculating moving averages
% Adds trials corresponding to shifting the windows which are averaged
%
% Inputs:
%   data = matrix (samples x channels x trials)
%   multiplier = integer; size of window to average across
%
% Outputs:
%   data_down = matrix (samples x channels x trials)
%       Note - samples dimension is reduced by factor of multiplier
%       Note - trials dimension is increased by factor of multiplier

sample_window = [0 multiplier-1];
down = movmean(data, sample_window, 1);

% Extract specific windows for each shift
sample_points = (1 : multiplier : size(data, 1));
if mod(sample_points(end), multiplier) > 0
    % truncate (round down to nearest multiple)
    sample_points = sample_points(1:end-1);
end
data_down = zeros(length(sample_points), size(data, 2), size(data, 3)*(multiplier-1)); % number of samples decrease, and number of trials increase
for shift = 0 : multiplier - 1
    data_down(:, :, (1:size(data, 3))+(size(data, 3)*shift)) = down(sample_points+shift, :, :);
end

end

