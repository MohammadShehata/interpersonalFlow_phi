function [diff_data] = binarise_diff(data)
% Binarises data based on gradient
%
% Inputs:
%   data = matrix of continuous data (samples x channels x trials)
%
% Outputs:
%   diff_data = matrix of binary data (samples-1 x channels x trials)
%       Note - samples is reduced by one (due to differencing)

diff_data = diff(data, 1, 1);

diff_data = diff_data > 0; % 1 if positive gradient

end

