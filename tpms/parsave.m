function [] = parsave(fname, data)
% For saving to files inside a parfor loop
%
% Inputs:
%   fname = string; output file location
%   save_data = data to save, can be struct, for saving multiple variables

save(fname, 'data');

end