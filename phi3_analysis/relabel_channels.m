function [mapped] = relabel_channels(map, map_unmap, labels)
% Relabels channel IDs given a label 1:1 map
% OR converts relabelled channels back to their original labels given the
% map
%
% Inputs:
%   map = vector (N x 1); where N is the total number of channels;
%       Each position corresponds to the original channel labels, and the
%       element held in that position gives the new label
%       Assumes channel labels are unique
%   map_unmap = integer (1 or 0);
%       1 = relabel the set of labels using those provided in map
%       0 = convert re-labelled channels back to their original labels
%   labels = matrix; each element holds a channel label
%       (original labels if map_unmap == 1)
%       (relabelled labels of map_unmap == 0)
%
% Outputs:
%   mapped = matrix; same size as labels, with new labels

%{
Expected channel-label map

scout_regions = [3 8 5 12 10 6 1 4 9 11 13 14 7 2];
scout_labels = {'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};

%}

mapped = zeros(size(labels));

if map_unmap == 1
    % Relabel based on map
    for label_c = 1 : numel(labels)
        mapped(label_c) = map(labels(label_c));
    end
else
    % Return to original labelling
    [~, map] = sort(map); % convert map into inverse map
    for label_c = 1 : numel(labels)
        mapped(label_c) = map(labels(label_c));
    end
end

end

