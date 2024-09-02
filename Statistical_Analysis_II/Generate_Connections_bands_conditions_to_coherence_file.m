% Load the existing Coherence.mat file
load('Coherence.mat');

% Define the bands
Coherence.bands = {'1 Delta Band: < 4Hz', '2 Theta Band: 4 to 7Hz', '3 Alpha Band: 8Hz to 12Hz', '4 Beta Band: 13Hz to 30Hz', '5 Gamma1 Band: 30Hz <'};

% Define the conditions
Coherence.Conditions = {'AZ', 'FTD', 'HC'};

% Define the connections
n = 68; % Number of scouts
Nc = 2346; % Brain scout interconnections

matrix = zeros(n);
count = 1;
for i = 1:n
    for j = 1:i
        matrix(j, i) = count;
        count = count + 1;
    end
end

region_mapping = {
    'bankssts L', 'bankssts R', 'caudalanteriorcingulate L', 'caudalanteriorcingulate R', ...
    'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', ...
    'entorhinal L', 'entorhinal R', 'frontalpole L', 'frontalpole R', ...
    'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', ...
    'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', ...
    'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', ...
    'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', ...
    'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', ...
    'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', ...
    'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', ...
    'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', ...
    'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', ...
    'precentral L', 'precentral R', 'precuneus L', 'precuneus R', ...
    'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', ...
    'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', ...
    'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', ...
    'temporalpole L', 'temporalpole R', 'transversetemporal L', 'transversetemporal R'
};

region_map = containers.Map(region_mapping, 1:length(region_mapping));

coordinates = zeros(Nc, 2);
for m = 1:Nc
    [row, col] = find(matrix == m);
    coordinates(m, :) = [row, col];
end
selfConnections = find((coordinates(:, 1) - coordinates(:, 2)) == 0);
nonselfConnections = findAnotinB(1:Nc, selfConnections);
region_names = cell(size(coordinates, 1), size(coordinates, 2));

for i = 1:size(coordinates, 1)
    for j = 1:size(coordinates, 2)
        region_names{i, j} = region_mapping{coordinates(i, j)};
    end
end
Connections = cell(1, Nc); % The list of all the possible connections between brain regions
for i = 1:size(coordinates, 1)
    Connections{i} = sprintf('%d-(%d) %s - (%d) %s', i, coordinates(i, 1), region_names{i, 1}, coordinates(i, 2), region_names{i, 2});
end

Coherence.Connections = Connections;

% Save the updated Coherence structure back to the file
save('Coherence.mat', 'Coherence');
