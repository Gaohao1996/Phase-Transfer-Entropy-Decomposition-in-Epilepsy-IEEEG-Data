close all;
clear;
clc;
%%  Load data (Time samples x Channel numbers)

% Epilepsy_data  = load('sz1_pre_clean.mat').pre_eeg_clean; 
Epilepsy_data  = load('sz1_ict_clean.mat').ict_eeg_clean; 

Fs = 400; 
M = size(Epilepsy_data,1); 
T_end = M/Fs;
T = (1/Fs):T_end;
Ch_num = size(Epilepsy_data,2);%channel num
Channels = 1:Ch_num;
%select freq band for certain oscillations; 
%  - 1: Delta: 0.5-4 Hz, 4s window
%  - 2: Theta: 4-8 Hz, 1.5s
%  - 3: Alpha: 8-13 Hz, 1s
%  - 4: Beta: 13 - 30 Hz, 0.4s
LFO_index = 2;
[HFs,LFs] = data_filter1(Epilepsy_data, Fs, [1/Fs,T_end], Channels,LFO_index); 
eeg_data = hilbert(LFs); 

%time window for PLV calculation (LFO)
window_length = ceil(1.5*Fs); 
t_samples = size(LFs,1);
overlap = 0.5 * window_length;
num_windows = floor((t_samples  - window_length) / (window_length - overlap))+1; 

eeg_data_timewindow = zeros(Ch_num,window_length,num_windows);

for i = 1:Ch_num
    for j = 1:num_windows
        start_sample = (j-1)*(window_length - overlap)+1;
        end_sample = start_sample+window_length-1;
   
        if end_sample <= t_samples
            window_sample = start_sample:end_sample;
            eeg_data_timewindow(i,:,j) = eeg_data(window_sample,i);
        end
%         window_sample = start_sample:end_sample;
%         eeg_data_timewindow(i,:,j) = eeg_data(window_sample,i);
    end
end

%calculate PLV
plv = PLV(eeg_data_timewindow);

%significance test of PLV
n_perm = 200; 
[plv_zscore, plv_pval] = plv_significance(eeg_data_timewindow, plv, n_perm);
any_nan = sum(isnan(plv_pval(:)));

%%False Discovery Rate(FDR) test
[adj_pval, sig_mask, sig_counts] = fdr_correct_plv_pval(plv_pval, 0.05);

%%weighted adjacent matrix for creating functional epilepsy network
adjacent_weighted = plv.*sig_mask; 

connection_12_11 = zeros(1,num_windows);
for i = 1:num_windows
    connection_12_11(i) = adjacent_weighted(11,12, i);
end
%% plot the heatmap between channels and  Graph structure visualization
%Select a time window for analyses
 t_index = 4;

%% Generate the heatmap plot for PLV
figure;
imagesc(plv(:, :, t_index));  % Heatmap plot
colorbar;
caxis([0 1]);  % Normalization
xlabel('Channel');
ylabel('Channel');
title(['PLV Heatmap in Theta band at Time Window ', num2str(t_index)]);

%% Graph visualization at a certain time window
G = graph(adjacent_weighted(:,:,t_index), 'upper');  %
node_position = load('layout_coords_76ch.mat').layout_coords; % load the node postion of the epilepsy network

%% Statistical analyses of the epilepsy network in seleted time window

%% Degree distribution
deg = degree(G);  
[degCounts, degValues] = groupcounts(deg);

figure;
plot(degValues, degCounts, 'o-');
xlabel('Degree');
ylabel('Number of Nodes');
title(['Degree Distribution at time window ', num2str(t_index)]);
grid on;

%% Degree Bar Chart
nodeIDs = 1:numnodes(G);
figure;
bar(nodeIDs, deg);
xlabel('Node Index');
ylabel('Degree');
title(['Node Degree Bar Chart at time window ', num2str(t_index)]);
grid on;

%% Edge account
edges_change = zeros(1,num_windows);

for t_index = 1: num_windows
    G1 = graph(adjacent_weighted(:,:,t_index), 'upper'); 
    edge_nums = numedges(G1);
    edges_change(t_index) = edge_nums;
end

figure;
plot(1:t_index, edges_change, '-o');
xlabel('Time Window Index');
ylabel('Number of Edges');
title('Edge Count of Theta Band functional network in serizure period');
grid on;

%% Graph visualization
SOZ_nodes =[11,12];
special_nodes = [1,2,3,4,9,10,17,18,19,26,71,72,73,74]; %potential EZ channels

%%Graph setting
n_nodes = size(adjacent_weighted, 1);
node_labels = arrayfun(@num2str, 1:n_nodes, 'UniformOutput', false);
node_strength = sum(G.adjacency > 0, 2); % sum(plv_sig > 0, 2);  % or sum(G.adjacency > 0, 2)
node_sizes = normalize_node_size(node_strength, 2, 8);
node_colors = repmat([0 0 0], n_nodes, 1);  % Regular nodes: Black color
node_colors(SOZ_nodes, :) = repmat([0 1 0], length(SOZ_nodes), 1);  % SOZ channel: Green color
node_colors(special_nodes, :) = repmat([0 0 1], length(special_nodes), 1);  % SOZ channel: Blue color

edge_weights = G.Edges.Weight;
edge_colors = map_strength_to_colormap(edge_weights, autumn);  % colormap
% edge_widths = normalize_edge_width(edge_weights);  

figure;
h = plot(G, ...
    'XData', node_position(:,1), ...
    'YData', node_position(:,2), ...
    'NodeColor', node_colors, ...
    'NodeLabel', node_labels, ...
    'MarkerSize', normalize_node_size(node_strength, 5, 10), ...
    'LineWidth', normalize_edge_width(G.Edges.Weight));


h.EdgeCData = G.Edges.Weight;
h.EdgeColor = 'flat';
colormap(flipud(autumn));  % 

cb = colorbar;
cb.Label.String = 'PLV Strength';

axis ij; 
axis off;
box off;
h.NodeFontSize = 8;
title(['Functional Network in Theta band at Time Window ', num2str(t_index)]);






