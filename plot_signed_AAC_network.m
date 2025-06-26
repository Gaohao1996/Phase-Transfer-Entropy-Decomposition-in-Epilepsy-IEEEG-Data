function plot_signed_AAC_network(AAC_matrix, node_position, SOZ_nodes, special_nodes, t_index)
% AAC_matrix: signed weighted matrix (NxN, values in [-1, 1])
% node_position: N x 2 matrix of node XY positions
% SOZ_nodes: indices of SOZ channels (green)
% special_nodes: indices of EZ candidate channels (blue)
% t_index: current time window index for display

% Create graph
G = graph(AAC_matrix,'upper');
n_nodes = size(AAC_matrix, 1);

% ----- Node settings -----
node_labels = arrayfun(@num2str, 1:n_nodes, 'UniformOutput', false);
node_strength = sum(abs(AAC_matrix) > 0, 2);  % degree or strength
node_sizes = normalize_node_size(node_strength, 5, 10);

node_colors = repmat([0 0 0], n_nodes, 1);  % Default: black
node_colors(SOZ_nodes, :) = repmat([0 1 0], length(SOZ_nodes), 1);     % Green
node_colors(special_nodes, :) = repmat([0 0 1], length(special_nodes), 1); % Blue

% ----- Edge settings -----
edge_weights = G.Edges.Weight;

% colormap: red-blue for signed edges
cmap = redbluecmap(256);  % assumes redbluecmap.m is available
norm_weights = (edge_weights + 1) / 2;  % [-1,1] -> [0,1]
color_idx = max(1, round(norm_weights * 255) + 1);  % avoid index 0
edge_colors = cmap(color_idx, :);

% Line width based on strength
line_widths = normalize_edge_width(edge_weights, 1, 3);

% ----- Plotting -----
figure;
h = plot(G, ...
    'XData', node_position(:,1), ...
    'YData', node_position(:,2), ...
    'NodeColor', node_colors, ...
    'NodeLabel', node_labels, ...
    'MarkerSize', node_sizes, ...
    'EdgeCData', edge_weights, ... % 自定义转换函数
    'LineWidth', line_widths);

colormap(redbluecmap(256));
colorbar;
caxis([-1 1]);  

axis ij; axis off; box off;
h.NodeFontSize = 8;
title(['AAC Functional Network at Time Window ', num2str(t_index)]);

end

