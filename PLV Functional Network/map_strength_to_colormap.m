function edge_colors = map_strength_to_colormap(weights, cmap)
% weights: edge weight vector
% cmap: N x 3 colormap (e.g. jet, hot, parula)
n_colors = size(cmap,1);
norm_w = (weights - min(weights)) / (max(weights) - min(weights) + eps);
idx = max(1, round(norm_w * (n_colors - 1)) + 1);
edge_colors = cmap(idx, :);
end

