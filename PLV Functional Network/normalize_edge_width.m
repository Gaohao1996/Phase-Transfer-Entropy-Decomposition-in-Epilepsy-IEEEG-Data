function widths = normalize_edge_width(weights, min_w, max_w)
if nargin < 2, min_w = 0.5; end
if nargin < 3, max_w = 3; end
weights = weights(:);
if max(weights) == min(weights)
    widths = ones(size(weights)) * mean([min_w max_w]);
else
    widths = (weights - min(weights)) / (max(weights) - min(weights));
    widths = widths * (max_w - min_w) + min_w;
end
end

