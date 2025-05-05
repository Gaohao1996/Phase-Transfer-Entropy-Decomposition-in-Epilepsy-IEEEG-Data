function sizes = normalize_node_size(values, min_s, max_s)
if nargin < 2, min_s = 5; end
if nargin < 3, max_s = 15; end
values = values(:);
if max(values) == min(values)
    sizes = ones(size(values)) * mean([min_s max_s]);
else
    sizes = (values - min(values)) / (max(values) - min(values));
    sizes = sizes * (max_s - min_s) + min_s;
end
end

