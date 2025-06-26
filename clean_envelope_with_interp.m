function [env_clean, outlier_idx, Nan_ratio] = clean_envelope_with_interp(env, threshold)
% 输入：
%   env         - 高频包络 (1D 数组)
%   threshold   - 异常检测阈值（建议值：4~6）

% 默认阈值
if nargin < 2
    threshold = 5;
end

% 1. 计算中位数和MAD
med_env = median(env);
mad_env = median(abs(env - med_env));
z_env = (env - med_env) / (mad_env * 1.4826);  % 转为z-score等效

% 2. 找到异常点（高于阈值）
outlier_idx = find(z_env > threshold);

% 3. 插值修复
env_clean = env;
env_clean(outlier_idx) = NaN;
Nan_ratio = sum(isnan(env_clean)) / length(env_clean);

% 找非异常点索引
valid_idx = find(~isnan(env_clean));
% 使用线性插值填补异常点
env_clean = fillmissing(env_clean, 'linear', 'EndValues', 'extrap');

end

