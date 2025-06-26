function [peak_freqs, win_lengths, window_duration, unified_win_length] = ...
    estimate_window_lengths(env_data, fs, nCycles, strategy)
% 估计 envelope 的主频及推荐滑窗长度
%
% 输入：
%   env_data   - [channels x time] 的 envelope 数据（已Hilbert变换）
%   fs         - 采样率
%   nCycles    - 推荐覆盖的周期数（默认7）
%   strategy   - 统一窗口选择策略（'median'（默认）, 'min', 'percentile25'）
%
% 输出：
%   peak_freqs          - 每个通道 envelope 主频（Hz）
%   win_lengths         - 每个通道滑窗长度（单位：点）
%   unified_win_length  - 推荐统一滑窗长度（点）

    if nargin < 3
        nCycles = 7;
    end
    if nargin < 4
        strategy = 'median';
    end

    [~, nChannels] = size(env_data);
    peak_freqs = zeros(nChannels, 1);
    win_lengths = zeros(nChannels, 1);

    for ch = 1:nChannels
        env = env_data(:, ch);

        % 计算 envelope 的功率谱
        [pxx, f] = pwelch(env, [], [], [], fs);

        % 仅考虑 0.5–40 Hz 范围（避免直流 & 高频噪声）
        freq_mask = (f >= 2) & (f <= 40);
        f_band = f(freq_mask);
        p_band = pxx(freq_mask);

        % 找主频（最大能量点）
        [~, idx_max] = max(p_band);
        peak_freq = f_band(idx_max);
        peak_freqs(ch) = peak_freq;

        % 滑窗长度 = n 个周期 × 周期长度 × 采样率
        period_sec = 1 / peak_freq;
        win_lengths(ch) = round(nCycles * period_sec * fs);
    end

    % 统一滑窗策略
    switch lower(strategy)
        case 'min'
            unified_win_length = max(win_lengths);
        case 'percentile25'
            unified_win_length = round(prctile(win_lengths, 25));
        otherwise  % default: median
            unified_win_length = round(median(win_lengths));
    end
    window_duration = unified_win_length/fs;
    fprintf('Recommended window length：%.2f secs（%d points）\n', ...
        unified_win_length/fs, unified_win_length);
end
