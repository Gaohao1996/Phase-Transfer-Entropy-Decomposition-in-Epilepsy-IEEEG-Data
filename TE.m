function TE_values = TE(X, Y, window_lengths, biascorrect, demeaned)
% test_TE_cmi 计算不同历史嵌入窗口长度下，从 X 到 Y 的传递熵
%
%   TE_values = test_TE_cmi(X, Y, window_lengths, biascorrect, demeaned)
%
% 输入参数：
%   X             - 源时间序列，样本在行，维度在列
%   Y             - 目标时间序列，样本在行，维度在列
%   window_lengths- 嵌入窗口长度的向量（例如 [1, 2, 3, 4]）
%   biascorrect   - 是否进行偏差校正（可选，默认 true）
%   demeaned      - 数据是否已经去均值（可选，默认 false）
%
% 输出参数：
%   TE_values     - 每个窗口长度对应的传递熵值（单位：bit）

    if nargin < 3
        error('至少需要三个输入参数：X, Y, 和 window_lengths');
    end
    if nargin < 4
        biascorrect = true;
    end
    if nargin < 5
        demeaned = false;
    end

    % 检查 X 和 Y 的样本数是否一致
    N = size(Y, 1);
    if size(X, 1) ~= N
        error('X 和 Y 的样本数不一致');
    end

    num_windows = length(window_lengths);
    TE_values = zeros(num_windows, 1);

    % 针对每个窗口长度计算传递熵
    for idx = 1:num_windows
        L = window_lengths(idx);
        % 有效样本数：从 t = L 到 t = N-1（因为需要 t+1）
        eff_samples = N - L;
        if eff_samples < 1
            error('窗口长度 L=%d 太大，导致有效样本数不足', L);
        end

        % 预分配嵌入矩阵
        % 每一行对应一个样本，列为嵌入向量（将历史时刻拼接为一行）
        Y_future_emb = zeros(eff_samples, size(Y, 2));
        Y_past_emb   = zeros(eff_samples, L * size(Y, 2));
        X_past_emb   = zeros(eff_samples, L * size(X, 2));

        % 构造嵌入：对 t = L, L+1, ..., N-1
        for t = L:(N-1)
            idx_eff = t - L + 1;
            % 未来状态：Y 在 t+1 时刻的值
            Y_future_emb(idx_eff, :) = Y(t+1, :);
            % 嵌入：历史 L 个时刻，顺序为从最旧到最新
            tempY = [];
            tempX = [];
            for k = 0:(L-1)
                tempY = [tempY, Y(t - L + 1 + k, :)];
                tempX = [tempX, X(t - L + 1 + k, :)];
            end
            Y_past_emb(idx_eff, :) = tempY;
            X_past_emb(idx_eff, :) = tempX;
        end

        % 计算传递熵： TE = I(Y_future; X_past | Y_past)
        TE_values(idx) = cmi_ggg(Y_future_emb, X_past_emb, Y_past_emb, biascorrect, demeaned);
    end
end
