function [best_L, best_delay, best_TE, best_AIC, AIC_matrix] = optimize_te_aic_delay(X, Y, Z, window_lengths, delays)
% 使用 AIC 优化传递熵或条件传递熵的阶数和延迟
%
% 输入：
%   X, Y        - 源和目标信号
%   Z           - 条件变量（可为空）
%   window_lengths - 嵌入阶数列表
%   delays         - 延迟列表
%
% 输出：
%   best_L, best_delay - 最优阶数与延迟
%   best_TE            - 最优 TE
%   best_AIC           - 对应的 AIC
%   AIC_matrix         - 所有 (L,d) 的 AIC 值

    best_AIC = Inf;
    best_L = NaN;
    best_delay = NaN;
    best_TE = NaN;

    num_L = length(window_lengths);
    num_d = length(delays);
    AIC_matrix = NaN(num_L, num_d);

    % 自动识别是否为条件 TE
    is_conditional = ~isempty(Z);
    if ~is_conditional
%         disp(' Using unconditional TE (Z is empty)');
        Z = [];  % 防止出错
    end

    for j = 1:num_d
        d = delays(j);

        for i = 1:num_L
            L = window_lengths(i);

            % 每次只计算一个 (L,d) 组合下的 TE
            try
                TE = conditional_TE_delay(X, Y, Z, L, d);
            catch
                warning('TE 计算失败：L=%d, delay=%d', L, d);
                continue;
            end

            if isnan(TE) || TE < 0
                continue;
            end

            N_eff = size(Y,1) - L - d;
            if N_eff <= 0
                continue;
            end

            logL = N_eff * TE * log(2);
            k = L * size(X,2);

            [aic, ~] = aicbic(logL, k, N_eff);

            if aic < best_AIC
                best_AIC = aic;
                best_L = L;
                best_delay = d;
                best_TE = TE;
            end

            AIC_matrix(i, j) = aic;
        end
    end
end
