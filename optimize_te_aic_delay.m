function [best_L, best_delay, best_TE, best_AIC, AIC_matrix] = optimize_te_aic_delay(X, Y, Z, window_lengths, delays)
% calculate the model order and delay for Transfer Entropy(TE) calculation
%
% input：
%   X, Y        - source/target
%   Z           - conditioned variables（[] when need to calculate TE other than conditional ）
%   window_lengths - model order
%   delays         - delay
%
% output：
%   best_L, best_delay - optimal order and delay
%   best_TE            - corresponding TE calue
%   best_AIC           - corresponding AIC value
%   AIC_matrix         -  AIC value

    best_AIC = Inf;
    best_L = NaN;
    best_delay = NaN;
    best_TE = NaN;

    num_L = length(window_lengths);
    num_d = length(delays);
    AIC_matrix = NaN(num_L, num_d);

    % detect if Z is empty
    is_conditional = ~isempty(Z);
    if ~is_conditional
%         disp(' Using unconditional TE (Z is empty)');
        Z = [];   
    end

    for j = 1:num_d
        d = delays(j);

        for i = 1:num_L
            L = window_lengths(i);

            try
                TE = conditional_TE_delay(X, Y, Z, L, d);
            catch
                warning('TE calculation error：L=%d, delay=%d', L, d);
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
