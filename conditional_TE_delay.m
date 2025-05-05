function TE_value = conditional_TE_delay(X, Y, Z, embedding_length, delay)
    if nargin < 5
        delay = 1;
    end
    % Calculate TE or CTE based on the fact that if Z is empty
    is_conditional = ~isempty(Z);
    if ~is_conditional
        Z = [];  % 
%         disp('Z is empty so the result is TEx→y');
    end

    N = size(Y, 1);
    if size(X,1) ~= N
        error('Length of X and Y do not match with each other');
    end
    if is_conditional && size(Z,1) ~= N
        error('length of Z and length of  X/Y do not match with each other');
    end

    L = embedding_length;
    eff_samples = N - L - delay;

    if eff_samples < 1
        error('embedding_length=%d 和 delay=%d lead to the empty of eff_samples（eff_samples=%d）', L, delay, eff_samples);
    end
 
    Y_future = Y(L+delay+1:N, :);

    Y_past = zeros(eff_samples, L * size(Y,2));
    X_past = zeros(eff_samples, L * size(X,2));
    if ~isempty(Z)
            Z_past = zeros(eff_samples, L * size(Z,2));
    else
            Z_past = [];
    end

    for l = 1:L
        Y_past(:, (l-1)*size(Y,2)+1 : l*size(Y,2)) = Y(L+delay - l + 1 : N - l, :);
        X_past(:, (l-1)*size(X,2)+1 : l*size(X,2)) = X(L - l + 1 : N - delay - l, :);
        if is_conditional
                Z_past(:, (l-1)*size(Z,2)+1 : l*size(Z,2)) = Z(L - l + 1 : N - delay - l, :);
       end
    end

    if is_conditional
            condition = [Y_past, Z_past];
    else
            condition = Y_past;
    end
    
    TE_value = gccmi_ccc(Y_future, X_past, condition);
end
