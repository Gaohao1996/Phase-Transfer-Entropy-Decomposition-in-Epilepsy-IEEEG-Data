function [AAC_matrix_all, AAC_matrix_connection, AAC_pvals] = compute_AAC(envelope_data, Fs, window_len, step_len, nSurrogates)
% envelope_data: N x Ch_num matrix (Hilbert envelope, already precomputed)
% Fs: sampling rate
% window_len: window length in seconds
% step_len: step size in seconds

[N, Ch_num] = size(envelope_data);
win_samples = round(window_len * Fs);
step_samples = round(step_len * Fs);
n_steps = floor((N - win_samples) / step_samples) + 1;

AAC_matrix_all = zeros(Ch_num, Ch_num, n_steps);
AAC_matrix_connection = zeros(Ch_num, Ch_num, n_steps);
AAC_thresh_all = zeros(Ch_num, Ch_num, n_steps);
AAC_pvals = NaN(Ch_num, Ch_num, n_steps);

surrogate_shift = round(win_samples / 2);
disp(surrogate_shift)
for k = 1:n_steps
    idx_start = (k - 1) * step_samples + 1;
    idx_end = idx_start + win_samples - 1;
    win_data = envelope_data(idx_start:idx_end, :);  % [samples x channels]

    % True AAC
    AAC = corrcoef(win_data);
    AAC(logical(eye(Ch_num))) = 0;
    AAC_matrix_all(:,:,k) = AAC;
    
    % Surrogates
    AAC_surrogates = zeros(Ch_num, Ch_num, nSurrogates);
    for s = 1:nSurrogates
        % 在循环外做一次 per-surrogate 全数据打乱
        surrogate_global = zeros(N,Ch_num);
        for ch = 1:Ch_num
            shift_amt = randi([surrogate_shift, N - surrogate_shift]);
            surrogate_global(:, ch) = circshift(envelope_data(:, ch), shift_amt);
        end
        % 在当前滑窗内提取打乱段
        shuffled_data = surrogate_global(idx_start:idx_end, :);
        AAC_surr = corrcoef(shuffled_data);
        AAC_surrogates(:,:,s) = AAC_surr;
    end

    % Compute empirical p-values
    for i = 1:Ch_num
        for j = i+1:Ch_num
            true_r = AAC(i,j);
            surr_r = squeeze(AAC_surrogates(i,j,:));
            p_val = (sum(surr_r >= true_r) + 1) / (nSurrogates + 1);
            AAC_pvals(i,j,k) = p_val;
            AAC_pvals(j,i,k) = p_val;
        end
    end

    % Optional: diagonal = NaN for FDR
    for d = 1:Ch_num
        AAC_pvals(d,d,k) = NaN;
    end
end

end




    

