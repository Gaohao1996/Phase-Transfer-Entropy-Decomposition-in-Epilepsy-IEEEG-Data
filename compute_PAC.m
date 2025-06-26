function [PAC_matrix_all, PAC_pvals, phase_win_all, env_win_all, time_slot] = compute_PAC(phase_data, envelope_data, Fs, window_len, step_len, nSurrogates)
% envelope_data: N x Ch_num matrix (Hilbert envelope, already precomputed)
% Fs: sampling rate
% window_len: window length in seconds
% step_len: step size in seconds

[N, Ch_num] = size(envelope_data);
win_samples = round(window_len * Fs);
step_samples = round(step_len * Fs);
n_steps = floor((N - win_samples) / step_samples) + 1;

phase_win = zeros(win_samples,Ch_num);
env_win = zeros(win_samples,Ch_num);
phase_win_all = zeros(win_samples,Ch_num,n_steps);
env_win_all = zeros(win_samples,Ch_num,n_steps);
time_slot = zeros(win_samples,n_steps);

PAC = zeros(Ch_num,Ch_num);
PAC_matrix_all = zeros(Ch_num, Ch_num, n_steps);
% PAC_matrix_connection = zeros(Ch_num, Ch_num, n_steps);
% PAC_thresh_all = zeros(Ch_num, Ch_num, n_steps);
PAC_pvals = NaN(Ch_num, Ch_num, n_steps);

surrogate_shift = round(win_samples / 2);
% disp(surrogate_shift)
for k = 1:n_steps
    % env and phase information in one window
    idx_start = (k - 1) * step_samples + 1;
    idx_end = idx_start + win_samples - 1;
    time_slot(:,k) = (idx_start/Fs):1/Fs:(idx_end/Fs);
    phase_win = phase_data(idx_start:idx_end, :);  % [samples x channels]
    env_win = envelope_data(idx_start:idx_end, :);
    phase_win_all(:,:,k) = phase_win;
    env_win_all(:,:,k) = env_win;
    % True PAC
    for i = 1: Ch_num
        for j = 1:Ch_num
            PAC(i,j) = compute_modulation_index(phase_win(:,i), env_win(:,j));
        end
    end
    PAC(logical(eye(Ch_num))) = 0;
    PAC_matrix_all(:,:,k) = PAC;
    
    % Surrogates
    PAC_surrogates = zeros(Ch_num, Ch_num, nSurrogates);
    % create surrogates
    for s = 1:nSurrogates
        surrogate_amt_global = zeros(N,Ch_num);
        % shuffle whole the envolope data by shiftting a step ([surrogate_shift, N - surrogate_shift])
        for ch = 1:Ch_num
            shift_amt = randi([surrogate_shift, N - surrogate_shift]);
            surrogate_amt_global(:, ch) = circshift(envelope_data(:, ch), shift_amt);
        end
        % extract the envolope data in current time window
        shuffled_amt_win = surrogate_amt_global(idx_start:idx_end, :);
        
        %surrogate PAC
        for i = 1: Ch_num
            for j = 1:Ch_num
                PAC_surrogates(i,j,s) = compute_modulation_index(phase_win(:,i), shuffled_amt_win(:,j));
            end
        end
    end

%     % Percentile threshold for this window
%     PAC_thresh = prctile(PAC_surrogates, 95, 3);  % Ch_num x Ch_num
%     PAC_thresh(logical(eye(Ch_num))) = 0;
%     PAC_thresh_all(:,:,k) = PAC_thresh;
% 
%     % Binary matrix: real AAC > surrogate threshold
%     PAC_binary = double(PAC > PAC_thresh);
%     PAC_matrix_connection(:,:,k) = PAC_binary;

    % Compute empirical p-values
    for i = 1:Ch_num
        for j = i+1:Ch_num
            true_r = PAC(i,j);
            surr_r = squeeze(PAC_surrogates(i,j,:));
            p_val = (sum(surr_r >= true_r) + 1) / (nSurrogates + 1);
            PAC_pvals(i,j,k) = p_val;
            PAC_pvals(j,i,k) = p_val;
        end
    end

    % Optional: diagonal = NaN for FDR
    for d = 1:Ch_num
        PAC_pvals(d,d,k) = NaN;
    end
end

end




    



