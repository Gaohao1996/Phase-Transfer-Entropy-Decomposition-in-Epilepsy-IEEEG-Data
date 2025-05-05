function [plv_zscore, plv_pval] = plv_significance(data, plv_real, n_perm)
% data: channels x samples x time windows
% fs: sampling rate (e.g., 400)
% win_size: window length in seconds (e.g., 1.5)
% step_size: step size in seconds (e.g., 0.1)
% n_perm: number of permutations (e.g., 200)

n_win = size(data,3);
ch_n = size(plv_real,1);
% plv_real = zeros(ch_n, ch_n, n_win);
plv_zscore = zeros(ch_n, ch_n, n_win);
plv_pval = zeros(ch_n, ch_n, n_win);

for w = 1:n_win
%     idx_start = (w-1)*step_pts + 1;
%     idx_end = idx_start + win_pts - 1;
    data_win = data(:,:,w);
    % 1. 提取相位
%     phase_data = angle(hilbert(data_win')'); % hilbert returns time x channels
%     
%     % 2. 计算真实PLV
%     for i = 1:ch_n-1
%         for j = i+1:ch_n
%             phase_diff = phase_data(i,:) - phase_data(j,:);
%             plv_real(i,j,w) = abs(mean(exp(1i*phase_diff)));
%             plv_real(j,i,w) = plv_real(i,j,w); % 对称
%         end
%     end

    % 3. 构造null分布
    plv_null = zeros(ch_n, ch_n, n_perm);
    for p = 1:n_perm
        shuffled_data = phase_randomize(data_win);
        shuffled_phase = angle(hilbert(shuffled_data')');
        
        for i = 1:ch_n-1
            for j = i+1:ch_n
                phase_diff = shuffled_phase(i,:) - shuffled_phase(j,:);
                plv_null(i,j,p) = abs(mean(exp(1i*phase_diff)));
                plv_null(j,i,p) = plv_null(i,j,p);
            end
        end
    end

    % 4. 计算z-score和p值
    plv_mean = mean(plv_null, 3);
    plv_std = std(plv_null, [], 3);
    plv_zscore(:,:,w) = (plv_real(:,:,w) - plv_mean) ./ plv_std;

    % p值：右尾检验
    for i = 1:ch_n-1
        for j = i+1:ch_n
            real_val = plv_real(i,j,w);
            null_vals = squeeze(plv_null(i,j,:));
            plv_pval(i,j,w) = mean(real_val <= null_vals); % = sum(null_vals >= real_val) / n_perm;
            plv_pval(j,i,w) = plv_pval(i,j,w);
        end
    end
end

end

