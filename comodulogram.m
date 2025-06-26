function [amp_freqs, phase_freqs, MI_matrix]= comodulogram(phase_all, amp_all, phase_freqs, amp_freqs, n_bins)
    % 输入：
    % phase_all: cell array，每个元素是一个低频相位序列（与 amp 对应）
    % amp_all:   cell array，每个元素是一个高频 envelope 序列
    % phase_freqs: vector of phase frequencies (e.g., 4:1:8)
    % amp_freqs:   vector of amplitude frequencies (e.g., 80:10:200)
    % n_bins: number of bins for phase histogram (e.g., 18)

    n_phase = length(phase_freqs);
    n_amp = length(amp_freqs);
    MI_matrix = zeros(n_phase, n_amp);
    
    for i = 1:n_phase
        for j = 1:n_amp
            phase = phase_all{i,j};
            amp = amp_all{i,j};
            MI_matrix(i,j) = compute_MI(phase, amp);
        end
    end
% 
%     % 可视化 comodulogram
%     figure;
%     imagesc(amp_freqs, phase_freqs, MI_matrix);
%     set(gca, 'YDir', 'normal');
%     xlabel('Amplitude Frequency (Hz)');
%     ylabel('Phase Frequency (Hz)');
%     title('PAC Comodulogram (Modulation Index)');
%     colorbar;
end

