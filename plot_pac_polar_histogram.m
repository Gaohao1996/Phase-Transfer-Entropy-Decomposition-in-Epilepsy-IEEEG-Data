function plot_pac_polar_histogram(lf_phase, hf_envelope, num_bins, use_log, modulate_ch, modulated_ch)
% 在笛卡尔坐标中绘制“极坐标风格”的PAC调制图
%
% 参数：
%   lf_phase      - 低频相位（弧度）
%   hf_envelope   - 高频 envelope 振幅
%   num_bins      - 相位 bin 数（默认 18）
%   use_log       - 是否对 envelope 做 log10 缩放（默认 false）

    if nargin < 3
        num_bins = 18;
    end
    if nargin < 4
        use_log = false;
    end

    % === Step 1: Bin 分布
    phase_edges = linspace(-pi, pi, num_bins+1);
    phase_centers = (phase_edges(1:end-1) + phase_edges(2:end)) / 2;

    amp_means = zeros(1, num_bins);
    for i = 1:num_bins
        idx = lf_phase >= phase_edges(i) & lf_phase < phase_edges(i+1);
        amp_means(i) = mean(hf_envelope(idx));
    end

    % === Step 2: 对数缩放
    if use_log
        amp_plot = log10(amp_means + eps);
        amp_plot = amp_plot - min(amp_plot);  % shift to min=0
    else
        amp_plot = amp_means;
    end

    % === Step 3: 绘图（笛卡尔模拟极坐标）
    figure;
    hold on;
    axis equal;
    axis off;

    % 绘制极坐标网格（手动）
    r_max = max(amp_plot) * 1.1;
    for r = linspace(r_max/3, r_max, 3)
        theta = linspace(-pi, pi, 300);
        [x, y] = pol2cart(theta, r);
        plot(x, y, 'k:', 'LineWidth', 0.5);
    end
    for ang = [0 pi/2 pi 3*pi/2]
        [x, y] = pol2cart([ang ang], [0 r_max]);
        plot(x, y, 'k:', 'LineWidth', 0.5);
    end

    % 极坐标标注
    text(0, r_max+0.05, '\pi/2', 'HorizontalAlignment','center');
    text(0, -r_max-0.05, '3\pi/2', 'HorizontalAlignment','center');
    text(-r_max-0.05, 0, '\pi', 'HorizontalAlignment','center');
    text(r_max+0.05, 0, '0', 'HorizontalAlignment','center');

    % === Step 4: 绘制每个扇形
    for i = 1:num_bins
        theta = linspace(phase_edges(i), phase_edges(i+1), 30);
        r = amp_plot(i) * ones(size(theta));
        [x, y] = pol2cart(theta, r);
        fill([0 x 0], [0 y 0], [0.4 0.4 0.4], 'EdgeColor', 'k', 'LineWidth', 1);
    end

    % === Step 5: 红色 modulation 向量
    vec = mean(hf_envelope .* exp(1i * lf_phase));
    r_vec = abs(vec);
    phi = angle(vec);
    [xv, yv] = pol2cart([0 phi], [0 r_vec]);
    plot(xv, yv, 'r', 'LineWidth', 2);
    title(sprintf('PAC Polar Modulation between channel %d and %d', modulate_ch, modulated_ch))
end
