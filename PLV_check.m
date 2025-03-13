close all;clear;clc;
Epilepsy_data  = readmatrix('sz1_ict.dat')'; %T: time samples X C: Channel numbers
Fs = 400; 
M = size(Epilepsy_data,1); 
T_end = M/Fs;
T = (1/Fs):T_end;
Ch_num = size(Epilepsy_data,2);%channel num
Channels = 1:Ch_num;
LFO_index = 2;
[HFs,LFs] = data_filter1(Epilepsy_data, Fs, [1/Fs,T_end], Channels,LFO_index); 
eeg_data = hilbert(HFs); 


%time window for PLV calculation (HFO)
window_length = ceil(0.2*Fs); %200ms
t_samples = size(HFs,1);
overlap = 0.5 * window_length;
num_windows = floor((t_samples  - window_length) / (window_length - overlap)); % 时间窗数量

eeg_data_timewindow = zeros(Ch_num,window_length,num_windows);

for i = 1:Ch_num
    for j = 1:num_windows
        start_sample = (j-1)*(window_length - overlap)+1;
        end_sample = start_sample+window_length-1;
       
        % 确保索引在范围内
        if end_sample <= t_samples
            window_sample = start_sample:end_sample;
            eeg_data_timewindow(i,:,j) = eeg_data(window_sample,i);
        end
%         window_sample = start_sample:end_sample;
%         eeg_data_timewindow(i,:,j) = eeg_data(window_sample,i);
    end
end

[nc, ns, nt] = size(eeg_data_timewindow ); 
ndat = eeg_data_timewindow./ abs(eeg_data_timewindow ); 
plv = zeros(nc, nc, nt); 
for t = 1: nt 
    plv(:,:, t) = abs(ndat(:,:, t) * ndat(:,:, t)') / ns; 
end

%%plot the heatmap between channels
t_index = 1;  % 选择要可视化的时间窗口
figure;
imagesc(plv(:, :, t_index));  % 绘制热图
colorbar;
caxis([0 1]);  % 归一化 PLV 取值范围
xlabel('Channel');
ylabel('Channel');
title(['PLV Heatmap at Time Window ', num2str(t_index)]);

reference_channel = 11; % 选择与目标通道计算PLV的参考通道
plv_time_series = squeeze(plv(reference_channel, :, :)); % 取出参考通道与所有通道的PLV随时间变化
figure;
imagesc(1:nt, 1:nc, plv_time_series); % 横轴时间窗口，纵轴通道
colorbar;
caxis([0 1]); % 归一化颜色范围
xlabel('Time Window');
ylabel('Channel');
title(['PLV Heatmap: Channel ', num2str(reference_channel), ' vs All Channels']);


%% Graph Visualization




