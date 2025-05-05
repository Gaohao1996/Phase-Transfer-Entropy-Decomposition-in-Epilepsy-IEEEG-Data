close all; 
clear;
clc;
% % Public epilepsy dataset (SOZ channel:11、12)
pre_eeg= load('sz1_pre_clean.mat').pre_eeg_clean;
ict_eeg  = load('sz1_ict_clean.mat').ict_eeg_clean;
% pre_eeg= load('sz5_pre_clean.mat').pre_eeg_raw;
% ict_eeg  = load('sz5_ict_clean.mat').ict_eeg_raw;
eeg_alltime = [pre_eeg; ict_eeg];

 %% Data Parameter
Fs = 400; %sample rate for this data in this dataset is 512 Hz % The amount of samples in 3 min
M = size(pre_eeg,1); %The amount of samples of the whole period
Ch_num = size(pre_eeg,2);%channel num
Channels = 1:Ch_num;
Samples = 1:M;
t_end = M/Fs;
t = (1/Fs):(1/Fs):t_end;
t_all = (1/Fs):(1/Fs):(t_end*2);  %20s

% Time domain visualization
% select_channel = 50; %
% figure;
% plot(pre_eeg(3601:4000,select_channel));
% xlabel('Time(s)')
% title(sprintf('Time series over the  pre-seizure 10s from channel : %d', select_channel))
% figure;
% plot(ict_eeg(1:400,select_channel));
% xlabel('Samples')
% title(sprintf('Time series samples in seizure period from channel : %d', select_channel))

%% filter the oscillation; 
%default freq band for HF is [80 150] from referrence;
%select freq band for certain oscillations; 
%  - 1: Delta: 0.5-4 Hz
%  - 2: Theta: 4-8 Hz
%  - 3: Alpha: 8-13 Hz
%  - 4: Beta: 13 - 30 Hz

LFO_index = 2;
[HFs,LFs] = data_filter1(ict_eeg, Fs, [1/Fs,t_end], Channels,LFO_index); 
[HFs_pre,LFs_pre] = data_filter1(pre_eeg, Fs, [1/Fs,t_end], Channels,LFO_index); 

%check if the filter works
%FFT
% N = length(HFs);
% EEG_spec = abs(fft(HFs));  
% freq = linspace(0, Fs/2, N/2);  
% % plot freq domain spectrum
% figure;
% plot(freq, EEG_spec(1:N/2));
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('Frequency domain spectrum');
% STFT
window_length = ceil(0.5 * Fs);  %
overlap = 0.5 * window_length; % 
nfft = 2^nextpow2(window_length); % 
num_windows = ceil((length(Samples) - window_length) / (window_length - overlap)); % window amounts
total_power_HF = zeros(Ch_num, num_windows); % 
signal_HF = zeros(1,M);
%HFOs PSD in seizure period
for ch = 1:Ch_num
    signal_HF = HFs(:,ch)'; 
    [S,F,T,P] = spectrogram(signal_HF, window_length, overlap, nfft, Fs, 'yaxis');
    total_power_HF(ch, 1:size(P,2)) = mean(10*log10(abs(P)), 1);% dB 
end
%LFOs PSD in pre-seizure period
window_length1 = ceil(1.5*Fs);  %
overlap1 = 0.5 * window_length1; % 
nfft1 = 2^nextpow2(window_length1); % 
num_windows1 = ceil((length(Samples) - window_length1) / (window_length1 - overlap1)); % window amounts
total_power_LF = zeros(Ch_num, num_windows1); % 
signal_LF = zeros(1,M);
for ch = 1:Ch_num
    signal_LF = LFs_pre(:,ch)'; 
    [S1,F1,T1,P1] = spectrogram(signal_LF, window_length1, overlap1, nfft1, Fs, 'yaxis');
    total_power_LF(ch, 1:size(P1,2)) = mean(10*log10(abs(P1)), 1);% dB 
end

 % Threshold setting for potential SOZ channels
threshold = mean(total_power_HF(:)) + 2* std(total_power_HF(:));
time2exceed_threshold = ones(Ch_num,1);
for i = 1:Ch_num
    for k = 1:size(total_power_HF,2)
            if total_power_HF(i,k)>= threshold
                break;
            end
    end
    time2exceed_threshold(i,1) = k; 
end
[rows, cols] = find(total_power_HF > threshold); 
potential_channel = unique(rows);
disp('Potential EZ channels:');
disp(potential_channel);

% plot PSD
figure;
imagesc(T, 1:Ch_num, total_power_HF);
xlabel('Time (s)');
ylabel('Channel');
title('HFO Power Spectrum Change Across Channels');
colorbar;
colormap jet;
axis xy;

figure;
imagesc(T1, 1:Ch_num, total_power_LF);
xlabel('Time (s)');
ylabel('Channel');
title('LFO Power Spectrum Change Across Channels');
colorbar;
colormap jet;
axis xy;


%% extract phase information and nomalisation
%phase data
Phase_eeg = angle(hilbert(LFs));

%Select source, target channels and mediators(conditioned variables)
% Information path test here: 
%     % Drivers          Mediators        Targets
%     1) 11         →        EZ             →     12 
%     2) 12         →        EZ             →     11
%     3) SOZ      →        EZ             →     non- EZ  
%     4) SOZ      →        EZ             →     Hippocamp  
%     5) SOZ      →        FC             →     Hippocamp
%     6) SOZ      →        Hippocamp     →    FC
%     7) Hippocamp    →    EZ       →    SOZ
%     8) Hippocamp    →    EZ       →    FC
%     9) FC          →    EZ    →    Hippocamp

EZ_channels = [1,2,3,4,9,10,11,12,17,18,19,26,71,72,73,74]; %potential EZ channels
SOZ_channel = [11 12];
deep_channel = [71 72 73 74 75 76]; % Deep hippocamp
stripe_channel = [65 66 67 68 69 70]; %Frontal cortex
source_channel = stripe_channel;
target_channel = deep_channel; %non-EZ nodes
source_target_channels = [source_channel,target_channel];
condition_channel = setdiff(EZ_channels,source_target_channels);%deep_channel;%setdiff(EZ_channels,source_target_channels);
%setdiff(special_nodes,source_target_channels);

%Extract phase information
source_phase = Phase_eeg(:,source_channel);
target_phase = Phase_eeg(:,target_channel);
condition_phase = Phase_eeg(:,condition_channel);


%% TE Decomposition test
% % test the path among all channels
x_channels = [source_channel, target_channel, condition_channel];
sc_num = length(source_channel);
tc_num = length(target_channel);
cc_num = length(condition_channel);
souce_index = 1:sc_num;
target_index = sc_num+1:sc_num+tc_num;
condition_index = sc_num+tc_num+1:length(x_channels);

%TE decomposition 
m = length(x_channels);
x_phase = Phase_eeg(:,x_channels);
%range of the parameter
model_orders =1:5;
delays = 50:100; %one circle for theta oscillation

%calculate the best order and delay for TEx→y 
GC_X = copnorm(x_phase(:,souce_index));
GC_Y = copnorm(x_phase(:,target_index));
[best_L, best_delay, ~, ~, ~] = optimize_te_aic_delay(GC_X, GC_Y, [], model_orders, delays);

[drivers_red, drivers_syn, g_red, g_syn]=TE_syn_red(x_phase, souce_index, target_index,best_L,best_delay);

%%surrogates
nsurr=99;
[g_red_surr,g_syn_surr]=CTE_surr(x_phase ,best_L, best_delay, drivers_red,drivers_syn,nsurr,souce_index, target_index,condition_index);

%%p value test
p_x = 1:cc_num+1;
p_values_red = ones(1, cc_num+1);
p_values_syn = ones(1, cc_num+1);
for t = 2:cc_num+1
    rank_red = sum(g_red_surr(:, condition_index(t-1)) <= g_red(t));  % 
    rank_syn = sum(g_syn_surr(:, condition_index(t-1)) >= g_syn(t));
    p_values_red(t) = (rank_red+1) / (size(g_red_surr, 1) + 1);
    p_values_syn(t) = (rank_syn+1) / (size(g_syn_surr, 1) + 1);
end

alpha = 0.05;
sig_idx_red = find(p_values_red < alpha);
sig_idx_syn = find(p_values_syn < alpha);


%% plot redundant information change in TE Decomposition
figure(1)
plot(g_red,'-*k');hold on
plot(2:cc_num+1,g_red_surr(:,condition_index),'-or');hold on
plot(p_x (sig_idx_red), g_red(sig_idx_red), 'g*', 'MarkerSize', 8, 'LineWidth', 2);
for i = 1:length(sig_idx_red)
    x_val_red = p_x(sig_idx_red(i));
    y_val_red = g_red(sig_idx_red(i));
    text(x_val_red, y_val_red, sprintf('%d', condition_channel(sig_idx_red(i)-1)), ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end
xlabel('Number of Conditioned Variables')
ylabel(sprintf('TE values with order %d and lag %d', best_L, best_delay))
title('Redundant information from FC to Hippocamp') 
%FC
%Hippocamp
legend('Raw signal','Surrogate Signal')

%% plot synergistic information change in TE Decomposition
figure(2)
plot(g_syn,'-*k');hold on
plot(2:cc_num+1,g_syn_surr(:,condition_index),'-or');
plot(p_x (sig_idx_syn), g_syn(sig_idx_syn), 'g*', 'MarkerSize', 8, 'LineWidth', 2);
for i = 1:length(sig_idx_syn)
    x_val_syn = p_x(sig_idx_syn(i));
    y_val_syn = g_syn(sig_idx_syn(i));  
    text(x_val_syn , y_val_syn, sprintf('%d', condition_channel(sig_idx_syn(i)-1)), ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end
xlabel('Number of Conditioned Variables')
ylabel(sprintf('TE values with order %d and lag %d', best_L, best_delay))
title('Synergistic information from FC to Hippocamp')
%FC
%Hippocamp
legend('Raw signal','Surrogate Signal')







