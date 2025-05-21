close all; 
clear;
clc;
% % Public epilepsy dataset (SOZ channel:75、76)
pre_eeg= load('sz5_pre_clean.mat').pre_eeg_raw;
ict_eeg  = load('sz5_ict_clean.mat').ict_eeg_raw;
% pre_eeg= load('sz2_pre_clean.mat').pre_eeg;
% ict_eeg  = load('sz2_ict_clean.mat').ict_eeg;
eeg_alltime = [pre_eeg; ict_eeg];


 %% Data Parameter
Fs = 400; %sample rate for this data in this dataset is 512 Hz % The amount of samples in 3 min
M = size(pre_eeg,1); %The amount of samples 
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


% check repeated 0 value in data 
eeg_check =  ict_eeg;
repeat_ratio = zeros(1,Ch_num);

for i = 1:Ch_num
    zero_ratio = sum(eeg_check(:,i) == 0) ./ M;  % x: samples x channels
    repeat_ratio(1,i) = zero_ratio;
    if any(zero_ratio > 0.1)
        warning('zero repeated ratio in this channels is too high(>0.1)');
    end
end

for i = 1:M
    time_zero_ratio = sum(eeg_check(i,:) == 0) ./Ch_num;
    if any(time_zero_ratio > 0.1)
        warning('zero repeated ratio at this time point is too high(>0.1)');
    end
end
% 
% Binary_matrix = ict_eeg ~= 0;
% %imagesc 的heatmap
% figure;
% imagesc(Binary_matrix);
% colormap([1 1 1; 0 0 1]); %
% colorbar;
% title('non-zero elements heatmap');
% xlabel('column');
% ylabel('row');


%% filter the oscillation
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
% N = length(LFs);
% EEG_spec = abs(fft(LFs));  
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
% figure;
% imagesc(T, 1:Ch_num, total_power_HF);
% xlabel('Time (s)');
% ylabel('Channel');
% title('HFO Power Spectrum Change Across Channels');
% colorbar;
% colormap jet;
% axis xy;
% 
% figure;
% imagesc(T1, 1:Ch_num, total_power_LF);
% xlabel('Time (s)');
% ylabel('Channel');
% title('LFO Power Spectrum Change Across Channels');
% colorbar;
% colormap jet;
% axis xy;

%% Adding jitter to the matrix for copula calculation
% ict_eeg = jitter_zero_in_embedding(ict_eeg);

%% extract phase information and nomalisation

%1) Select source, target channels and mediators(conditioned variables)
% Information path test here: 
%     % Drivers          Mediators        Targets
%     1) 75         →        76            →     Grid (1-64)

%     2) 76         →        75            →    Grid

%     3) SOZ ([75 76])      →        FC (65-70)            →     Grid

%     4) SOZ    →       EZ ([1,2,3,4,9,10,11,12,17,18,19,26])       →     Non-EZ (other channels in Grid) 

EZ_channels = [1,2,3,4,9,10,11,12,17,18,19,26]; %potential EZ channels within the grid
Grid_channels = 1:64;
SOZ_channel = [75 76];
deep_channel = [71 72 73 74 75 76]; % Deep hippocamp
stripe_channel = [65 66 67 68 69 70]; %Frontal cortex
source_channel = SOZ_channel;
target_channel = Grid_channels;%setdiff(Grid_channels,EZ_channels); %Grid_channels;
source_target_channels = [source_channel,target_channel];
condition_channel = stripe_channel;%setdiff(EZ_channels,source_target_channels);%deep_channel;

x_channels = [source_channel, target_channel, condition_channel];
sc_num = length(source_channel);
tc_num = length(target_channel);
cc_num = length(condition_channel);
souce_index = 1:sc_num;
target_index = sc_num+1:sc_num+tc_num;
condition_index = sc_num+tc_num+1:length(x_channels);

%% Testing different types of coupling between oscillations
% 1) Raw ieeg data (contains all knids of coupling types)
x = ict_eeg(:,x_channels); %over the whole spectrum
% x_LFs = LFs(:,x_channels); % only the theta band 

% 2) Phase data 
%(Phase-Phase conpling, PPC, which is suitable for low frequency oscillations)
% x = angle(hilbert(LFs)); %Theta Band Signal

%3) Amplitude data
%(Amplitude-Amplitude coupling, AAC, which is suitable for high frequency oscillations)
% x = 

%4) Phase and Amplitude data 
%(Phase-Amplitude coupling, PAC, which is the common coupling way in epilepsy and other activities)




%% TE decomposition 
%Setting the range of two parameters
% model_orders =1:5; % for theta band
model_orders = 1:3; %for raw IEEG data


% delays = 50:100; %one circle for theta oscillation (125:250ms)
delays = 4:4:120; % 10-300ms (for raw IEEG data)

%calculate the best order and delay for TEx→y 
GC_X = copnorm(x(:,souce_index));
GC_Y = copnorm(x(:,target_index));
[best_L, best_delay, ~, ~, ~] = optimize_te_aic_delay(GC_X, GC_Y, [], model_orders, delays);
[drivers_red, drivers_syn, g_red, g_syn]=TE_syn_red(x, souce_index, target_index,best_L,best_delay);

%%surrogates
nsurr=99;
[g_red_surr,g_syn_surr]=CTE_surr(x ,best_L, best_delay, drivers_red,drivers_syn,nsurr,souce_index, target_index,condition_index);

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
figure
plot(g_red,'-*k');hold on
plot(2:cc_num+1,g_red_surr(:,condition_index),'-or');hold on
plot(p_x (sig_idx_red), g_red(sig_idx_red), 'g*', 'MarkerSize', 8, 'LineWidth', 2);
% xlim([0 3])
for i = 1:length(sig_idx_red)
    x_val_red = p_x(sig_idx_red(i));
    y_val_red = g_red(sig_idx_red(i));
    text(x_val_red, y_val_red, sprintf('%d', x_channels(drivers_red(condition_index(sig_idx_red(i)-1)))), ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end
xlabel('Number of Conditioned Variables')
ylabel(sprintf('TE values with order %d and lag %d', best_L, best_delay))
title('Redundant information from SOZ to Grid') 
legend('Raw signal','Surrogate Signal')

%% plot synergistic information change in TE Decomposition
figure
plot(g_syn,'-*k');hold on
plot(2:cc_num+1,g_syn_surr(:,condition_index),'-or');
plot(p_x (sig_idx_syn), g_syn(sig_idx_syn), 'g*', 'MarkerSize', 8, 'LineWidth', 2);
% xlim([0 3])
for i = 1:length(sig_idx_syn)
    x_val_syn = p_x(sig_idx_syn(i));
    y_val_syn = g_syn(sig_idx_syn(i));  
    text(x_val_syn , y_val_syn, sprintf('%d', x_channels(drivers_syn(condition_index(sig_idx_syn(i)-1)))), ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end
xlabel('Number of Conditioned Variables')
ylabel(sprintf('TE values with order %d and lag %d', best_L, best_delay))
title('Synergistic information from SOZ to Grid')
legend('Raw signal','Surrogate Signal')






