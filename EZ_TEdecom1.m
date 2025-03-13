close all; 
clear;
clc;
% Public epilepsy dataset (SOZ channel:11、12)
pre_eeg =  readmatrix('sz1_pre.dat')';
ict_eeg = readmatrix('sz1_ict.dat')';
eeg_alltime = [pre_eeg; ict_eeg];

 %% seizure data visualization
Fs = 400; %sample rate for this data in this dataset is 512 Hz % The amount of samples in 3 min
M = size(pre_eeg,1); %The amount of samples of the whole period
Ch_num = size(pre_eeg,2);%channel num
Channels = 1:Ch_num;
SOZ_channel = [11 12];
Samples = 1:M;
T_end = 10;
T = (1/Fs):(1/Fs):T_end;
T_all = (1/Fs):(1/Fs):(T_end*2);  %20s

%Time domain visionlization
% select_channel = 11; %
% figure;
% plot(T_all,eeg_alltime(:,select_channel));
% title(sprintf('Time series of the channel: %d', select_channel));

figure;
plot(T_all,eeg_alltime);
title('Time series of the data');
%% filter the oscillation; 
%default freq band for HF is [80 150] from referrence;
%select freq band for certain oscillations; 
%  - 1: Delta: 0.5-4 Hz
%  - 2: Theta: 4-8 Hz
%  - 3: Alpha: 8-13 Hz
%  - 4: Beta: 13 - 30 Hz
LFO_index = 2;
[HFs,LFs] = data_filter1(ict_eeg, Fs, [1/Fs,T_end], Channels,LFO_index); 

% [HFs_all,LFs_all] = data_filter1(eeg_alltime, Fs, [1/Fs,(T_end*2)], Channels,LFO_index); 
%Frequency domain check if the filter works
%FFT
% N = length(HFs);
% EEG_spec = abs(fft(HFs));  
% freq = linspace(0, Fs/2, N/2);  
% 
% % plot freq domain spectrum
% figure;
% plot(freq, EEG_spec(1:N/2));
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('Frequency domain spectrum');


%%  PSD changing rate among channels over the whole time period
% 1. all time period
window_length = 1 * Fs;  %
overlap = 0.5 * window_length; % 
nfft = 2^nextpow2(window_length); % 
num_windows = floor((length(Samples) - window_length) / (window_length - overlap)); % window amounts
total_power_HF = zeros(Ch_num, num_windows); % 
total_power_LF = zeros(Ch_num, num_windows); % 

%% PSD over different channels
signal_HF = zeros(1,M);
signal_LF = zeros(1,M);
for ch = 1:Ch_num
    signal_HF = HFs(:,ch)'; 
    [S,F,T,P] = spectrogram(signal_HF, window_length, overlap, nfft, Fs, 'yaxis');
    total_power_HF(ch, 1:size(P,2)) = mean(10*log10(abs(P)), 1);% dB 
end

for ch = 1:Ch_num
    signal_LF = LFs(:,ch)'; 
    [S1,F1,T1,P1] = spectrogram(signal_LF, window_length, overlap, nfft, Fs, 'yaxis');
    total_power_LF(ch, 1:size(P1,2)) = mean(10*log10(abs(P1)), 1);% dB 
end

% % % Threshold setting for potential SOZ channels
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
disp('Potential seizure onset channels:');
disp(potential_channel);

% plot PSD
figure;
imagesc(T, 1:Ch_num, total_power_HF);
xlabel('Time (s)');
ylabel('Channel');
title('HFO Power Spectrum Change Across Channels (Fixed)');
colorbar;
colormap jet;
axis xy;

figure;
imagesc(T_all, 1:Ch_num, total_power_LF);
xlabel('Time (s)');
ylabel('Channel');
title('LFO Power Spectrum Change Across Channels (Fixed)');
colorbar;
colormap jet;
axis xy;


%% Dectect the Phase-Amplitude Coupling(PAC) within the potential SOZ channels
% 1. if SOZ channels are modulated by other channels with PAC
% Parameter for calculating Modulation Index(MI)
HFs_SOZ = HFs(:,potential_channel); 
hf_signal = zeros(size(HFs_SOZ,1),1);
low_phase = angle(hilbert(LFs)); %extract phase range from low freq oscillation
high_amplitude = abs(hilbert(HFs_SOZ));
low_phase_index =  zeros(size(low_phase,1),1);
high_amplitude_index = zeros(size(high_amplitude,1),1);
n_bins = 18;
MI_comparision = zeros(Ch_num-1,1); %calculate the MI between SOZ channl and all other channels
MI = zeros(length(potential_channel),1);
PAC_channel = zeros(length(potential_channel),1); %detect the channel ID that will have the most PAC influence on SOZ channel


%%significance test
num_surrogates = 1000;  
MI_surrogates = zeros(1, num_surrogates);
Significance = zeros(length(SOZ_channel),1); 

for i = 1:length(potential_channel)
        index = setdiff(Channels,potential_channel(i)); %all other channels except the current potential SOZ channel
        for j = 1:Ch_num-1
            low_phase_index = low_phase(:,index(j));
            high_amplitude_index = high_amplitude(:,i); % the ith SOZ channel 
            MI_comparision(j)= compute_modulation_index(low_phase_index, high_amplitude_index);  % 计算 Modulation Index
        end
        [maxvalue, maxindex] = max(MI_comparision);
        MI(i) = maxvalue;
        PAC_channel(i) = index(maxindex);
        %significance test
        hf_signal = HFs_SOZ(:,i);  
        for k = 1:num_surrogates
            shuffled_hf = hf_signal(randperm(size(hf_signal,1)));  % randomize HFO
            shuffled_hf_amplitude = abs(hilbert(shuffled_hf));
            MI_surrogates(k) = compute_modulation_index(low_phase(:,PAC_channel(i)), shuffled_hf_amplitude);
        end
        %% threshold
        p_value = 0.01;  
        MI_threshold = prctile(MI_surrogates, (1 - p_value) * 100);  
        if MI(i) > MI_threshold
            Significance(i) = 1;  %if MI value is significant to believe; 1:yes, 0: no
        end
end


figure;
scatter(potential_channel, MI, 50, 'b', 'filled'); 
xlabel('channel ID'); 
ylabel('Modulation Index'); 
title('MI between other channels and potential EZ channels'); 
grid on; 

figure;
scatter(potential_channel, Significance, 50, 'R', 'filled'); 
xlabel('channel ID'); 
ylabel('significance(1: PAC  0: No PAC)'); 
title('Significance test for PAC'); 
grid on; 

values = PAC_channel(:); % 
unique_vals = unique(values); % 
counts = histc(values, unique_vals); % 
% 
disp('Potential PAC channel and their Counts:');
disp([unique_vals, counts]);

% 2. if SOZ channels have PAC modulation on other channels 
% LFs_SOZ = LFs(:,potential_channel);
% low_phase = angle(hilbert(LFs_SOZ)); %extract phase range from low freq oscillation
% lf_signal = zeros(size(low_phase,1),1);
% high_amplitude = abs(hilbert(HFs));
% low_phase_index =  zeros(size(low_phase,1),1);
% high_amplitude_index = zeros(size(high_amplitude,1),1);
% n_bins = 18;
% Index = zeros(Ch_num-1,1);
% MI_comparision = zeros(Ch_num-1,1); %calculate the MI for comparsion, and select the maximum
% MI = zeros(length(potential_channel),1);
% PAC_channel = zeros(length(potential_channel),1); %detect the channel ID that will have the most PAC influence on SOZ channel
% 
% % % significance test parameter
% num_surrogates = 1000;  % 
% MI_surrogates = zeros(1, num_surrogates);
% Significance = zeros(length(SOZ_channel),1); 
% 
% for i = 1:length(potential_channel)
%         index = setdiff(Channels,potential_channel(i)); %all other channels except the current SOZ channel
%         for j = 1:Ch_num-1
%             low_phase_index = low_phase(:,i);
%             high_amplitude_index = high_amplitude(:,index(j)); % the ith SOZ channel 
%             MI_comparision(j)= compute_modulation_index(low_phase_index, high_amplitude_index);  % 计算 Modulation Index
%         end
%         [maxvalue, maxindex] = max(MI_comparision);
%         MI(i) = maxvalue;
%         PAC_channel(i) = index(maxindex);
% % %         significance test
%         hf_signal = HFs(:,PAC_channel(i));
%         for k = 1:num_surrogates
%             shuffled_hf = hf_signal(randperm(size(hf_signal,1))); 
%             shuffled_hf_amplitude = abs(hilbert(shuffled_hf));
%             MI_surrogates(k) = compute_modulation_index(low_phase_index, shuffled_hf_amplitude);
%         end
% %         threshold
%         p_value = 0.05;  % 
%         MI_threshold = prctile(MI_surrogates, (1 - p_value) * 100); 
%         if MI(i) > MI_threshold
%             Significance(i) = 1;  %Which means  p_value < 0.05
%         end
% end
% 
% values = PAC_channel(:); % 
% unique_vals = unique(values); % 
% counts = histc(values, unique_vals); % 
% % 
% disp('Potential PAC channel and their Counts:');
% disp([unique_vals, counts]);
% 
% 
% figure;
% scatter(potential_channel, MI, 50, 'b', 'filled'); 
% xlabel('channel ID'); 
% ylabel('Modulation Index'); 
% title('MI between potential EZ channels and other channels'); 
% grid on; 
% 
% figure;
% scatter(potential_channel, Significance, 50, 'R', 'filled'); 
% xlabel('channel ID'); 
% ylabel('significance(1: PAC  0: No PAC)'); 
% title('Significance test for PAC'); 
% grid on; 

%% Causality investigation(TE)
%extract phase information from channels and turn them into Gaussian Coupula
%test the time lag between channels
% max_model_order =10;
% lag_samples = 50;
% source_channel = 11;
% target_channel = 9;
% source_data = hilbert(HFs(1:(end-lag_samples),source_channel));
% target_data = hilbert(HFs((lag_samples+1):end,target_channel));
% source_phase = angle(source_data);
% target_phase = angle(target_data);
% source_phase_CDF = ctransform(source_phase);
% target_phase_CDF = ctransform(target_phase);
% source_phase_GauC = copnorm(target_phase_CDF);
% target_phase_GauC = copnorm(target_phase_CDF);
% 
% %test the best model order
% model_order = 1:max_model_order;
% TE_values = TE(source_phase_GauC, target_phase_GauC, model_order);
% 
% figure;
% plot(TE_values);
% xlabel('Model order')
% ylabel('TE value')
% title(['TE value between Channel', num2str(source_channel), ' and Channel', num2str(target_channel)]);










