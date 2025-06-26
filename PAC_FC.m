close all;
clear;
clc;
%%  Load data (Time samples x Channel numbers)
% Epilepsy_data  = load('sz1_pre_clean.mat').pre_eeg_clean; 
Epilepsy_data  = load('sz4_ict_clean.mat').ict_eeg_raw; 
Fs = 400; 
M = size(Epilepsy_data,1); 
T_end = M/Fs;
% T = (1/Fs):T_end;
T = (1/Fs):(1/Fs):T_end;
Ch_num = size(Epilepsy_data,2);%
Channels = 1:Ch_num;

%select freq band for certain oscillations; 
% Low Frequency band
%  - 1: Delta: 0.5-4 Hz
%  - 2: Theta: 4-8 Hz
%  - 3: Alpha: 8-13 Hz
%  - 4: Beta: 13 - 30 Hz
% High Frequency band
%  - 1: Gamma: 30 - 80 Hz
%  - 2: Ripple: 80 - 200 Hz  

LFO_index = 1;
HFO_index = 2;
[HFs,LFs] = data_filter(Epilepsy_data, Fs, [1/Fs,T_end], Channels, LFO_index, HFO_index); 

% Extract phase and envolope data
phase_data = angle(hilbert(LFs));
env_data_raw = abs(hilbert(HFs));
env_data = zeros(size(env_data_raw));
%Measure the amount of interpolation values for each channel
Nan_ratio_all = zeros(Ch_num,1); 

%% process the epileptogenic spikes (avoid artificial PAC modulation later)
for i = Channels
    [env_interp, ~, Nan_ratio] = clean_envelope_with_interp(env_data_raw(:,i), 5);
    env_data(:,i) = env_interp;
    Nan_ratio_all(i) = Nan_ratio;
    if Nan_ratio>0.15
        warning('Too much interpolation for sipkes in channel %d', i)
    end
end

%Compare the two envlope
% test_channel = 10;
% figure
% plot(T,env_data_raw(:,test_channel),'r')
% hold on 
% plot(T,env_data(:,test_channel),'b')
% legend('raw envolopre', 'envolope')


%% Modulation Index(MI) Calculation (PAC matrix)
%moving windows setting
window_length = 4; %Delta - 6s ;1.5s-Theta; 4s - Alpha 1.2s;  Beta 1s
overlap = 0.5; %overrate = 0.5
step_len = (1-overlap)*window_length;

nSurrogates = 100;    
[PAC_matrix_all, PAC_pvals, phase_win_all, env_win_all, time_slot] = compute_PAC(phase_data, env_data, Fs, window_length, step_len, nSurrogates);
any_nan = sum(isnan(PAC_pvals(:)));
num_windows = size(PAC_matrix_all,3);

%%False Discovery Rate(FDR) test
[adj_pval, sig_mask, sig_counts] = fdr_correct(PAC_pvals, 0.05);
PAC_weighted = PAC_matrix_all .* sig_mask;

% [adj_pval, sig_mask, sig_counts] = fdr_correct(AAC_pvals, 0.05);
% AAC_weighted_signed = AAC_matrix_all .* sig_mask;


%% Graph statistics
%select a time window
t_index =1;
G = graph(PAC_weighted(:,:,t_index), 'upper');  %
node_position = load('layout_coords_76ch.mat').layout_coords; % load the node postion of the epilepsy network

% Statistical analyses of the epilepsy network in seleted time window
deg = degree(G);  
[degCounts, degValues] = groupcounts(deg);

%1) Degree Distribution
figure;
plot(degValues, degCounts, 'o-');
xlabel('Degree');
ylabel('Number of Nodes');
title(['Degree Distribution at time window ', num2str(t_index)]);
grid on;

%2) Degree Bar Chart
nodeIDs = 1:numnodes(G);
figure;
bar(nodeIDs, deg);
xlabel('Node Index');
ylabel('Degree');
title(['Node Degree Bar Chart at time window ', num2str(t_index)]);
grid on;

%% Significant Comodulogram Visualization
figure;
imagesc(Channels, Channels, PAC_weighted(:,:,t_index));
set(gca, 'YDir', 'normal');
xlabel('Channels(Amplitude provided)');
ylabel('Channels (Phase provided)');
title(['PAC Comodulogram (Modulation Index) at time window ', num2str(t_index)])
colorbar;

%3) Edge account
edges_change = zeros(1,num_windows);

for t = 1: num_windows
    G1 = graph(PAC_weighted(:,:,t), 'upper'); 
    edge_nums = numedges(G1);
    edges_change(t) = edge_nums;
end

figure;
plot(1:t, edges_change, '-o');
xlabel('Time Window Index');
ylabel('Number of Edges');
title('Edge Count of PAC modulation over all channels');
grid on;


%% Check the modulation in time domain
Phase_modulate = phase_win_all(:,:,t_index);
Envolope_modulated = env_win_all(:,:,t_index);
modulate_ch = 10;
modulated_ch = 11;
Phase_modulate_ch = Phase_modulate(:,modulate_ch); 
Envolope_modulated_ch = Envolope_modulated(:,modulated_ch);

%  Phase-envelope figure of two channels in one time window

figure
plot(time_slot(:,t_index),Phase_modulate_ch,'k')
hold on;
plot(time_slot(:,t_index),Envolope_modulated_ch,'r')
title(sprintf('PAC modulation in time window %d between channel %d and %d', t_index, modulate_ch, modulated_ch))
legend('Phase','Envolope')
 
% polar histogram
plot_pac_polar_histogram(Phase_modulate_ch, Envolope_modulated_ch, 18, false, modulate_ch, modulated_ch);





