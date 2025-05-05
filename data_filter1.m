function [HFs,LFs] = data_filter1(signal,Fs, time_period, channels, LFO_index)
% INPUTS:
%          - signal: multichannel EEG recordings, LxN (number of channels x
%            number of samples).
%          - Fs: sampling frequency
%          - time_period: certain time period for analysis
%          - channels: certain channels for analysis
%          - LFO_index: select which LFO for analysis
% OUTPUTS:
%          - HFs:  orignal signal filtered by high-frequency bandpass
%          filter
%          - LFs:  orignal signal filtered by low-frequency bandpass
%          filter
%% i) Information of signal
start_time = time_period(1);
end_time = time_period(2);
samples = (start_time*Fs):(end_time*Fs);
EEG_orig = signal(samples,channels);


%% iii) Filtering setting
%% Delta Oscillation filter design
%passband and stop band setting for different oscillations
Delta_pass = [1.5 3];  
Delta_stop = [0.5 4];  
Theta_pass = [5 7];  
Theta_stop = [4 8];  
Alpha_pass = [9 12];  
Alpha_stop = [8 13];  
Beta_pass = [15 28];  
Beta_stop = [13 30];  

% % FIR filter design
if  LFO_index ==1
    d_LF = designfilt('bandpassfir', 'FilterOrder', 200, ...
                   'StopbandFrequency1', Delta_stop(1), ...
                   'PassbandFrequency1', Delta_pass(1), ...
                   'PassbandFrequency2', Delta_pass(2), ...
                   'StopbandFrequency2', Delta_stop(2), ...
                   'SampleRate', Fs);

    %% Theta Oscillation filter design
    elseif LFO_index == 2
    d_LF = designfilt('bandpassfir', 'FilterOrder', 300, ...
                   'StopbandFrequency1', Theta_stop(1), ...
                   'PassbandFrequency1', Theta_pass(1), ...
                   'PassbandFrequency2', Theta_pass(2), ...
                   'StopbandFrequency2', Theta_stop(2), ...
                   'SampleRate', Fs);


    %% Alpha Oscillation filter design
    elseif LFO_index == 3
    d_LF = designfilt('bandpassfir', 'FilterOrder', 200, ...
                   'StopbandFrequency1', Alpha_stop(1), ...
                   'PassbandFrequency1', Alpha_pass(1), ...
                   'PassbandFrequency2', Alpha_pass(2), ...
                   'StopbandFrequency2', Alpha_stop(2), ...
                   'SampleRate', Fs);

    %% Beta Oscillation filter design
    else
    d_LF = designfilt('bandpassfir', 'FilterOrder', 400, ...
                   'StopbandFrequency1', Beta_stop(1), ...
                   'PassbandFrequency1', Beta_pass(1), ...
                   'PassbandFrequency2', Beta_pass(2), ...
                   'StopbandFrequency2', Beta_stop(2), ...
                   'SampleRate', Fs);
end
LFs = filtfilt(d_LF, EEG_orig);

%[80 150] freq band
% HFO_pass = [85 145];  
% HFO_stop = [75 155];  
% 
% % FIR bandpass filter for HFO
% d_HF = designfilt('bandpassfir', 'FilterOrder', 300, ...
%                'StopbandFrequency1', HFO_stop(1), ...
%                'PassbandFrequency1', HFO_pass(1), ...
%                'PassbandFrequency2', HFO_pass(2), ...
%                'StopbandFrequency2', HFO_stop(2), ...
%                'SampleRate', Fs);

%high pass for freq > 80Hz
Hp_cutoff = 85;  %
d_HF = designfilt('highpassfir', 'FilterOrder', 150, ...  
               'StopbandFrequency', Hp_cutoff-10, ...  
               'PassbandFrequency', Hp_cutoff, ...
               'SampleRate', Fs);
           
HFs = filtfilt(d_HF, EEG_orig);





