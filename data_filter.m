function [HFs,LFs] = data_filter(signal,Fs, time_period, channels, LFO_index, HFO_index)
% INPUTS:
%          - signal: multichannel EEG recordings, LxN (number of channels x
%            number of samples).
%          - Fs: sampling frequency
%          - time_period: certain time period for analysis
%          - channels: certain channels for analysis
%          - LFO_index: select which LFO for analysis
%          - HFO_index: select which LFO for analysis
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


%% Filtering setting
%% LFO filter design
%passband and stop band setting for different oscillations  (% 3 cycles for orders is commonly used ) 
% n_cycles = 3; 

%1.5 cycles for delta or the length of filter will long than samples 
n_cycles1 = 1.5; 
f_low_delta = 0.5;
order_delta = round(n_cycles1 * Fs / f_low_delta);  
% f_low_theta = 4;
% order_theta = round(n_cycles * Fs / f_low_theta);  

%% FIR filter design and padding to avoid edge artificials
if  LFO_index ==1
    d_LF = designfilt('bandpassfir', ...
    'FilterOrder', order_delta, ...
    'CutoffFrequency1', 0.5, ...
    'CutoffFrequency2', 4, ...
    'SampleRate', Fs);
    pad_len = order_delta;
    EEG_pad_delta = [ ...
    flipud(EEG_orig(1:pad_len, :)); ...  %left padding
    EEG_orig; ...   %data
    flipud(EEG_orig(end-pad_len+1:end, :)) ... %right padding
];
    Filtered_EEG_delta = filtfilt(d_LF, EEG_pad_delta);
    LFs = Filtered_EEG_delta(pad_len+1:end-pad_len,:); %cut off the padding data

    %% Theta Oscillation filter design
elseif LFO_index == 2
    d_LF = designfilt('bandpassfir', ...
        'FilterOrder', 600, ...
        'CutoffFrequency1', 4, ...
        'CutoffFrequency2', 8, ...
        'SampleRate', Fs);
    pad_len = 600;
    EEG_pad_theta = [ ...
    flipud(EEG_orig(1:pad_len, :)); ...  %left padding
    EEG_orig; ...   %data
    flipud(EEG_orig(end-pad_len+1:end, :)) ... %right padding
];
    Filtered_EEG_theta = filtfilt(d_LF, EEG_pad_theta);
    LFs = Filtered_EEG_theta(pad_len+1:end-pad_len,:); %cut off the padding data


    %% Alpha Oscillation filter design
elseif LFO_index == 3
    d_LF = designfilt('bandpassfir', ...
        'FilterOrder', 400, ...  
        'CutoffFrequency1', 8, ...
        'CutoffFrequency2', 13, ...
        'SampleRate', Fs); %order = 400 based on testing
    pad_len = 400;
    EEG_pad_alpha = [ ...
    flipud(EEG_orig(1:pad_len, :)); ...  %left padding
    EEG_orig; ...   %data
    flipud(EEG_orig(end-pad_len+1:end, :)) ... %right padding
];
    Filtered_EEG_alpha = filtfilt(d_LF, EEG_pad_alpha);
    LFs = Filtered_EEG_alpha(pad_len+1:end-pad_len,:); %cut off the padding data
        
        %% Beta Oscillation filter design
elseif LFO_index == 4
    d_LF = designfilt('bandpassfir', ...
        'FilterOrder', 400, ...
        'CutoffFrequency1', 13, ...
        'CutoffFrequency2', 30, ...
        'SampleRate', Fs); %order = 400 based on testing
    pad_len = 400;
    EEG_pad_beta = [ ...
    flipud(EEG_orig(1:pad_len, :)); ...  %left padding
    EEG_orig; ...   %data
    flipud(EEG_orig(end-pad_len+1:end, :)) ... %right padding
];
    Filtered_EEG_beta = filtfilt(d_LF, EEG_pad_beta);
    LFs = Filtered_EEG_beta(pad_len+1:end-pad_len,:); %cut off the padding data
end

%% HFO filter design
Gamma_pass = [27 78];  
Gamma_stop = [25 80];
%high pass for freq > 80Hz
Hp_cutoff = 85;  %
if HFO_index == 1
    d_HF = designfilt('bandpassfir', 'FilterOrder', 500, ...
        'StopbandFrequency1', Gamma_stop(1), ...
        'PassbandFrequency1', Gamma_pass(1), ...
        'PassbandFrequency2', Gamma_pass(2), ...
        'StopbandFrequency2', Gamma_stop(2), ...
        'SampleRate', Fs);
    pad_len = 500;
    EEG_pad_gamma = [ ...
    flipud(EEG_orig(1:pad_len, :)); ...  %left padding
    EEG_orig; ...   %data
    flipud(EEG_orig(end-pad_len+1:end, :)) ... %right padding
];
    Filtered_EEG_gamma = filtfilt(d_HF, EEG_pad_gamma);
    HFs = Filtered_EEG_gamma(pad_len+1:end-pad_len,:); %cut off the padding data
    
elseif HFO_index == 2
    d_HF = designfilt('highpassfir', 'FilterOrder', 200, ...
        'StopbandFrequency', Hp_cutoff-10, ...
        'PassbandFrequency', Hp_cutoff, ...
        'SampleRate', Fs); 
    pad_len = 200;
    EEG_pad_HFO = [ ...
    flipud(EEG_orig(1:pad_len, :)); ...  %left padding
    EEG_orig; ...   %data
    flipud(EEG_orig(end-pad_len+1:end, :)) ... %right padding
];
    Filtered_EEG_HFO = filtfilt(d_HF, EEG_pad_HFO);
    HFs = Filtered_EEG_HFO(pad_len+1:end-pad_len,:); %cut off the padding data

end
            
% HFs = filtfilt(d_HF, EEG_orig);

end






