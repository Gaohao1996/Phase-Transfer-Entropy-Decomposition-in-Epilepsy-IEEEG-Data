%% 7. Calculate Modulation Index (MI) 
function MI = compute_modulation_index(lf_phase, hf_amplitude)
  
    num_bins = 18;  
    phase_bins = linspace(-pi, pi, num_bins+1);
    amplitude_means = zeros(1, num_bins);


    for i = 1:num_bins
        idx = (lf_phase >= phase_bins(i)) & (lf_phase < phase_bins(i+1));
        amplitude_means(i) = mean(hf_amplitude(idx));
    end

    amplitude_probs = amplitude_means / sum(amplitude_means);

 
    H = -sum(amplitude_probs .* log(amplitude_probs + eps));  
    Hmax = log(num_bins);
    MI = (Hmax - H) / Hmax;
end
