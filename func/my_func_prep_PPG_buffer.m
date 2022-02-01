function [PPG_buffer,fs_PPG] = my_func_prep_PPG_buffer(PPG_raw_buffer,fs_PPG)
    % Resample PPG to 50 Hz.
    if fs_PPG ~= 50 % Hz
        PPG_down = resample(PPG_raw_buffer,50,fs_PPG);
        fs_PPG = 50;
    else
        PPG_down = PPG_raw_buffer;
    end

    PPG_buffer = PPG_down(:); % Make sure PPG is column vector
    % Standardizing PPG in sub-function.
    PPG_buffer = my_func_standardizing_PPG(PPG_buffer);

    % Filter signal.
    [b, a] = butter(6,[0.5 20]/(fs_PPG/2)); % Bandpass filter.
    PPG_buffer = filtfilt(b, a, PPG_buffer);
end