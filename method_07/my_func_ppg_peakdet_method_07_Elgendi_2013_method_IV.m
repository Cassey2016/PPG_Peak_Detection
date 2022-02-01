function output_Elgendi_4_2013 = my_func_ppg_peakdet_method_07_Elgendi_2013_method_IV(raw_PPG,fs_PPG)
% =========================================================================
% This is my implementation of the method IV in this paper:
% Elgendi, Mohamed, et al. 
% "Systolic peak detection in acceleration photoplethysmograms measured from 
% emergency responders in tropical conditions." PLoS One 8.10 (2013): e76585.
%
% Implemented by Dong Han on 03/02/2020.
%
% Please cite our paper if you used this code:
% Han, Dong, Syed K. Bashar, Jesús Lázaro, Fahimeh Mohagheghian, 
% Andrew Peitzsch, Nishat Nishita, Eric Ding, Emily L. Dickson, 
% Danielle DiMezza, Jessica Scott, Cody Whitcomb, Timothy P. Fitzgibbons, 
% David D. McManus, and Ki H. Chon. 2022. 
% "A Real-Time PPG Peak Detection Method for Accurate Determination of 
% Heart Rate during Sinus Rhythm and Cardiac Arrhythmia" 
% Biosensors 12, no. 2: 82. https://doi.org/10.3390/bios12020082 
%
% Please cite our paper if you used our code. Thank you.
% =========================================================================
%% pre-processing - bandpass filtering
    [b, a] = butter(2,[0.5 8]/(fs_PPG/2)); % 2nd order bandpass filter 0.5-8Hz;
    filtered_PPG = filtfilt(b, a, raw_PPG); % zero-phase filter.
    filtered_PPG = filtered_PPG ./ std(filtered_PPG); % normalizing data is very important for my peak detection.
    filtered_PPG = filtered_PPG - mean(filtered_PPG);
    
    debugging_plot_flag = false; % only for plotting debugging figures. 
    
    % clip the signal by keeping the signal above zero.
    % I do not want to do this, so i will move all signal above zero.
    S_n = filtered_PPG;
% ---- Not following the paper to clip signal but move all signal above zero:
%     if min(S_n) < 0
%         Z_n = S_n - min(S_n); % elevate signal above zero.
%     else
%         % the minimum of S_n is still above zero, so do nothing.
%         Z_n = S_n;
%     end
% ---- Following the paper: only keep the positive value:
    Z_n = S_n;
    Z_n(Z_n < 0) = 0; 
%% pre-processing - squaring
    y_n = (Z_n).^2; % element-wise power.
%% feature extraction - generating potential blocks using two moving averages
    W_1 = round(0.111 * fs_PPG); % mentioned as the paper by brute-force search.
    % first moving average:
%     MA_peak = y_n; % for the beginning and ending signal, use the original signal.
%     for nn = 1+round(W_1/2):length(raw_PPG)-round(W_1/2)
%         temp_range = (nn-round(W_1/2)):(nn+round(W_1/2));
%         MA_peak(nn) = sum(y_n(temp_range))/W_1;
%     end
    MA_peak = movmean(y_n,W_1);
    
    % second moving average: 
    W_2 = round(0.667 * fs_PPG);
%     MA_beat = y_n;
%     for nn = 1+round(W_2/2):length(raw_PPG)-round(W_2/2)
%         temp_range = (nn-round(W_2/2)):(nn+round(W_2/2));
%         MA_beat(nn) = sum(y_n(temp_range))/W_2;
%     end
    MA_beat = movmean(y_n,W_2);
%% classification - thresholding
    beta = 0.02; % from the paper, by brute force search.
    z_bar = mean(y_n);
    alpha = beta * z_bar; % offset level.
    THR_1 = MA_beat + alpha;

    Blocks_Of_Interest = zeros(size(MA_peak)); % I initial it as zero.
    for nn = 1:length(MA_peak)
        if MA_peak(nn) > THR_1(nn) % I think it is THR_1(nn).
            Blocks_Of_Interest(nn) = 0.1;
        else
            % since I inital block of interest as zero, so I do not need to
            % assign zero again.
        end
    end
    
    % searh for onset and offset of each block.
    count_blocks = 0;
    block_onset = NaN(size(MA_peak));
    block_offset = NaN(size(MA_peak));
    if any(Blocks_Of_Interest > 0) % there is a block exist.
        for nn = 1:length(MA_peak)
            if nn == 1 && Blocks_Of_Interest(nn) > 0
               % the first point is a block;
               count_blocks = count_blocks + 1; % since the block start from zero, I have to add the counter first.
               block_onset(count_blocks,1) = nn;
            elseif nn == length(MA_peak) && Blocks_Of_Interest(nn) > 0
                % end with a block:
                % no need to add count_blocks;
                block_offset(count_blocks,1) = nn;
            else
                if nn > 1
                    if Blocks_Of_Interest(nn-1) == 0 && Blocks_Of_Interest(nn) > 0 % a jump means a new block.
                        count_blocks = count_blocks + 1;
                        block_onset(count_blocks,1) = nn;
                    elseif Blocks_Of_Interest(nn-1) > 0 && Blocks_Of_Interest(nn) == 0 % a drop means the end of previous block.
                        block_offset(count_blocks,1) = nn;
                    end
                end
            end
        end
    else
        % there is no block existed. Check why.
%        keyboard;
       HR_Elgendi_4_2013 = 0; % there is no peak location.
        S_peaks = 1;
    output_Elgendi_4_2013 = struct('filtered_PPG_Elgendi_4_2013',S_n,...
        'PPG_peak_loc_Elgendi_4_2013',S_peaks,...
        'HR_Elgendi_4_2013',HR_Elgendi_4_2013);
    return
    end
    
    block_onset(isnan(block_onset)) = []; % remove extra elements.
    block_offset(isnan(block_offset)) = []; % remove extra elements.
    if size(block_onset,1) ~= size(block_offset,1)
        % not same number of onset and offset, check here.
        keyboard;
    end
    
    if size(block_onset,1) ~= count_blocks
        keyboard;
    end
    S_peaks = NaN(count_blocks,1);
    THR_2 = W_1;
    
    for jj = 1:count_blocks
        block_idx = [block_onset(jj,1):block_offset(jj,1)];
        [~,I] = max(y_n(block_idx));
        S_peaks(jj,1) = block_onset(jj,1) + I - 1;
    end
    
    if debugging_plot_flag
        figure;
        plot(filtered_PPG);hold on;
        plot(S_peaks,y_n(S_peaks),'r.','markersize',10);
        plot(y_n);
        plot(MA_peak,'k:');
        plot(MA_beat,'r--');
        plot(THR_1,'g.-');
        plot(Blocks_Of_Interest*max(y_n)*10,'color',[0.5,0.5,0.5]); % grey color. I want to make block more obvious.
        
        legend('filtered PPG','peaks', 'squared PPG with clip to zero', 'MA peak', 'MA beat','THR 1', 'Blocks of Interest');
    end
    
    if isempty(S_peaks)
        HR_Elgendi_4_2013 = 0; % there is no peak location.
        S_peaks = 1;
    else
        HR_Elgendi_4_2013 = 60 * fs_PPG ./ diff(S_peaks); % calculate the HR.
    end
    
    output_Elgendi_4_2013 = struct('filtered_PPG_Elgendi_4_2013',S_n,...
        'PPG_peak_loc_Elgendi_4_2013',S_peaks,...
        'HR_Elgendi_4_2013',HR_Elgendi_4_2013);
end