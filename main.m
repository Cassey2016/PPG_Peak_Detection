% =========================================================================
% Below functions are the implementation for the comparison methods in
% paper:
% Han, Dong, Syed K. Bashar, Jesús Lázaro, Fahimeh Mohagheghian, 
% Andrew Peitzsch, Nishat Nishita, Eric Ding, Emily L. Dickson, 
% Danielle DiMezza, Jessica Scott, Cody Whitcomb, Timothy P. Fitzgibbons, 
% David D. McManus, and Ki H. Chon. 2022. 
% "A Real-Time PPG Peak Detection Method for Accurate Determination of 
% Heart Rate during Sinus Rhythm and Cardiac Arrhythmia" 
% Biosensors 12, no. 2: 82. https://doi.org/10.3390/bios12020082 
%
% Please cite our paper if you used our implementation code. Thank you.
% Author: Dong Han (dong.han@uconn.edu), 01/31/2022.
% =========================================================================

% -------------------------------------------------------------------------
% Input: 
% PPG_raw_buffer: should be 30-sec segment.
% fs_PPG_raw: the sampling frequency of the PPG_raw_buffer.
% -------------------------------------------------------------------------
%% Preparation of PPG signal:
addpath('.\func')
[PPG_buffer,fs_PPG] = my_func_prep_PPG_buffer(PPG_raw_buffer,fs_PPG_raw);

%% Method 1: implemented method 1-a
V_max_flag = true; % true == upper peak detection.
addpath('.\method_01_and_02');
output_upper_Shin_2009 = my_peak_compare_Shin_2009(PPG_buffer,fs_PPG,V_max_flag); % Implementation of Shin 2009 paper.

%% Method 2: implemented method 1-b
V_max_flag = false; % false == lower peak detection.
output_lower_Shin_2009 = my_peak_compare_Shin_2009(PPG_buffer,fs_PPG,V_max_flag); % Implementation of Shin 2009 paper.

%% Method 3 & 4: implemented method 2, it has two output peaks in "output_Elgendi_1_2013"
delta = 0.5; % it was 0.1 as mentioned in the paper. But I think 0.5 works better (0.5 is in the billauer's website).
addpath('.\method_03_and_04');
[output_Elgendi_1_2013] = my_Elgendi_2013_method_I_peakdet(PPG_buffer, delta, fs_PPG);

%% Method 5: first derivative and adaptive thresholding method in Li et al. [4] and Elgendi's paper [3]
abpsig = resample(PPG_buffer,fs_abpsig,fs_PPG_buffer); % upsampling it to 125 Hz.
addpath('.\method_05');
[output_Elgendi_2_2013] = my_func_ppg_peakdet_method_05_Elgendi_2013_method_II(abpsig,fs_abpsig);

%% Method 6: implemented method 4
fs_abp = 250; % Hz.
abp = resample(PPG_buffer,fs_abp,fs_PPG); % upsampling it to 125 Hz.
addpath('.\method_06');
[output_Elgendi_3_2013] = my_Elgendi_2013_method_III_peakdet(abp,fs_abp);
            
%% Method 7: event-related moving averages with dynamic threshold method in Elgendi et al.'s paper [3] 
addpath('.\method_07');
[output_Elgendi_4_2013] = my_func_ppg_peakdet_method_07_Elgendi_2013_method_IV(-PPG_raw_buffer,fs_PPG_raw);

%% Method 8 & 9: peak detection on Stationary Wavelet Transform of PPG signal
fs_swt = 125; % Hz.
PPG_swt = resample(PPG_buffer,fs_swt,fs_PPG); % upsampling it to 125 Hz.
addpath('.\method_08_and_09');
[output_Vadrevu_1_2019,output_Vadrevu_2_2019] = my_Vadrevu_2019_peakdet(PPG_swt,fs_swt);