function [output_Elgendi_3_2013] = my_Elgendi_2013_method_III_peakdet(raw_PPG,fs_PPG)
% -------------------------------------------------------------------------
% This peak detection function was mentioned in this paper:
% Elgendi, Mohamed, et al. 
% "Systolic peak detection in acceleration photoplethysmograms measured from 
% emergency responders in tropical conditions." PLoS One 8.10 (2013): e76585.
% 
    [r,ssf,my_avg0,A] = my_revise_run_wabp(raw_PPG,fs_PPG);
% -------------------------------------------------------------------------
    if isempty(r)
        HR_Elgendi_3_2013 = 0; % there is no peak location.
        r = 1;
    else
        HR_Elgendi_3_2013 = 60 * fs_PPG ./ diff(r); % calculate the HR.
    end
    A = [A;0;0;0;]; % add zero
    A(1:6) = A(7); % first six plots are all high amplitude.
    output_Elgendi_3_2013 = struct('PPG_peak_loc_Elgendi_3_2013',r,...
        'HR_Elgendi_3_2013',HR_Elgendi_3_2013,...
        'filtered_PPG_Elgendi_3_2013',A,...
        'thres_Elgendi_3_2013',my_avg0);
end