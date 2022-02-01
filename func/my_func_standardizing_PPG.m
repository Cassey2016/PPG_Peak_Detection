function PPG_buffer = my_func_standardizing_PPG(PPG_buffer)
% Standardizing PPG into zero-mean and uni-variance.
    var_sig_PPG = var(PPG_buffer);
    if var_sig_PPG == 0
        univar_sig_PPG = PPG_buffer;
    else
        univar_sig_PPG = sqrt(1/var_sig_PPG) * PPG_buffer;
    end
    zeromean_sig_PPG = univar_sig_PPG - mean(univar_sig_PPG);
    PPG_buffer = zeromean_sig_PPG; % univariance for PPG 30 sec segment
end