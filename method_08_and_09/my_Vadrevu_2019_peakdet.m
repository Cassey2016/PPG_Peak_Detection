function [output_Vadrevu_1_2019,output_Vadrevu_2_2019] = my_Vadrevu_2019_peakdet(PPG_buffer,fs_PPG)
% =========================================================================
% This is my implementation of this paper:
%
% Vadrevu, Simhadri, and M. Sabarimalai Manikandan. 
% "A robust pulse onset and peak detection method for automated PPG signal 
% analysis system." IEEE Transactions on Instrumentation and Measurement 
% 68.3 (2018): 807-817.
%
% Implemented by Dong Han on 05/03/2020.
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
    debug_flag = false; % decide to plot the paper figure or not.
    %% A. Stationary Wavelet Transform of PPG signal.

    % first, for the input length, you can know the maximum wavelet
    % decomposition level you can get:
    TYPE = '1D'; % extension method.
    MODE = 'zpd'; % zero extension.
    X = PPG_buffer;

    % based on your input signal length, you have to extend your input signal
    % to MATLAB suggested length.
    LEN = 45;%18; % 18 for fs_PPG 50, 45 for fs_PPG 125; for 30 sec input.
    YEXT = wextend(TYPE,MODE,X,LEN); % required by swt.
    sig = YEXT;
    % s = PPG_buffer;
    sLen = length(sig);
    wname = 'bior1.5';
    L = wmaxlev(sLen,wname);

    % [swa,swd] = swt(s,3,'bior1.5'); % the author mentioned wavelet biorthogonal 1.5 (bior1.5)
    [swa,swd] = swt(sig,L,wname); % the author mentioned wavelet biorthogonal 1.5 (bior1.5)
    s1 = swd(3,:) + swd(4,:); 
    s1 = s1(:); % make sure it is column vector.
    s2 = swd(5,:) + swd(6,:) + swd(7,:);
    s2 = s2(:); % make sure it is column vector.

    if debug_flag
        % if you want to debug the result.
        figure;
        t_plot = [1:length(sig)]'./fs_PPG; %
        subplot(5,1,1);
        plot(t_plot,sig);
        xlim([0 t_plot(end)])
        ylabel('Orig');
        title('Fig.3 in TIM 2019 paper');

        subplot(5,1,2)
        plot(t_plot,(swd(1,:) + swd(2,:)))
        xlim([0 t_plot(end)])
        ylabel('s_0');

        subplot(5,1,3);
        plot(t_plot,s1);
        xlim([0 t_plot(end)])
        ylabel('s_1');

        subplot(5,1,4);
        plot(t_plot,s2);
        xlim([0 t_plot(end)])
        ylabel('s_2');

        subplot(5,1,5);
        plot(t_plot,swa(7,:));
        xlim([0 t_plot(end)])
        ylabel('a_7');
    end
    %% B. Multiscale Sum and Products:
    p = s1 .* s2;
    p = p(:);

    if debug_flag
        % if you want to debug the result.
        figure;
        ax(1) = subplot(4,1,1);
        plot(t_plot,sig);
        xlim([0 t_plot(end)])
        ylabel('Orig');
        title('Fig.4 in TIM 2019 paper');

        ax(2) = subplot(4,1,2);
        p1 = swd(1,:) .* swd(2,:) .* swd(3,:) .* swd(4,:) .* swd(5,:) .* swd(6,:) .* swd(7,:); 
        plot(t_plot,p1);
        xlim([0 t_plot(end)])
        ylabel('p_1');

        ax(3) = subplot(4,1,3);
        p1 = swd(3,:) .* swd(4,:) .* swd(5,:) .* swd(6,:) .* swd(7,:); 
        plot(t_plot,p1);
        xlim([0 t_plot(end)])
        ylabel('p_2');

        ax(4) = subplot(4,1,4);
        plot(t_plot,p);
        xlim([0 t_plot(end)])
        ylabel('p');

        linkaxes(ax,'x');
    end
    %% C. Shannon Entropy Envelope Extraction
    eta = 0.01 + std(p);
    p_tilda = abs(p);
    p_tilda(p_tilda < eta) = 0;
    p_tilda = p_tilda(:);

    % normalize p_tilda:
    norm_p_tilda = (p_tilda - min(p_tilda)) ./ (max(p_tilda) - min(p_tilda));
    norm_p_tilda = norm_p_tilda(:);

    se = NaN(size(norm_p_tilda));

    for tttt = 1:size(norm_p_tilda,1)
        if norm_p_tilda(tttt) == 0
            % from MATLAB page: https://www.mathworks.com/help/wavelet/ref/wentropy.html
            %  log(0) = 0
            % 0log(0) = 0.
            se(tttt) = 0;
        else
            se(tttt) = -1 * norm_p_tilda(tttt) .* log(norm_p_tilda(tttt));
        end
    end

    % % method 1: CONV twice:
    filt_Len = floor(0.2 * fs_PPG); % 0.4 is better. 05/04/2020.
    % h = ones(filt_Len,1)/filt_Len;   % A third-order filter has length 4
    % s = conv(se,h,'same'); % return the same size as se
    % s = conv(s,h,'same'); % conv twice

    % method 2: FILTFILT.
    % for 4020, ii = 2, PPG is zero.
    if any(isnan(se))
        % any sample is NaN.
        new_se = se;
        new_se(isnan(new_se)) = [];
        if isempty(new_se)
            % nothing left after removing NaN.

                HR_Vadrevu_1_2019 = 0; % there is no peak location.
                onset_zx = 1;

                HR_Vadrevu_2_2019 = 0; % there is no peak location.
                peak_zx = 1;

                filter_PPG = PPG_buffer;
            output_Vadrevu_1_2019 = struct('filtered_PPG_Vadrevu_2019',filter_PPG,...
                'PPG_peak_loc_Vadrevu_1_2019',onset_zx,...
                'HR_Vadrevu_1_2019',HR_Vadrevu_1_2019);

                output_Vadrevu_2_2019 = struct('filtered_PPG_Vadrevu_2019',filter_PPG,...
                'PPG_peak_loc_Vadrevu_2_2019',peak_zx,...
                'HR_Vadrevu_2_2019',HR_Vadrevu_2_2019);
            return
        else
            % part of data is NaN, maybe I should fill zeros in it?
            keyboard;
        end
    end
    b = ones(filt_Len,1);
    a = -1;
    s = filtfilt(b, a, se); % -> AC component


    %% D. Pulse Peak and Onset Determination.
    % 1. Gaussian derivative kernel:
    sigma_1 = floor(0.05 * fs_PPG); % 0.05 mentioned in the paper.
    M = floor(2 * fs_PPG); % 2 mentioned in the paper.
    g = gausswin(M,sigma_1); % size should be 250 if Fs = 125.
    h_d = diff(g); % g(m+1) - g(m).
    z = conv(s,h_d,'same');

    % % My conv function did not work.
    % temp_z = zeros(size(s,1),1);
    % for nnnn = 1:size(s,1)
    %     for mmmm = 1:size(g,1)-1
    %         if (nnnn-mmmm+1 > 0)
    % %             h_d(mmmm) = g(mmmm+1) - g(mmmm);
    %         temp_z(nnnn) = temp_z(nnnn) + s(mmmm) * h_d(nnnn-mmmm+1);
    %         end
    %     end
    % end

    DownZCi = @(v) find(v(1:end-1) >= 0 & v(2:end) < 0); % Returns Down Zero-Crossing Indices. https://www.mathworks.com/matlabcentral/answers/267222-easy-way-of-finding-zero-crossing-of-a-function
    zx = DownZCi(z); % negative zero crossing point.

    % peak correction algorithm for onset:
    search_intv = floor(0.1 * fs_PPG / 2); % w/2
    onset_zx = NaN(size(zx));
    for zz = 1:size(zx,1)
        temp_zx = zx(zz);
        if temp_zx - search_intv > 0 % not exceed signal limit.
            if temp_zx + search_intv <= size(sig,1)
                temp_PPG = sig(temp_zx - search_intv : temp_zx + search_intv);
                [~,I] = min(temp_PPG);
                if isempty(I) ~= 1
                    adj_loc = temp_zx - search_intv + I - 1;
                else
                    % no local minimum.
                    adj_loc = temp_zx; 
                end
                onset_zx(zz) = adj_loc;
            else
                % right interval exceed signal length.
                onset_zx(zz) = zx(zz);
            end
        else
            % left interval exceed index 1.
            onset_zx(zz) = zx(zz);
        end
    end

    % find peak:
    peak_zx = NaN(size(onset_zx,1)-1,1); % one sample smaller.
    for zz = 2:size(onset_zx,1)
        temp_onset_1 = onset_zx(zz-1);
        temp_onset_2 = onset_zx(zz);
        temp_PPG = sig(temp_onset_1:temp_onset_2);
        [~,I] = max(temp_PPG);
        if isempty(I) ~= 1
           peak_zx(zz-1) = temp_onset_1 + I - 1; % peak is one sample size smaller.
        else
            peak_zx(zz-1) = onset_zx(zz);
        end
    end

    % prepare to output signal:
    filter_PPG = z(LEN+1:end-LEN);
    remove_left = find(onset_zx < LEN+1);
    if isempty(remove_left) ~= 1
        onset_zx(remove_left) = [];
    end
    remove_right = find(onset_zx > size(z,1) - LEN);
    if isempty(remove_right) ~= 1
        onset_zx(remove_right) = [];
    end
    onset_zx = onset_zx - LEN; % shifted.

    remove_left = find(peak_zx < LEN+1);
    if isempty(remove_left) ~= 1
        peak_zx(remove_left) = [];
    end
    remove_right = find(peak_zx > size(z,1) - LEN);
    if isempty(remove_right) ~= 1
        peak_zx(remove_right) = [];
    end
    peak_zx = peak_zx - LEN;

    if debug_flag
        % if you want to debug the result.
        figure;
        ax(1) = subplot(7,1,1);
        plot(t_plot,sig);
        xlim([0 t_plot(end)])
        ylabel('Orig');
        title('Fig.5 in TIM 2019 paper');

        ax(2) = subplot(7,1,2);
        plot(t_plot,p);
        xlim([0 t_plot(end)])
        ylabel('p');

        ax(3) = subplot(7,1,3);
        plot(t_plot,norm_p_tilda);
        xlim([0 t_plot(end)])
        ylabel('p_th');

        ax(4) = subplot(7,1,4);
        plot(t_plot,se);
        xlim([0 t_plot(end)])
        ylabel('se');

        ax(5) = subplot(7,1,5);
        plot(t_plot,s);
        xlim([0 t_plot(end)])
        ylabel('s');

        ax(6) = subplot(7,1,6);
        plot(t_plot,z);
        hold on;
        plot(t_plot(zx),z(zx),'ro');
        xlim([0 t_plot(end)]);
        ylabel('z');


        ax(7) = subplot(7,1,7);
        plot(t_plot,sig);
        hold on;
        plot(t_plot(onset_zx),sig(onset_zx),'go');
        plot(t_plot(peak_zx),sig(peak_zx),'ro');
        xlim([0 t_plot(end)])
        ylabel('orig with peak');

        linkaxes(ax,'x');
    end

    if isempty(onset_zx)
        HR_Vadrevu_1_2019 = 0; % there is no peak location.
        onset_zx = 1;
    else
        HR_Vadrevu_1_2019 = 60 * fs_PPG ./ diff(onset_zx); % calculate the HR.
    end
    
    if isempty(peak_zx)
        HR_Vadrevu_2_2019 = 0; % there is no peak location.
        peak_zx = 1;
    else
        HR_Vadrevu_2_2019 = 60 * fs_PPG ./ diff(peak_zx); % calculate the HR.
    end
    
    output_Vadrevu_1_2019 = struct('filtered_PPG_Vadrevu_2019',filter_PPG,...
        'PPG_peak_loc_Vadrevu_1_2019',onset_zx,...
        'HR_Vadrevu_1_2019',HR_Vadrevu_1_2019);

        output_Vadrevu_2_2019 = struct('filtered_PPG_Vadrevu_2019',filter_PPG,...
        'PPG_peak_loc_Vadrevu_2_2019',peak_zx,...
        'HR_Vadrevu_2_2019',HR_Vadrevu_2_2019);
    
end
