function [output_Elgendi_1_2013] = my_Elgendi_2013_method_I_peakdet(raw_PPG, delta, fs_PPG)
% -------------------------------------------------------------------------
% Dong add this on 02/25/2020, based on this paper:
% Elgendi, Mohamed, et al. 
% "Systolic peak detection in acceleration photoplethysmograms measured from 
% emergency responders in tropical conditions." PLoS One 8.10 (2013): e76585.
% 
% (1): bandpass filter (0.5-8Hz)
    [b, a] = butter(6,[0.5 8]/(fs_PPG/2)); % bandpass filter 0.5-10Hz, changed from 0.5-20 to 0.5-9 Hz at 11/21/2018
    raw_PPG = filtfilt(b, a, raw_PPG); % -> AC component
    raw_PPG = raw_PPG ./ std(raw_PPG); % normalizing data is very important for my peak detection.
    raw_PPG = raw_PPG - mean(raw_PPG);
    
    debugging_plot_flag = false; % only for plotting debugging figures. 
% -------------------------------------------------------------------------
% Below code is copied from: http://billauer.co.il/peakdet.html
% PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%      
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.
    maxtab = [];
    mintab = [];

    raw_PPG = raw_PPG(:); % Just in case this wasn't a proper vector

    % if nargin < 3
      x = (1:length(raw_PPG))';
    % else 
    %   x = x(:);
    %   if length(raw_PPG)~= length(x)
    %     error('Input vectors v and x must have same length');
    %   end
    % end

    if (length(delta(:)))>1
      error('Input argument DELTA must be a scalar');
    end

    if delta <= 0
      error('Input argument DELTA must be positive');
    end

    mn = Inf; mx = -Inf;
    mnpos = NaN; mxpos = NaN;

    lookformax = 1;

    for i=1:length(raw_PPG)
      this = raw_PPG(i);
      if this > mx, mx = this; mxpos = x(i); end
      if this < mn, mn = this; mnpos = x(i); end

      if lookformax
        if this < mx-delta
          maxtab = [maxtab ; mxpos mx];
          mn = this; mnpos = x(i);
          lookformax = 0;
        end  
      else
        if this > mn+delta
          mintab = [mintab ; mnpos mn];
          mx = this; mxpos = x(i);
          lookformax = 1;
        end
      end
    end



    if isempty(maxtab)
        HR_Elgendi_1_max_2009 = 0; % there is no peak location.
        peak_loc_max = 1;
    else
        peak_loc_max = maxtab(:,1);
        HR_Elgendi_1_max_2009 = 60 * fs_PPG ./ diff(peak_loc_max); % calculate the HR.
    end

    if isempty(mintab)
        HR_Elgendi_1_min_2009 = 0; % there is no peak location.
        peak_loc_min = 1;
    else
        peak_loc_min = mintab(:,1);
        HR_Elgendi_1_min_2009 = 60 * fs_PPG ./ diff(peak_loc_min); % calculate the HR.
    end

    output_Elgendi_1_2013 = struct('filtered_PPG_Elgendi_1_2013',raw_PPG,...
        'PPG_peak_loc_Elgendi_1_max_2013',peak_loc_max,...
        'PPG_peak_loc_Elgendi_1_min_2013',peak_loc_min,...
        'HR_Elgendi_1_max_2013',HR_Elgendi_1_max_2009,...
        'HR_Elgendi_1_min_2013',HR_Elgendi_1_min_2009);


    if debugging_plot_flag % debugging plot
        figure;
        plot(x,raw_PPG);hold on;
        plot(peak_loc_max,raw_PPG(peak_loc_max),'ro');
        plot(peak_loc_min,raw_PPG(peak_loc_min),'go');
    end
end