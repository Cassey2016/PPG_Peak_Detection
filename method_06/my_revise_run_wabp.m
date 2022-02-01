function [r,ssf,my_avg0,A] = my_revise_run_wabp(abp,fs_abp)
% Below was copied from Erick Andres Perez Alday's Github repository 
% " physionetchallenges / matlab-classifier-2020 ": 
% https://github.com/physionetchallenges/matlab-classifier-2020/blob/master/Tools/PhysioNet-Cardiovascular-Signal-Toolbox-master/Tools/BP_Tools/run_wabp.m
% WABP  ABP waveform onset detector.
%   r = run_wabp(abp) obtains the onset time (in samples) 
%       of each beat in the ABP waveform.
%
%   In:   ABP (125Hz sampled)
%   Out:  Onset sample time
% 
%   Usage:
%   - ABP waveform must have units of mmHg
%
%   Written by James Sun (xinsun@mit.edu) on Nov 19, 2005.  This ABP onset
%   detector is adapted from Dr. Wei Zong's wabp.c.
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

% Dong changed: input should be 250 Hz for filtering.
%% Input checks
    % if nargin ~=1
    %     error('exactly 1 argment needed');
    % end

    if size(abp,2)~=1
        error('Input must be a <nx1> vector');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % scale physiologic ABP
    offset   = 1600;
    scale    = 20;
    Araw     = abp*scale-offset;

    % LPF
    A = filter([1 0 0 0 0 -2 0 0 0 0 1],[1 -2 1],Araw)/24+30;
    A = (A(4:end)+offset)/scale; % Takes care of 4 sample group delay

    % ------- Dong changed this: -------
    A = A ./ std(A); % normalizing data is very important for my peak detection.
    A = A - mean(A);

    % Slope-sum function
    dypos          = diff(A);
    dypos(dypos<0) = 0;
    % ssf            = [0; 0; conv(ones(16,1),dypos)];
    w = 16/125*fs_abp; % 125 Hz to 250 Hz.
    ssf            = [0; 0; conv(ones(w,1),dypos)];

    % Decision rule
    first_8sec = 8*fs_abp;
    % avg0       = sum(ssf(1:1000))/1000;   % average of 1st 8 seconds (1000 samples) of SSF
    avg0       = sum(ssf(1:first_8sec))/first_8sec;
    Threshold0 = 3*avg0;                  % initial decision threshold

    % ignoring "learning period" for now
    lockout    = 0;    % lockout >0 means we are in refractory
    timer      = 0;
    % z          = zeros(100000,1);
    z          = zeros(fs_abp*800,1);
    counter    = 0;

    % Dong: copied from wabp.c, 02/27/2020. % Dong change here. 02/27/2020.
    TmDEF = 0.25; %5;% Dong change here. 02/27/2020.
    max_min_thres = 0.1; %10;% Dong change here. 02/27/2020.
    my_avg0 = zeros(size(abp));% Dong change here. 02/27/2020.
    step_adjust_thres = 0.025; % it was 0.1 % Dong change here. 02/27/2020.
    % for t = 50:length(ssf)-17
    for t = round(0.4*fs_abp):length(ssf)-w-1
        lockout = lockout - 1;
        timer   = timer   + 1;      % Timer used for counting time after previous ABP pulse

        if (lockout<1) & (ssf(t)>avg0+TmDEF) %(ssf(t)>avg0+5)  % Not in refractory and SSF has exceeded threshold here % Dong change here. 02/27/2020.
            timer = 0;
            maxSSF = max(ssf(t:t+w));  % Find local max of SSF
            minSSF = min(ssf(t-w:t));  % Find local min of SSF
            if maxSSF > (minSSF + max_min_thres) %(minSSF + 10)% Dong change here. 02/27/2020.
                onset = 0.01*maxSSF ;  % Onset is at the time in which local SSF just exceeds 0.01*maxSSF

                tt       = t-w:t;
                dssf     = ssf(tt) - ssf(tt-1);
                BeatTime = find(dssf<onset,1,'last')+t-w-1;
                counter  = counter+1;

                if isempty(BeatTime)
                    counter = counter-1;
                else
                    z(counter) = BeatTime;
                end
                Threshold0 = Threshold0 + step_adjust_thres*(maxSSF - Threshold0);  % adjust threshold
                avg0 = Threshold0 / 3;        % adjust avg

                lockout = round(32/125*fs_abp);   % lock so prevent sensing right after detection (refractory period)
            end
        end

        if timer > round(312/125*fs_abp)  % Lower threshold if no pulse detection for a while
            Threshold0 = Threshold0 - 0.1; %Threshold0 - 1; % Dong change here. 02/27/2020.
            avg0       = Threshold0/3;
        end
        my_avg0(t,1) = avg0+TmDEF; % % Dong change here. 02/27/2020.
    end
    r = z(find(z))-2;
end