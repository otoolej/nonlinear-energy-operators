%-------------------------------------------------------------------------------
% do_bandpass_filtering: Bandpass filtering using Butterworth (high-pass) and Epilletic
% IIR zero-phase filters (according to reference [1])
%
% Syntax: x=do_bandpass_filtering(x,Fs,LP_fc,HP_fc)
%
% Inputs: 
%     x  - input signal 
%     Fs - sampling frequency
%     LP_fc - low-pass filter cutoff (ignore if =Fs/2 or [])
%     HP_fc - high-pass filter cutoff (ignore if =0 or [])
%
% Outputs: 
%     x - filtered signal
%
% Example:
%     
% [1] Palmu, K., Stevenson, N., Wikström, S., Hellström-Westas, L., Vanhatalo, S., &
% Palva, J. M. (2010). Optimization of an NLEO-based algorithm for automated detection of
% spontaneous activity transients in early preterm EEG. Physiological measurement, 31(11),
% N85–93.

% John M. O' Toole, University College Cork
% Started: 06-04-2014
%
% last update: Time-stamp: <2014-04-06 01:41:51 (otoolej)>
%-------------------------------------------------------------------------------
function x=do_bandpass_filtering(x,Fs,LP_fc,HP_fc)
HP_order=1;
LP_order=6;

if(HP_fc==0),    HP_fc=[]; end
if(LP_fc==Fs/2), LP_fc=[]; end


if(~isempty(HP_fc) && HP_fc>0)
    [b,a]=butter(HP_order,HP_fc/(Fs/2),'high');    
    x=filtfilt(b,a,x);
end

if(~isempty(LP_fc) && LP_fc>0)
    [b,a]=ellip(LP_order,3,50,LP_fc/(Fs/2));
    x=filtfilt(b,a,x);
end
