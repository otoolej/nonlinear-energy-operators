%-------------------------------------------------------------------------------
% cal_freqweighted_energy: Use a 'nonlinear energy operator' to calculate an
% instantaneous frequency-weighted energy measure
%
% Syntax: [x_nleo]=cal_freqweighted_energy(x,Fs,wlength_ma,bandpass_filter_params)
%
% Inputs: 
%     x           - input signal
%     Fs          - sampling frequency
%     method      - either 'teager','arg',...
%     wlength_ma  - length (in samples) of output moving-average filter; set to [] to
%                   turn off
%     bandpass_filter_params - band-pass filter values; set to [] to turn off
%
% Outputs: 
%     [x_nleo] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 04-04-2014
%
% last update: Time-stamp: <2014-04-06 02:36:45 (otoolej)>
%-------------------------------------------------------------------------------
function [x_nleo]=cal_freqweighted_energy(x,Fs,method,wlength_ma, ...
                                                 bandpass_filter_params)
if(nargin<2 || isempty(Fs)), Fs=1; end
if(nargin<4 || isempty(wlength_ma)), wlength_ma=[]; end
if(nargin<5 || isempty(bandpass_filter_params)), bandpass_filter_params=[]; end


%---------------------------------------------------------------------
% 1. bandpass filter the signal (use zero-phase filter)
%---------------------------------------------------------------------
if(~isempty(bandpass_filter_params))
    LP_fc=bandpass_filter_params(1);  HP_fc=bandpass_filter_params(2);     
    x=do_bandpass_filtering(x,Fs,HP_fc,LP_fc);    
end


%---------------------------------------------------------------------
% 2. calculate frequency-weighted energy
%---------------------------------------------------------------------
switch method
  case 'teager'
    l=0; p=0; q=1; s=-1;
    x_nleo=general_nleo(x,l,p,q,s);    
    
  case 'agarwal'
    l=1; p=2; q=0; s=3;
    x_nleo=general_nleo(x,l,p,q,s);
    
  case 'envelope_diff'
    [x_nleo]=discrete_Hilbert_diff_operator(x);
    
  case 'palmu'
    l=1; p=2; q=0; s=3;
    x_nleo=general_nleo(x,l,p,q,s);
    x_nleo=abs(x_nleo);
    
  case 'abs_teager'
    l=0; p=0; q=1; s=-1;
    x_nleo=general_nleo(x,l,p,q,s);    
    x_nleo=abs(x_nleo);
    
  case 'env_only'
    x_nleo=abs( hilbert(x) ).^2;
    
    
  otherwise 
    error('which method?');
end


%---------------------------------------------------------------------
% 3. smooth with MA filter
%---------------------------------------------------------------------
if(~isempty(wlength_ma))
    x_filter=ma_filter(x_nleo,wlength_ma);
    x_nleo=circshift(x_filter,-ceil(wlength_ma/2));
end



