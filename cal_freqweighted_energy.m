%-------------------------------------------------------------------------------
% cal_freqweighted_energy: Use a 'nonlinear energy operator' to calculate an
% instantaneous frequency-weighted energy measure
%
% Syntax: [x_nleo]=cal_freqweighted_energy(x,Fs,wlength_ma,bandpass_filter_params)
%
% Inputs: 
%     x           - input signal
%     Fs          - sampling frequency
%     method      - either 'teager','agarwal', 'envelope_diff', 'palmu', 'abs_teager',
%                   or 'env_only'
%     wlength_ma  - length (in samples) of output moving-average filter; set to [] to
%                   turn off
%     bandpass_filter_params - band-pass filter values; set to [] to turn off
%
% Outputs: 
%     [x_nleo] - output operator
%
% Example:
%
%     % generate two sinusoidal signals:
%     N=256; n=0:N-1;
%     w1=pi/(N/32); ph1=-pi+2*pi*rand(1,1);  a1=1.3;
%     w2=pi/(N/8); ph2=-pi+2*pi*rand(1,1);  a2=3.1;
%     x1=a1.*cos(w1.*n + ph1);  x2=a2.*cos(w2.*n + ph2);
%     x=x1+x2;
%
%     % compute instantaneous energy:
%     x_env_diff=cal_freqweighted_energy(x,1,'envelope_diff');
%     x_teager  =cal_freqweighted_energy(x,1,'teager');
%
%     % plot:
%     figure(1); clf; 
%     subplot(211); plot(x); ylabel('amplitude');
%     subplot(212); hold all; plot(x_env_diff,'-'); plot(x_teager,'--');
%     ylabel('energy');
%     legend('envelope-derivative','Teager-Kaiser');
%


% John M. O' Toole, University College Cork
% Started: 04-04-2014
%
% last update: Time-stamp: <2014-06-05 14:50:34 (otoolej)>
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
    x_nleo=abs( (x) ).^2;
    
    
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



