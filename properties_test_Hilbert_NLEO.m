%-------------------------------------------------------------------------------
% properties_test_Hilbert_NLEO: Verify properties for 'envelope--derivative' 
% operator
%
% Syntax: properties_test_Hilbert_NLEO(prop_numb)
%
% Inputs: 
%     prop_numb - either: 0,1,2,4
%
% Example:
%     >> properties_test_Hilbert_NLEO(1);
%

% John M. O' Toole, University College Cork
% Started: 25-03-2014
%
% last update: Time-stamp: <2014-06-09 09:52:24 (otoolej)>
%-------------------------------------------------------------------------------
function []=properties_test_Hilbert_NLEO(prop_numb,PRINT_PLOTS)
if(nargin<1 || isempty(prop_numb)), prop_numb=4; end
if(nargin<2 || isempty(PRINT_PLOTS)), PRINT_PLOTS=0; end


% for plotting:
FANCY_PLOT=1;
FONT_NAME='Arial';
FONT_SIZE=12;

nleo_parameters;


if(PRINT_PLOTS), FANCY_PLOT=1; end
method_teager_str=sprintf('%s%s%s','Teager','-','Kaiser');
method_hilbert_str=sprintf('%s%s%s','envelope','-','derivative');
leg_pos=[];

Fs=1;

%---------------------------------------------------------------------
% 0. define some simple test signals of the form:  a₁cos(ω₁t + φ₁)
%---------------------------------------------------------------------
N=256;
n=0:N-1;
w1=pi/(N/32); w2=pi/(N/8);
ph1=-pi+2*pi*rand(1,1); 
ph2=-pi+2*pi*rand(1,1);
a1=1.3; a2=3.1;
x1=a1.*cos(w1.*n + ph1); 
x2=a2.*cos(w2.*n + ph2);


x=[]; extact=[];
switch prop_numb
    
  case 0
    %---------------------------------------------------------------------
    % concatentated (in time) signals: a₁cos(ω₁t + φ₁) PLUS a₂cos(ω₂t + φ₂) 
    %---------------------------------------------------------------------
    Nh=ceil(N/2);
    a1=1; a2=1;
    w1=pi/(N/32); w2=pi/(N/16);
    
    x{1}=[a2.*cos(w2.*(Nh:N)) a1.*cos(w1.*(1:Nh))];
    extact{1}=[ones(1,Nh).*(a2^2)*sin(w2)^2 ones(1,Nh).*(a1^2)*sin(w1)^2];
    extact_tk{1}=[ones(1,Nh).*(a2^2)*sin(w2)^2 ones(1,Nh).*(a1^2)*sin(w1)^2];    
    leg_pos=[0.57 0.38];        
    
    % for plotting:    
    ysig_ticks=[-1 0 1];     ysig_limits=[-1.4 1.4];
    ynleo_ticks=[0 0.1 0.2]; ynleo_limits=[0 0.27];

    
  case 1
    %---------------------------------------------------------------------
    % stationary signal: a₁cos(ω₁t + φ₁)
    %---------------------------------------------------------------------
    x{1}=x2./max(x2);
    extact{1}=ones(1,N).*(a1^2)*(w2^2);    
    extact_tk{1}=ones(1,N).*(a1^2)*sin(w2)^2;        
    
    % for plotting:    
    ysig_ticks=[-1 0 1];     ysig_limits=[-1.4 1.4];
    ynleo_ticks=[0 0.1 0.2]; ynleo_limits=[0 0.27];
    leg_pos=[0.55 0.3];    
    
    
  case 2
    %---------------------------------------------------------------------
    % signal with amplitude modulation: e^{rt} a₁cos(ω₁t + φ₁)
    %---------------------------------------------------------------------
    r1=0.005;
    x{1}=exp(-r1.*n).*x1;
    extact{1}=(a1^2).*exp(-2*r1.*n).*(sin(w1)^2+r1^2);
    extact_tk{1}=(a1^2).*exp(-2*r1.*n).*sin(w1)^2;    
    
    % for plotting:    
    ysig_ticks=[-1 0 1];     ysig_limits=[-1.4 1.4];
    ynleo_ticks=[0 0.1 0.2]; ynleo_limits=[0 0.27];
    leg_pos=[0.55 0.33];    
    
  case 3
    %---------------------------------------------------------------------
    % signal with frequency modulation: a₁cos(φ(t))
    %---------------------------------------------------------------------

    if_law=0.1+0.3*sin(n.*pi/N);
    ph=cumsum(if_law);
    x{1}=a1.*cos( ph + ph1 );
    extact{1}=(a1^2).*sin(if_law).^2;
    
    extact_tk{1}=extact{1};
    
    % for plotting:
    ysig_ticks=[-1 0 1];     ysig_limits=[-1.4 1.4];
    ynleo_ticks=[0 0.1 0.2]; ynleo_limits=[0 0.27];
    leg_pos=[0.35 0.2];

    
  case 4
    %---------------------------------------------------------------------
    % sum of two signals: a₁cos(ω₁t + φ₁) + a₂cos(ω₂t + φ₂) 
    %---------------------------------------------------------------------
    x{1}=x1+x2;
    
    extact{1}=(a1^2)*(w1^2) + (a2^2)*(w2^2) + ...
              2*a1*a2*w1*w2.*cos( n.*(w1-w2)+ph1-ph2 );
    
    extact_tk{1}=(a1^2)*(w1^2) + (a2^2)*(w2^2) + ...
        a1*a2*(1-cos(w1+w2)).*cos( n.*(w1-w2)+ph1-ph2 );

    % for plotting:    
    ysig_ticks=[-3 0 3]; ysig_limits=[-4.5 4.5];
    ynleo_ticks=[0 1 2]; ynleo_limits=[-0.3 1.2];
    leg_pos=[0.57 0.38];    
    
  otherwise
    error('input argument either: 0,1,2,3, or 4');
end


if(FANCY_PLOT)
    subplot=@(m,n,p) subtightplot(m,n,p,[0.035,0.2],[0.15,0.1],[0.1,0.1]);
end


%---------------------------------------------------------------------
% 1. generate operators and plot
%---------------------------------------------------------------------
for n=1:length(x)
      
    % generate the NLEO functions:
    xn_hilbert{n}=cal_freqweighted_energy(x{n},Fs,'envelope_diff');
    xn_teager{n}=cal_freqweighted_energy(x{n},Fs,'teager');

    % plot:
    figure(n); clf; hold all;
    hax(1)=subplot(2,1,1); hold all; hp(1)=plot(x{n},'color',lcolor{1}); 
    npart=4:N-4;
    hax(2)=subplot(2,1,2); hold all;  hp(2)=plot(npart,xn_hilbert{n}(npart), ...
                                                 '--','color',lcolor{3});  
    hp(3)=plot(npart,xn_teager{n}(npart),'-','color',lcolor{2});  
    
    h=legend(method_hilbert_str,method_teager_str);
    if(~isempty(leg_pos))
        cur_pos=get(h,'position');
        set(h,'position',[leg_pos(1) leg_pos(2) cur_pos(3) cur_pos(4)]);
% $$$     legend(h,'boxoff');
    end


    if(FANCY_PLOT)
        set(hax(1),'ytick',ysig_ticks); 
        set(hax(1),'xticklabel',{[]});         
        set(hax(1),'xlim',[0 N]);                        
        set(hax(1),'ylim',ysig_limits); 
        set(hp,'linewidth',1.4);
        ylabel(hax(1),'amplitude');

        set(hax(2),'ytick',ynleo_ticks);         
        set(hax(2),'ylim',ynleo_limits);        
        set(hax(2),'xlim',[0 N]);
        xlabel(hax(2),'time (samples)'); 
        ylabel(hax(2),'energy');        
        set_gca_fonts(FONT_NAME,FONT_SIZE,hax(1));
        set_gca_fonts(FONT_NAME,FONT_SIZE,hax(2));               
        
        
        if(PRINT_PLOTS)
            print2eps([PIC_DIR 'testing_property_no' num2str(prop_numb) '.eps']);
        end
    end
end

