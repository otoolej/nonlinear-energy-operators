%-------------------------------------------------------------------------------
% bias_of_estimators: Estimate bias of operators by approximating expectation 
% operator. Using 10,000 iterations of white Gaussian noise.
%
% Syntax: []=bias_of_estimators(PRINT_PLOTS)
%
% Inputs: 
%      PRINT_PLOTS - print or not, either 0 or 1; 
%
% Outputs: 
%     [] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 28-03-2014
%
% last update: Time-stamp: <2014-06-06 18:20:27 (otoolej)>
%-------------------------------------------------------------------------------
function []=bias_of_estimators(PRINT_PLOTS)
if(nargin<1 || isempty(PRINT_PLOTS)), PRINT_PLOTS=0; end


%---------------------------------------------------------------------
% 0. set parameters
%---------------------------------------------------------------------
nleo_parameters;
FONT_NAME='Arial';
FONT_SIZE=12;

% number of iterations (L) and signal length (N) of white Gaussian noise:
N=64; L=10000;
xx=randn(L,N);

% set to 1 to take absolute value of the operator:
ABS_NLEO=0;
% set to 1 to up-sampling signal before analysis:
RESAMPLE_X=1;

Nup=N*4; Nup2=ceil(N*6.33);
nplain_teag=zeros(L,Nup);
nenv_diff=zeros(L,Nup);
nplain_got=zeros(L,Nup2);


%---------------------------------------------------------------------
% 1. iterate analysis and estimate Expectation operator by taking the 
%    mean value
%---------------------------------------------------------------------
for n=1:L
    
    % up-sample by a factor of 4 for the Teager--Kaiser and envelope-derivate 
    % operators:
    if(RESAMPLE_X)
        xxf=resample(xx(n,:),ceil(4*N),N);
    else
        xxf=xx(n,:);
    end
    xxf=xxf./std(xxf);

    % compute the Teager--Kaiser and envelope-derivate operators:
    nenv_diff(n,:)=cal_freqweighted_energy(xxf,1,'envelope_diff');
    nplain_teag(n,:)=cal_freqweighted_energy(xxf,1,'teager');

    
    % up-sample by a factor of 6.33 for the Agarwal--Gotman operator:
    if(RESAMPLE_X)
        xxf=resample(xx(n,:),ceil(6.33*N),N);
        xxf=xxf./std(xxf);        
    end

    % compute the Agarwal--Gotman operator:    
    nplain_got(n,:)=cal_freqweighted_energy(xxf,1,'agarwal');
end
if(ABS_NLEO)
    nplain_teag=abs(nplain_teag);
    nplain_got=abs(nplain_got);        
end


%---------------------------------------------------------------------
% 2. PLOT the expectation values (i.e. mean-values)
%---------------------------------------------------------------------
Nteag=size(nplain_teag,2);
npart=4:Nteag-4;
Ngot=size(nplain_got,2);
npart_got=8:Ngot-8;

figure(1); clf; hold all;
hp(1)=plot( npart_got, mean(nplain_got(:,npart_got)),'-.');
hp(3)=plot( npart, mean(nplain_teag(:,npart)));
hp(2)=plot( npart, mean(nenv_diff(:,npart)),'--');


line(xlim,[1 1],'linestyle','-','linewidth',3.5,'color',[1 1 1].*0.5);


xlabel('time (samples)'); ylabel('energy');
h=legend('Agarwal-Gotman','envelope-derivative','Teager-Kaiser','\sigma^2');
legend(h,'location','best');
legend('boxoff');

if(RESAMPLE_X)
    set(gca,'ylim',[0.24 0.37]); 
    set(gca,'ytick',[0.25 0.3 0.35]);    
    set(gca,'ylim',[0 1.2]); 
    set(gca,'ytick',[0:0.2:1]);    
    
    
    set(gca,'xlim',[0 400]);    
    set(gca,'xtick',[0:100:400]);
    set(hp,'linewidth',1.4);
    set(hp(1),'color',lcolor{1}); set(hp(2),'color',lcolor{2});
    set(hp(3),'color',lcolor{3});    
    dd=get(findobj(gcf,'tag','legend'),'position');
    set(findobj(gcf,'tag','legend'),'position',[0.6, 0.45, dd(3:4)])
    
    set_gca_fonts(FONT_NAME,FONT_SIZE,gca);
    
end


if(PRINT_PLOTS)
    print2eps([PIC_DIR 'bias_resampled.eps']);
end


% $$$ return;
% $$$ figure(2); clf; hold all;
% $$$ subplot(221); plot( nplain_teag.' );
% $$$ subplot(222); plot( nplain_got.' );
% $$$ subplot(223); plot( nenv_diff.' );
% $$$ subplot(224); plot( xxf );
% $$$ 
% $$$ figure(3); clf; hold all;
% $$$ pfft(xxf,1);

    
