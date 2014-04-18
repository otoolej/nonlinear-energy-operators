%-------------------------------------------------------------------------------
% plot_eeg_examples: Plot some examples of the EEG
%
% Syntax: []=plot_eeg_examples(PRINT_PLOTS)
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
% Started: 04-04-2014
%
% last update: Time-stamp: <2014-04-18 11:09:17 (otoolej)>
%-------------------------------------------------------------------------------
function []=plot_eeg_examples(PRINT_PLOTS)
if(nargin<1 || isempty(PRINT_PLOTS)), PRINT_PLOTS=0; end


nleo_parameters;


if(~exist(DATA_DIR,'dir'))
    error('EEG data is needed for this function.');
end


b=load([D_NLEO_SIM 'eeg_segments_normterms.mat']);
x_normterms=b.eeg_data{1}; x_mask=b.anno{1}; Fs=b.Fs{1};
x_normterms=do_bandpass_filtering(x_normterms,Fs,LOW_PASS_DATA,[]);

b=load([D_NLEO_SIM 'eeg_segments_preterms.mat']);
x_preterms=b.eeg_data{7}; x_mask_pre=b.anno{7}; Fs=b.Fs{1};
x_preterms=do_bandpass_filtering(x_preterms,Fs,LOW_PASS_DATA,[]);


figure(1); clf; hold all;
nn=(0:30*Fs);
xx=x_normterms(nn+1);
t=nn./Fs;    
hp(1)=plot(t,xx);
xm=x_mask(nn+1);
xm(xm==0)=NaN;
hm=plot(t,xm-110,'linewidth',7,'color',[1 1 1].*0.3);

plot_bias=350;
xx=x_preterms(nn+1);
hp(2)=plot(t,xx-plot_bias,'b');
xm=x_mask_pre(nn+1);
xm(xm==0)=NaN;
hm=plot(t,xm+min(xx)-20-plot_bias,'linewidth',7,'color',[1 1 1].*0.3);


set(gca,'xtick',[0:10:50]);
yl=ylim;
set(gca,'ylim',[yl(1) 70]);
set(gca,'ytick',[-(plot_bias+100) -plot_bias -(plot_bias-100) -70 0 70]);
set(gca,'yticklabel',{'-100','0','100','-70', '0', '70'});


xlabel('time (seconds)');
ylabel('voltage (\muV)');
set_gca_fonts('Arial',12,gca);

set(hp,'linewidth',0.6);
if(PRINT_PLOTS)
    %    set(gca,'Position',[0 0 1 1]);
% $$$     set(gcf, 'Units','centimeters', 'Position',[0 0 10 6]);
% $$$     set(gcf, 'PaperPositionMode','auto');
    
    
    print2eps([PIC_DIR 'eeg_preterms_normterms_examples.eps']); %,'-depsc2');
end

