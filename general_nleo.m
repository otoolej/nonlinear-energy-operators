%-------------------------------------------------------------------------------
% general_nleo: 'General' NLEO expression: Ψ(n)=x(n-l)x(n-p)-x(n-q)x(n-s)
%                for l+p=q+s  (and [l,p]≠[q,s], otherwise Ψ(n)=0)
%
% Syntax: x_nleo=general_nleo(x,l,p,q,s)
%
% Inputs: 
%     x       - input signal
%     l,p,q,s - parameters for the operator; integer values and must satisfy: 
%                 l+p=q+s and [l,p]≠[q,s]
%
% Outputs: 
%     x_nleo - output nonlinear energy operator
%
% Example:
%     % generate test signals:
%     N=256; n=0:N-1;
%     w1=pi/(N/32); ph1=-pi+2*pi*rand(1,1);  a1=1.3;
%     x1=a1.*cos(w1.*n + ph1);  
%
%     % compute instantaneous energy:
%     x_nleo=general_nleo(x1,1,2,0,3);
%
%     % plot:
%     figure(1); clf; 
%     subplot(211); plot(x1); ylabel('amplitude');
%     subplot(212); plot(x_nleo,'-'); ylim([0 0.5])
%     ylabel('energy');


% John M. O' Toole, University College Cork
% Started: 28-01-2014
%-------------------------------------------------------------------------------
function x_nleo=general_nleo(x,l,p,q,s)
if(nargin<2 || isempty(l)), l=1; end
if(nargin<3 || isempty(p)), p=2; end
if(nargin<4 || isempty(q)), q=0; end
if(nargin<5 || isempty(s)), s=3; end


if( (l+p)~=(q+s) && any(sort([l p])~=sort([q s])) )
    warning('Incorrect parameters for NLEO.');
end


N=length(x);
x_nleo=zeros(1,N);

iedges=abs(l)+abs(p)+abs(q)+abs(s)+1;
n=(1+iedges):(N-iedges);


x_nleo(n)=x(n-l).*x(n-p) - x(n-q).*x(n-s); 
