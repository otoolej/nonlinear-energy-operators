%-------------------------------------------------------------------------------
% general_nleo: 'General' NLEO expression: Ψ(n)=x(n-l)x(n-p)-x(n-q)x(n-s)
%                for l+p=q+s  (and [l,p]≠[q,s], otherwise Ψ(n)=0)
%
% Syntax: x_nleo=general_nleo(x,l,p,q,s)
%
% Inputs: 
%     x,l,p,q,s - 
%
% Outputs: 
%     x_nleo - 
%
% Example:
%     
%

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
