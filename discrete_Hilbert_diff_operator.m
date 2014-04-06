%-------------------------------------------------------------------------------
% discrete_Hilbert_diff_operator: Discrete-implementation of the 'frequency-weighted
% envelope' operator.  Use central-finite method for differentiation and discrete Hilbert
% transform (of length-N) for Hilbert transform
%
% Γ[x(n)] = y(n)² + H[y(n)]², where y(n) is the derivative of x(n) and H[] is the
% Hilber transform
%
% Syntax: [x_nleo]=discrete_Hilbert_diff_operator(x)
%
% Inputs: 
%     x - input signal
%
% Outputs: 
%     [x_nleo] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 03-04-2014
%
% last update: Time-stamp: <2014-04-04 16:35:12 (otoolej)>
%-------------------------------------------------------------------------------
function [x_nleo]=discrete_Hilbert_diff_operator(x)


%---------------------------------------------------------------------
% 1. insert zeros at start and end:
%---------------------------------------------------------------------
Nstart=length(x);
if(rem(length(x),2)~=0), x=[x 0]; end
N=length(x);


nl=2:N-1;
xx=zeros(1,N);

%---------------------------------------------------------------------
% 2. calculate the Hilber transform
%---------------------------------------------------------------------
h=discrete_hilbert(x);


%---------------------------------------------------------------------
% 3. implement with the central finite difference equation
%---------------------------------------------------------------------
xx(nl)=(x(nl+1).^2 + x(nl-1).^2 + h(nl+1).^2 + h(nl-1).^2)./4 - ...
       (x(nl+1).*x(nl-1) + h(nl+1).*h(nl-1))./2;

x_nleo=xx(3:end-2);
x_nleo=[0 0 xx(3:end-2) 0 0];


x_nleo=x_nleo(1:Nstart);


%---------------------------------------------------------------------
% DEBUG: testing and compare 
%---------------------------------------------------------------------
% $$$ if(DBplot)
% $$$     xtest=zeros(1,N); xd=zeros(1,N);
% $$$     xd(nl)=(x(nl+1)-x(nl-1))./2;
% $$$     xtest=xd.^2 + discrete_hilbert(xd).^2;
% $$$ 
% $$$     xtest=xtest(3:end-2);
% $$$     
% $$$     figure(11); clf; hold all;
% $$$     plot(xx); plot(xtest); plot(xx-xtest);
% $$$     dispEE(xx,xtest);
% $$$     dispSize(xx,xtest);
% $$$ end



function x_hilb=discrete_hilbert(x)
%---------------------------------------------------------------------
% Discrete Hilbert transform
%---------------------------------------------------------------------
N=length(x); Nh=ceil(N/2);
k=0:N-1;

H=-j.*sign(Nh-k).*sign(k);

%H=fft( ifft(H).*padWin(shiftWin(hamming(3)),N).' );

x_hilb=ifft( fft(x).*H );


