function[xf]=My_FBSE_EWT(x,fs,L,P)
Fs=fs;x=x(:);
f=x;N=length(f);params.SamplingRate = Fs; 
params.globtrend = 'none';
params.degree=8; % degree for the polynomial interpolation
params.reg = 'gaussian';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;
params.detect = 'scalespace';
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans
params.N = 10; % maximum number of bands
params.completion = 0; % choose if you want to force to have params.N modes
params.log=1;
Bound=1;   % Display the detected boundaries on the spectrum
subresf=1;
boundaries=[0.16 0.315 0.510 1.2 3.0]; %physionet
ff=fft(f);
mfb=EWT_Meyer_FilterBank(boundaries,length(ff));
ewt=cell(length(mfb),1);
for k=1:length(mfb)
    ewt{k}=real(ifft(conj(mfb{k}).*ff));
    ewt{k}=ewt{k}(1:end);
end
Bound=1;
xxx=(linspace(0,1,round(length(mfb{1,1}))))*Fs;
% if Bound==1 %Show the boundaries on the spectrum
%     div=1;
%     if (strcmp(params.detect,'adaptive')||strcmp(params.detect,'adaptivereg'))
%         Show_EWT_Boundaries(abs(fft(f)),boundaries,div,params.SamplingRate,InitBounds);
%     else
%         Show_EWT_Boundaries(abs(fft(f)),boundaries,div,params.SamplingRate);
%     end
% end
xa=ewt{1,1};xd1=ewt{2,1};
xd2=ewt{3,1};xd3=ewt{4,1};
xd4=ewt{5,1};xd5=ewt{6,1};
% lambda = 600;       % lambda : TV regularization parameter

deg = 1;            % deg : degree of polynomial
NB=(N-L)/(L-P)+1 ;      % This is the number of blocks - it should be an
                     % integer, otherwise the data will be truncated
mu0 = 50;
mu = 0.5;
Nit = 300;
p = 0.7;
E = 1e-8;
lambda = 900;
[s1, s2, cost] = lopatv_Lp(xa, L, P, deg, lambda, Nit, mu0, mu, p, E);
kk=xa-s1-s2;
xfn=kk+xd1+xd2+xd3+xd4+0*xd5;
xf=(x-xfn)';
