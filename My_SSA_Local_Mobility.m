function[a]=My_SSA_Local_Mobility(x,M,fs,f)
x=x(:)';N=length(x);
[th]=My_signal_complexity(sin(2*pi*[0:M-1]*f/fs));
K=N-M+1;
for k=1:K
    X(:,k)=x(1,k:k+M-1);
end
C=X*X';
[U,D]=svd(C);
for k=1:M
[sm(k)]=My_signal_complexity(U(:,k));
end
arg=find(sm<=th);
a=diag_avg(U(:,1:2)*U(:,1:2)'*X,M,N);
