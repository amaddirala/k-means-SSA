function[a]=My_SSA1(x,M,TSSA)
x=x(:)';N=length(x);
K=N-M+1;
for k=1:K
    X(:,k)=x(1,k:k+M-1);
end
C=X*X';
[U,D]=svd(C);en=diag(D)/sum(diag(D));
arg=find(en>=TSSA);
a=diag_avg(U(:,arg)*U(:,arg)'*X,M,N);