function[dec]=my_ssa_dec(x,M,fs,f,order)
N=length(x);
K=N-M+1;
for k=1:K
    X(:,k)=x(1,k:k+M-1);
end
C=X*X';
[U,D]=svd(C);
for k=1:M
    [p1,f1]=pburg(U(:,k),order,2^10,fs);
    [max1,arg1]=max(p1);temp(k)=f1(arg1);
end
arg=find(temp>f);
dec(1,:)=diag_avg(U(:,arg)*U(:,arg)'*X,M,N);%%% high frequency component
dec(2,:)=x-dec(1,:);%%%low frequency component
