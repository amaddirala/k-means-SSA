function[y]=diag_avg(Xr,L,N);
Frec=zeros(N,1); % reconstructed time series
% diagonal averaging
K=N-L+1;X=Xr;
L0 = min(L,K);
K0 = max(L,K);

if L<K 
    Y=X;
else
    Y=X';
end

for k=1:L0-1
    d=0;
    for m=1:k
        d = d + Y(m,k-m+1);
    end
    Frec(k) = 1/k * d;
end

for k=L0:K0
    d=0;
    for m=1:L0
        d = d + Y(m,k-m+1);
    end
    Frec(k) = 1/L0 * d;
end

for k=K0+1:N
    d=0;
    for m=k-K0+1:N-K0+1
        d = d + Y(m,k-m+1);
    end
    Frec(k) = 1/(N-k+1) * d;
end
y=Frec;

y=y(:)';