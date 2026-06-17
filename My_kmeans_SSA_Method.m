function[acap,fd_val,dec]=My_kmeans_SSA_Method(x,M,clusters,Th,Tssa)
%%%% Input arguments %%%%%%%%
%%%% x----> input single channel EEG signal %%%%%
%%%% M----> window length, you can set 400 to 500ms%%%
%%%% clusters-----> Number of clusters. (2 to 4 clusters good enough)
%%%% Th----> Threshold to select eye-blink associated component (step-5)
%%%% resulted by step 3 and 4 in the paper
%%%% Tssa---> Threshold for grouping step of SSA %%%%%
%%%% output arguments %%%%%%
%%%% acap ----> estimated eye-blink artifact
%%%% Written by Ajay Kumar Maddirala and please contact me at maddirala@knu.ac.kr (or) a.maddirala3@gmail.com%%%%%%%%
%%%% Note that this method can work efficiently, where the eye-blink artifact is strong enough in the contaminated EEG signal %%%%%
%%%% Please cite the paper if you use this code %%%%
%%%% How to cite: Maddirala, A.K., Veluvolu, K.C. Eye-blink artifact removal from single channel EEG with k-means and SSA.
%%%% Sci Rep 11, 11043 (2021). https://doi.org/10.1038/s41598-021-90437-7
x=x(:)';N=length(x);K=N-M+1;
for k=1:K
    X(:,k)=x(1,k:k+M-1);
   features(k,1)=((sum(X(:,k).^2)));%%%%% feature-1 %%%%%%
   [features(k,2)]= HjorthParameters(X(:,k));%%%%% feature-2 %%%%%%
   features(k,3)=kurtosis(X(:,k));%%%%% feature-3 %%%%%%
   max1=max(X(:,k));min1=min(X(:,k));features(k,4)=((max1)-abs(min1)).^2;%%%%% feature-4 %%%%%%
end
%lables = kmeans(features,clusters);
lables = kmeans(features,clusters,'start','plus');
%%%%%% implementation of equation (2) in the paper.
for k=1:clusters
temp=zeros(M,K);
arg=find(lables==k);temp(:,arg)=X(:,arg);
dec(k,:)=diag_avg(temp,M,N);
fd_val(k,1)= fd(dec(k,:));
end
%%%%%% implementation of equation (2) in the paper.

arg1=find(fd_val<=Th);
if length(arg1)==1
a=dec(arg1,:);arg=find((a>0|a<0));dummy=zeros(1,length(x));dummy(1,arg)=1;
temp=x.*dummy;[acap]=My_SSA1(temp,M,Tssa);
elseif length(arg1)>1
a=sum(dec(arg1,:));arg=find((a>0|a<0));
dummy=zeros(1,length(x));dummy(1,arg)=1;
temp=x.*dummy;[acap]=My_SSA1(temp,M,Tssa);
else
acap=zeros(1,length(x));
end
