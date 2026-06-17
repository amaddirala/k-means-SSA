function[a]=My_SSA_ICA(x,M,levels,fs,order)
temp=x;f=fs/4;
for k=1:levels
[X]=my_ssa_dec(temp,M,fs,f,order);dec(k,:)=X(1,:);
temp=X(2,:);f=f/2;
end
dec(k+1,:)=temp;%%%% decomposed components from signals x; they are input to 
[ics,A,W]=fastica(dec,'approach','defl','g','tanh');
[r,c]=size(ics);
for k=1:r
S(k,:)=sum(A(:,k)*ics(k,:));
end
[a,s,in]=My_Plots_function(x,dec);
