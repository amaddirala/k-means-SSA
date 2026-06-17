function[a]=My_EEMD_ICA(x,nstd,ne)
x=x(:)';N=length(x);
[modes]=eemd(x,nstd,ne);modes=modes';modes=modes(2:end,:);
[ics,A,W]=fastica(modes,'approach','symm','g','tanh');
[r,c]=size(ics);figure
for k=1:r
    s(k,:)=sum(A(:,k)*ics(k,:));fd_val(k,1)= fd(s(k,:));
    subplot(ceil(r/2),2,k);plot(x);hold on;plot(s(k,:));title([num2str(k)],'FontWeight','normal')
end
set(gcf, 'Position',[200 100 1500 800]);
in=input('Enter argument of artifact component:');
if length(in)==1
    a=s(in,:);
else
    a=sum(s(in,:));
end

% 
