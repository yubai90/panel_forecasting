function X_new=genx(X,betatilde)

% construct the new X matrix used for the second step in the MCMC algorithm
[T,k]=size(X);
betat=reshape(betatilde,k,T);
X_new=zeros(T,k);

for i=1:T
    X_new(i,:)=X(i,:).*betat(:,i)';
end

X_new=[X X_new];