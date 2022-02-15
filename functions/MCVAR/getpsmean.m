function [mu,Sig] = getpsmean(y,X,ht,invA)
% Compute some quantities used for drawing posterior mean for VAR
% coefficients

[T,NG]=size(y);
mu=0;
Sig=0;
for t=1:T
    Sigmatemp=invA*sparse(1:NG,1:NG,ht(t,:))*invA';
    Sigmainvt=Sigmatemp\speye(NG);
    mutemp=X(t,:)'*y(t,:)*Sigmainvt;
    mu=mu+mutemp;
    Sigtemp=kron(Sigmainvt,X(t,:)'*X(t,:));
    Sig=Sig+Sigtemp;
end
mu=mu(:);
