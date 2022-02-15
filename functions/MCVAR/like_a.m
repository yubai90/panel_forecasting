function loglike = like_a(a,b,c,lambda,beta0,hyp1,hyp2)
% compute the log-likelihood function of the posterior of a 

k=size(beta0,1);

prior_th=a*(-k*log(2)+k*log(b)/2-k*log(c)/2+0.5*sum(log(lambda))+0.5*sum(log(beta0.^2)))...
    +1.25*k*log(a)+k*a*log(a)/2-k*log(gamma(a+1))...
    +sum(log(besselk(a-1/2,sqrt(lambda*b*a/c).*abs(beta0))));

prior_b=-log(beta(a,c))+a*(log(a)+log(b/(2*c)))-log(a)-(a+c)*log(1+a*b/(2*c));

prior_a=(hyp1-1)*log(2*a)+(hyp2-1)*log(1-2*a);

cv=log(a)+log(0.5-a);

loglike=prior_th+prior_b+prior_a+cv;

