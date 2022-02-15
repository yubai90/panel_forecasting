function loglike = like_c(a,b,c,tau,para,hyp1,hyp2)
% compute the log-likelihood function of the posterior of a 

k=size(para,1);

prior_pa=k*log(gamma(c+0.5))-k*log(gamma(c+1))+k*log(c)/2-...
    (c+0.5)*(sum(log(4*c*tau+para.^2*b*a))-sum(log(4*c*tau)));

prior_b=-log(beta(a,c))-(a-1)*log(c)-(a+c)*log(1+a*b/(2*c));

prior_c=(hyp1-1)*log(2*c)+(hyp2-1)*log(1-2*c);

cv=log(c)+log(0.5-c);

loglike=prior_pa+prior_b+prior_c+cv;