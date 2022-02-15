function loglike = like_lv(a,kappa,para)
% marginal prior for the local variance parameter

d=size(para,1);
t1=0.5*d*(a+0.5)*log(kappa*a);
t2=a*d*log(2);
t3=d*log(gamma(a));
t4=a*sum(log(abs(para)));
t5=sum(log(besselk(a-0.5,sqrt(kappa*a)*abs(para))));
loglike=t1-t2-t3+t4+t5;