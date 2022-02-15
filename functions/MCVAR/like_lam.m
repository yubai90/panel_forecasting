function loglike=like_lam(b,lambda,kappa,omega)

% Compute log-likelihood of posterior of lambda, by conditioning on the
% local variance parameters

p=size(omega,1);
loglike=-b*lambda+p*lambda*log(lambda*kappa/2)-p*log(gamma(lambda))...
    +(lambda-1)*sum(log(omega))-lambda*kappa*sum(omega)/2;