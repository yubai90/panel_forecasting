% This function constructs the conditional prior of the VAR coefficients
% and the impact matrix for multi-country VARs with Minnesota type prior

function Vbeta = getVthetaMC(idx_kappa1,idx_kappa2,idx_kappa3,kappa,C,sig2)
np = length(idx_kappa1);
k_beta = length(C);

Vbeta = zeros(k_beta,1);


Vbeta(1:np+1:end) = kappa(3)*sig2;          % intercepts
Vbeta(idx_kappa1) = kappa(1)*C(idx_kappa1); % own lags
Vbeta(idx_kappa2) = kappa(2)*C(idx_kappa2); % other lags in the country
Vbeta(idx_kappa3) = kappa(4)*C(idx_kappa3); % other lags in the country

