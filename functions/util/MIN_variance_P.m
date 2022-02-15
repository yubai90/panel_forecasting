function omegab=MIN_variance_P(N,G,p,data_endo,theta)

% obtain omegab, a Minnesota covariance matrix, for unrestricted VAR, imposing additional shrinakge 
% on other countries' variables
% theta(1): the overall shrinkage hyperparameter
% theta(2): cross variable shrinkage
% theta(3): intercept
% theta(4): lagged shrinkage

NG=N*G;                % total number of variables
theta_p=0.5;           % overall shrinkage parameters for other countries' variables

% first obtain the residual variance of individual (pooled) autoregressive models
arvar=panelarloop(data_endo,1);

K=NG*p+1;  % # of parameters in each equation

Pi_pv=zeros(NG*(K-1),1); co=0; sigma_const=zeros(NG,1);
omega=[];
for i=1:NG
    ind=ceil(i/G);                      % country index
    sigma_const(i)=arvar(i,i)*theta(3); % this sets the prior variance on the intercept 
    for l=1:p
        for j=1:NG
            co=co+1;
            if (i==j)
                Pi_pv(co)=theta(1)/(l^theta(4)); % prior variance, own lags
            elseif (ind-1)*G+1<=j<=ind*G
                Pi_pv(co)=(arvar(i,i)/arvar(j,j)*theta(1)*theta(2)/(l^theta(4))); % prior variance, cross-lags
            else
                Pi_pv(co)=(arvar(i,i)/arvar(j,j)*theta(1)*theta_p*theta(2)/(l^theta(4))); % impose additional shrinkage for other countries
            end
        end
    end
    omega=[omega;sigma_const(i);Pi_pv((K-1)*(i-1)+1:((K-1)*(i-1)+K-1))];
end

omegab= sparse(1:NG*K,1:NG*K,omega);


