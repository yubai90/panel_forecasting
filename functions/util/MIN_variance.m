function omegab=MIN_variance(G,p,data_endo,theta)

% obtain omegab, a Minnesota covariance matrix 
% theta(1): the overall shrinkage hyperparameter
% theta(2): cross variable shrinkage
% theta(3): intercept
% theta(4): lagged shrinkage

% first obtain the residual variance of individual (pooled) autoregressive models
arvar=panelarloop(data_endo,1);

K=G*p+1;  % # of parameters in each equation

Pi_pv=zeros(G*(K-1),1); co=0; sigma_const=zeros(G,1);
omega=[];
for i=1:G
    sigma_const(i)=arvar(i,i)*theta(3); % this sets the prior variance on the intercept 
    for l=1:p
        for j=1:G
            co=co+1;
            if (i==j)
                Pi_pv(co)=theta(1)/(l^theta(4)); % prior variance, own lags
            else
                Pi_pv(co)=(arvar(i,i)/arvar(j,j)*theta(1)*theta(2)/(l^theta(4))); % prior variance, cross-lags
            end
        end
    end
    omega=[omega;sigma_const(i);Pi_pv((K-1)*(i-1)+1:((K-1)*(i-1)+K-1))];
end

omegab= sparse(1:G*K,1:G*K,omega);


