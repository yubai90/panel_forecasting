function [laglength, IC, logL] = VARlag(ENDO,maxlag,const,flag)
% =======================================================================
% Determine VAR lag length wit Information criterion
% const: 0 no const; 1 const ; 2 const&trend;
% flag: 1 for AIC criterion, 2 for BIC criterion

% OUTPUT
%	- laglength: preferred lag lenghth according to Akaike information criterion
%   - IC: vector [maxlag x 1] of (Akaike or Bayesian) information criterion
%   - logL: vector [maxlag x 1] of loglikelihood
% =======================================================================

% Note: 
% (1) The determinant of the residual covariance is computed without adjusting for the degrees of freedom, as in Eviews. 
% (2) The log likelihood value is computed assuming a multivariate normal
% (Gaussian) distribution.

%% Check inputs

% Check if ther are constant, trend, both, or none
if ~exist('const','var')
    const = 1;
end

%% Compute log likelihood and Bayesian criterion
logL = zeros(maxlag,1);
BIC  = zeros(maxlag,1);
AIC  = zeros(maxlag,1);
for i=1:maxlag
    X = ENDO(maxlag+1-i:end,:);
    aux = VARmodel(X,i,const);
    NOBS = aux.nobs;
    NEQS = aux.nvar;
    NTOTCOEFF = aux.ncoeff;
    RES = aux.residuals;
    SIGMA = (1/(NOBS)).*(RES)'*(RES); 
    logL(i) = -(NOBS/2)* (NEQS*(1+log(2*pi)) + log(det(SIGMA)));
    BIC(i)  = -2*(logL(i)/NOBS) + 2*(NTOTCOEFF*log(NOBS)/NOBS);
    AIC(i) = -2*(logL(i)/NOBS) + 2*(NTOTCOEFF/NOBS);
end
   if flag==1
       laglength = find(AIC==min(AIC));
       IC=AIC;
   elseif flag==2
      laglength = find(BIC==min(BIC));
      IC=BIC;
   end


