function VAR = VARmodel(ENDO,nlag,const)
% =======================================================================
% Perform vector autogressive (VAR) estimation with OLS 
% const: 0 no constant; 1 constant; 2 constant and trend;
% OUTPUT
%   - VAR: structure including VAR estimation results
%   - VARopt: structure including VAR options (see VARoption)
% =======================================================================
% Written by: Yu Bai
% April 2017

%% Check inputs
[nobs, nvar] = size(ENDO);

% Check if ther are constant, trend, both, or none
if ~exist('const','var')
    const = 1;
end

%% create data matrices
nobse         = nobs - nlag;
VAR.nlag          = nlag;
VAR.nobs      = nobse;
VAR.nvar      = nvar;
ncoeff        = nvar*nlag+const; 
VAR.ncoeff = ncoeff;

% Create independent vector and lagged dependent matrix
Y = ENDO(nlag+1:end,:);
if const==0
    X=[];
    for j=0:nlag-1
        X = [ENDO(j+1:nobs-nlag+j,:), X];
    end 
elseif const==1
    X=[];
     for j=0:nlag-1
        X = [ENDO(j+1:nobs-nlag+j,:), X];
     end
    X=[ones(nobs-nlag,1) X];
elseif const==2
    X=[];
     for j=0:nlag-1
        X = [ENDO(j+1:nobs-nlag+j,:), X];
     end
     trend=1:size(X,1);
    X=[ones(nobs-nlag,1) trend' X];
end
        
%% Compute the matrix of coefficients & VCV
Ft = (X'*X)\(X'*Y);
VAR.Ft = Ft;
SIGMA = (1/(nobse-ncoeff))*(Y-X*Ft)'*(Y-X*Ft); 
VAR.sigma = SIGMA;
VAR.residuals = Y - X*Ft;

%% companion form of the coefficient matrix
F = Ft';
VAR.F=F;
Fcomp = [F(:,1+const:nvar*nlag+const); eye(nvar*(nlag-1)) zeros(nvar*(nlag-1),nvar)];
VAR.Fcomp = Fcomp;

