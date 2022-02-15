function  [alpha_OLS_vec, sigma_OLS,err] =OLS_PVAR(Yraw,N,G,p)

% Create VAR data matrices
[Traw, NG] = size(Yraw);
if NG ~= N*G; error('wrong specification of N and G'); end  % Check dimensions
Ylag = mlag2(Yraw,p);
k = p*NG;             % number of coefficients in each equation
T=Traw-p;
X = Ylag(p+1:Traw,:); % VAR data (RHS) matrix on the original specification    
x = kron(eye(NG),X);  % VAR data (RHS) matrix on the SURE model
% Correct time-series dimension of Y due to taking lags
Y = Yraw(p+1:Traw,:); 
y = Y(:);  

% OLS results
alpha_OLS_vec = (x'*x)\(x'*y);
alpha_OLS_mat = (X'*X)\(X'*Y);
err=Y - X*alpha_OLS_mat;
SSE = err'*err;
sigma_OLS = SSE./(T-(k-1));
