function arvar=panelarloop(data_endo,p)

% preliminaries
[T,G]=size(data_endo);
arvar=zeros(G,G);

% obtain the (pooled) variance for each endogenous variables

% loop over endogenous variables
for ii=1:G
   
    % define x and y
    y=data_endo(:,ii);
    xlag=mlag2(data_endo(:,ii),p);
    xreg=[ones(T-p,1) xlag(p+1:end,1)];
    yreg=y(p+1:end,1);
    
    % obtain the OLS estimator
    B=(xreg'*xreg)\(xreg'*yreg);

    % obtain the vector of residuals
    eps=yreg-xreg*B;

    % obtain the variance of the series;
    arvar(ii,ii)=(eps'*eps)/(T-1-2);

end
