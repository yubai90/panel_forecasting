function PAI=triang_new(Y,X,N,K,T,invA_,sqrt_ht,V)
% This is an improved triangular algorithm for the case when K is much
% greater than T, by Bhattacharya et al. 2016.
% Notice that gains may be small when K is almost equal to T and there are
% no gains when K<T.

PAI=zeros(K,N);
for j=1:N
    
    % iteratively compute previous equation residuals and rescale the data.
    % *** This is the sub-loop giving O(N^2) computational gains *** 
    Yr= Y(:,j);  
    if j==1; resid=zeros(T,N);
    else
        for l=1:j-1; Yr= Yr-invA_(j,l)*(sqrt_ht(:,l).*resid(:,l)); end
    end
    Yr= Yr./sqrt_ht(:,j); Xr= X./repmat(sqrt_ht(:,j),1,K);
    
    % index to select equation j coefficients from general specification
    index=K*(j-1)+1:(K*(j-1)+K);
    
    % perform the draw of PAI(equation j)|Pai(equations 1,...,j-1) 
    % following the algorithm proposed by Bhattacharya et al. 2016
    u=chol(V(index,index))*randn(K,1);   
    v=Xr*u+randn(T,1);
    w=(Xr*V(index,index)*Xr'+speye(T))\(Yr-v);
    pai_j=u+V(index,index)*Xr'*w;
   % V_post = (V(index,index) + Xr'*Xr)\eye(K); % equation 33
   % b_post = V_post*(iVb_prior(index) + Xr'*Yr); % equation 32
   % pai_j  = b_post + chol(V_post,'lower')*randn(K,1); % draw
    
    % The 2 lines below perform a draw from pai_j using a formula similar
    % to that in footnote 5, but equation by equation. These 2 lines are
    % based on inv(Chol_precision)*inv(Chol_precision)'=chol(V_post,'lower')*chol(V_post,'lower')'
    % Using these 2 lines will give further (but minor) speed improvements. 
    % Chol_precision = chol(iV(index,index) + Xr'*Xr);
    % pai_j = Chol_precision\((Chol_precision'\(iVb_prior(index) + Xr'*Yr)) + randn(K,1));
    
    % save residuals and present draw of pai_j
    resid(:,j)=Yr-Xr*pai_j; % store residuals to be removed in next equation
    PAI(:,j) = pai_j;       % store the draw of PAI  
    
end
            
