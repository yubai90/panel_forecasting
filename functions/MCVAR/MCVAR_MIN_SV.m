function Y_pred = MCVAR_MIN_SV(data,nsave,nburn,N,G,L,h,theta)
% Multi-country VAR-SV with Minnesota prior, we have one more
% hyperparameter to control the prior tightness for cross-country shrinkage

% create data matrices
[Traw,n]=size(data);
T=Traw-L;
X_lag=mlag2(data,L);
X=[ones(T,1) X_lag(L+1:Traw,:)];
Y=data(L+1:Traw,:);
K=n*L+1;        %  # of parameters in each equation 
m=n*(n-1)/2;    %  # of parameters in A
ntot=nsave+nburn;

% prior on conditional mean coefficient: use Minnesota setup, prior mean is
% set to 0 (variables in growth rates)
sig2=diag(panelarloop(data,1));
[C,idx_kappa1,idx_kappa2,idx_kappa3] = get_C_MC(N,G,L,sig2);

%N=7;G=3;
%OMEGA_pai=MIN_variance(N,G,L,data,theta); % This two lines are used as unrestricted VAR with additional shrinkage for other countries' variables


% prior mean and precision for A
a0=zeros(m,1);
iVa=0.1*eye(m);

% prior for PHI
d_PHI = n+2;                
s_PHI = 0.01*eye(n); 

% prior on initial states
Vol_0mean = zeros(n,1);   
Vol_0var  = 100*eye(n);   

% define something before MCMC loop  
Y_pred = zeros(nsave,n,h);


% initialize the Markov Chain
A_= eye(n);                                     
[~,~,ARresid]=OLS_PVAR(data,1,n,1);
htemp=mean(ARresid'.^2,2);
sqrt_ht=sqrt(repmat(htemp',T,1));       
Vol_states=2*log(sqrt_ht);                              
PHI_=0.0001*eye(n);    

comp=[eye(n*(L-1)),zeros(n*(L-1),n)]; 
check_stationarity=1;  

% | ----- MCMC starts here ----- |
for irep=1:ntot
    % step 1: Draw from the conditional posterior of PAI
    % This is the new step tp achieve O(N^2) computational gains
    Vbeta = getVthetaMC(idx_kappa1,idx_kappa2,idx_kappa3,theta,C,sig2);
    OMEGA_pai=sparse(1:n*K,1:n*K,Vbeta);
    iV=diag(1./diag(OMEGA_pai)); 
    
    if irep==1
        PAI=X\Y;
    end
    stationary=0;
    while stationary==0
        PAI=CTA(Y,X,n,K,A_,sqrt_ht,iV,zeros(K*n,1),PAI); 
        
        if (check_stationarity==0 || max(abs(eig([PAI(2:K,:)' ; comp]))) < 1); stationary = 1; end
    end 
      
    RESID = Y - X*PAI;
    beta=PAI(:);
    
    % Step 2: Draw theta1 and theta2
    tmpc1 = sum(beta(idx_kappa1).^2./C(idx_kappa1));
    tmpc2 = sum(beta(idx_kappa2).^2./C(idx_kappa2));
    tmpc3 = sum(beta(idx_kappa3).^2./C(idx_kappa3));
    theta(1) = gigrnd(1-n*L/2,2*(1/0.04),tmpc1,1);
    theta(2) = gigrnd(1-(G-1)*n*L/2,2*(1/0.04^2),tmpc2,1);   
    theta(4) = gigrnd(1-(n-G)*n*L/2,2*(1/0.02^2),tmpc3,1); 
    
    
    % step 3: Draw covariance
    count_E = 0;
    E = zeros(T*n,m);
    for jj=1:n-1
        E(jj+1:n:end,count_E+1:count_E+jj) = -RESID(:,1:jj);        
        count_E = count_E+jj;
    end
    iD = sparse(1:T*n,1:T*n,reshape(1./(sqrt_ht.^2)',1,T*n));
    Ka = iVa + E'*iD*E;
    a_hat = Ka\(iVa * a0 + E'*iD*reshape(RESID',T*n,1));   
    a = a_hat + chol(Ka,'lower')'\randn(m,1);  
    A_=chofac(n,a);
    invA_=A_\speye(n);
    
    % Step 4: Draw volatility state: KSC algorithm
    Vol_states  = KSC(log((RESID*A_').^2 + 1e-6),Vol_states,PHI_,Vol_0mean,Vol_0var);
    sqrt_ht  = exp(Vol_states/2); 
    
    % Step 5: Draw volatility variance
    eta  = Vol_states(2:end,:) - Vol_states(1:end-1,:); 
    temp = chol(inv(s_PHI + eta'*eta),'lower')*randn(n,T + d_PHI);
    PHI_ = (temp*temp')\speye(n);
    
    % simulate y path: use draws from posterior distribution
    if irep>nburn
        Y_pred(irep-nburn,:,:)=get_paths_VARSV(PAI(:),h,L,n,PHI_,invA_,Vol_states,Y);
    end
end

    