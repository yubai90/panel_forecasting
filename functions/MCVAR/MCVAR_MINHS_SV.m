function Y_pred = MCVAR_MINHS_SV(data,nsave,nburn,N,G,L,h)
% Multi-country VAR-SV with Minnesota type Horseshoe prior, we have one more
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
theta=[0.003 0.0015 100 0.001]; 
xi=10*ones(4,1);xia=10;thetaa=0.01;
cf=[1;4;9;16];  %  quadratic lag shrinkage factor


% prior on conditional mean coefficient: use Minnesota setup, prior mean is
% set to 0 (variables in growth rates), initialized at convential Minnesota
% prior
sig2=0.01*ones(n+1,n);gam2=0.01*ones(n+1,n);
sig2a=0.01*ones(m,1);gam2a=0.01*ones(m,1);
[C,idx_kappa1,idx_kappa2,idx_kappa3] = get_C_MC_H(N,G,L,sig2);
Vbeta = getVthetaMC_H(idx_kappa1,idx_kappa2,idx_kappa3,theta,C);
idx_K1=reshape(idx_kappa1,L,n);
idx_K2=reshape(idx_kappa2,L*(G-1),n);
idx_K3=reshape(idx_kappa3,G*L*(N-1),n);
OMEGA_pai=sparse(1:n*K,1:n*K,Vbeta);
iV=diag(1./diag(OMEGA_pai)); 

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
    
           
    % step 2: Draw covariance
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

    
    % step 3: Draw local variance parameters
    % intercept
    for j=1:n
        sig2(1,j)=gigrnd(0,2*gam2(1,j),PAI(1,j)^2/theta(3),1);
        gam2(1,j)=gamrnd(1,1/(sig2(1,j)+1));
    end
    % own lags
    for j=1:n
        id_temp=idx_K1(:,j);
        sig2(j+1,j)=gigrnd(0.5-0.5*L,2*gam2(j+1,j),sum(beta(id_temp).^2.*cf)/theta(1),1);
        gam2(j+1,j)=gamrnd(1,1/(sig2(j+1,j)+1));
    end
    % cross-variable lags: own country
    for j=1:n
        inc=ceil(j/G);           % country index
        idxg=(inc-1)*G+2:inc*G+1;
        id_temp2=idx_K2(:,j);
        if rem(j,G)==1
            idxg2=idxg(2:3);
        elseif rem(j,G)==2
            idxg2=idxg([1 3]);
        elseif rem(j,G)==0
            idxg2=idxg(1:2);
        end
        for k=1:G-1
            id_temp21=id_temp2(k:G-1:end);
            sig2(idxg2(k),j)=gigrnd(0.5-0.5*L,2*gam2(idxg2(k),j),sum(beta(id_temp21).^2.*cf)/theta(2),1);
            gam2(idxg2(k),j)=gamrnd(1,1/(sig2(idxg2(k),j)+1));
        end
    end
    % cross-variable lags: foreign country
    for j=1:n
        inc=ceil(j/G); 
        idxtem=1:n+1;
        idxg=[1 (inc-1)*G+2:inc*G+1];
        idxtem(idxg)=[];
        id_temp3=idx_K3(:,j);
        for k=1:G*(N-1)
            id_temp31=id_temp3(k:G*(N-1):end);
            sig2(idxtem(k),j)=gigrnd(0.5-0.5*L,2*gam2(idxtem(k),j),sum(beta(id_temp31).^2.*cf)/theta(4),1);
            gam2(idxtem(k),j)=gamrnd(1,1/(sig2(idxtem(k),j)+1));
        end
    end
    % free elements in A
    for j=1:m
        sig2a(j)=gigrnd(0,2*gam2a(j),a(j)^2/thetaa,1);
        gam2a(j)=gamrnd(1,1/(sig2a(j)+1));
    end
    
    % Step 4: Draw thetas
    C = get_C_MC_H(N,G,L,sig2);
    tmpc1 = sum(beta(idx_kappa1).^2./C(idx_kappa1));
    tmpc2 = sum(beta(idx_kappa2).^2./C(idx_kappa2));
    tmpc3 = sum(beta(idx_kappa3).^2./C(idx_kappa3));
    
    theta(1) = gigrnd(0.5-n*L/2,2*xi(1),tmpc1,1);
    xi(1)=gamrnd(1,1/(theta(1)+1));
    
    theta(2) = gigrnd(0.5-(G-1)*n*L/2,2*xi(2),tmpc2,1); 
    xi(2)=gamrnd(1,1/(theta(2)+1));
    
    theta(4) = gigrnd(0.5-(n-G)*n*L/2,2*xi(4),tmpc3,1); 
    xi(4)=gamrnd(1,1/(theta(4)+1));
    % Here, we also updtate intercept
    tmpc4= sum(PAI(1,:).^2./sig2(1,:));
    theta(3)= gigrnd(0.5-n/2,2*xi(3),tmpc4,1);
    xi(3)=gamrnd(1,1/(theta(3)+1));
    % We also need to update global shrinkage parameters for a
    tmpca= sum(a.^2./sig2a);
    thetaa= gigrnd(0.5-m/2,2*xia,tmpca,1);
    xia=gamrnd(1,1/(thetaa+1));
    
    % Update the variance parameters
    Vbeta = getVthetaMC_H(idx_kappa1,idx_kappa2,idx_kappa3,theta,C);
    OMEGA_pai=sparse(1:n*K,1:n*K,Vbeta);
    iV=diag(1./diag(OMEGA_pai)); 
    iVa=sparse(1:m,1:m,1./sig2a);
    
    
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
