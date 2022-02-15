function Y_pred =MCVAR_NGG_SV(data,reps,burnin,N,G,L,h)
% We estimate the model using the triangular algorithm in Carriero, Clark
% and Marcellino (2019, JoE) with corrections.

% compute and define some preliminaries 
[Traw,NG]=size(data);
T=Traw-L;
K_eq=NG*L+1;

nm1=NG*L;nm2=(G-1)*L*NG;   % # of parameters related to own-lags and cross-variable lags
%nm=G*L+1;                 % # of parameters in each country (for each equation, CSH)
nd=G*L*(N-1);              % # of parameters in each country (for each equation, DI)
m = NG*(NG-1)/2;           % # of free elements in A
alp=zeros(m,1);
fsize=reps-burnin;
Y_pred = zeros(fsize,NG,h);

% Construct data matrix
Y=data(L+1:Traw,:);
Y_lag=mlag2(data,L);
X=[ones(T,1) Y_lag(L+1:Traw,:)];

% index used to select coefficients
index_dom=zeros(N,G*L+1);
index_for=zeros(N,(NG-G)*L);
for i=1:N
    for j=1:L
        index=(j-1)*NG+2:j*NG+1;
        index_dom(i,(j-1)*G+2:j*G+1)=(i-1)*G+(j-1)*NG+2:i*G+(j-1)*NG+1;
        index((i-1)*G+1:i*G)=[];
        index_for(i,(j-1)*(NG-G)+2:j*(NG-G)+1)=index;
    end
end

index_dom(:,1)=[];
index_for=index_for(:,2:end);

% Initialize hyperparameters: notet that for local variance parameters we
% only need to define 2 K*NG matrices and we will get grouped parameters
% from the function get_grhyp.m. We also need to define separatelt for the
% free elements in A.
Omega=0.0001*ones(K_eq,NG);Gamma=ones(K_eq,NG);
Omega_a=ones(m,1);Gamma_a=ones(m,1);
lambda=0.1*ones(N+4,1);c=0.5*ones(N+4,1);kappa=10*ones(N+4,1);d=ones(N+4,1);
lambda_loc=0.1*ones(N+4,1);c_loc=0.5*ones(N+4,1);
phi=2*c./(lambda.*kappa);
Sigma_b= get_priorV(Omega,Gamma,phi(1:end-1),index_dom,index_for,N,G,L);
%V=sparse(1:K_eq*NG,1:K_eq*NG,Sigma_b(:)); 
iV=sparse(1:K_eq*NG,1:K_eq*NG,1./Sigma_b(:)); 
Sigma_a=phi(end)*Omega_a./Gamma_a;     % prior variance for a

% prior for lambda and c
alpha_h=2;
beta_h=1;

% prior on initial states for volatilities
Vol_0mean = zeros(NG,1);   
Vol_0var  = 100*eye(NG); 

% initialize the volatility
ARresid=zeros(Traw-1,NG);
for i=1:NG
    [~,~,ARresid(:,i)]=OLS_PVAR(data(:,i),1,1,1);
end
htemp=mean(ARresid'.^2,2);
sqrt_ht=sqrt(repmat(htemp',T,1));    
Vol_states=2*log(sqrt_ht);                              
PHI_=0.0001*eye(NG);  

% prior for PHI
d_PHI = NG+2;                
s_PHI = 0.01*eye(NG); 

% Initialize for Adaptive RWM 
scaleg=ones(N+4,1);
z=10*ones(N+4,1);z_loc=10*ones(N+4,1);
scaleg2=ones(N+4,1);
z2=10*ones(N+4,1);z_loc2=10*ones(N+4,1);
alpha=zeros(N+4,1);
alpha2=zeros(N+4,1);

A_= eye(NG);                                     

comp=[eye(NG*(L-1)),zeros(NG*(L-1),NG)]; 


%  ========= | MCMC starts here |=============
for irep=1:reps
    
    % step 1: Draw VAR mean coefficients
   if irep==1
        PAI=X\Y;
    end
    if irep<=burnin
        PAI=CTA(Y,X,NG,K_eq,A_,sqrt_ht,iV,zeros(K_eq*NG,1),PAI); 
    else 
        stationary=0;
        while stationary==0
            
            % This is the only new step (triangular algorithm).
            PAI=CTA(Y,X,NG,K_eq,A_,sqrt_ht,iV,zeros(K_eq*NG,1),PAI); 
            
            if max(abs(eig([PAI(2:K_eq,:)' ; comp]))) < 1 
                stationary = 1; 
            end
        end 
    end
   % PAI=triang(Y,X,NG,K_eq,T,invA_,sqrt_ht,iV,iVb_prior); 
    RESID = Y - X*PAI;
    
    % step 2: Draw free elements in A
    count = 0;
    ivaprior=sparse(1:m,1:m,1./Sigma_a);
    for ii = 2:NG
       
        % weighted regression to get Z'Z and Z'z (in Cogley-Sargent 2005 notation)
        y_spread_adj=RESID(:,ii)./sqrt_ht(:,ii);
        X_spread_adj=[]; for vv=1:ii-1;  X_spread_adj=[X_spread_adj RESID(:,vv)./sqrt_ht(:,ii)]; end  %#ok<AGROW>
        ZZ=X_spread_adj'*X_spread_adj; Zz=X_spread_adj'*y_spread_adj;

        % computing posteriors moments
        Valpha_post = (ZZ + ivaprior(count+1:count+ii-1,count+1:count+ii-1))\speye(ii-1);
        alpha_post  = Valpha_post*Zz;
        % draw and store 
        alphadraw   = alpha_post+chol(Valpha_post,'lower')*randn(ii-1,1);
        a1=-1*alphadraw;
        A_(ii,1:ii-1)= a1';
    
        alp(count+1:count+ii-1)=a1;
        count = count + ii-1;
        
    end
    invA_= inv(A_); 

   
    % Before updating other parameters, let us first get blocked parameters
    [Omega_int,Omega_olag,Omega_cvlag,Omega_di] = get_grhyp(Omega,index_dom,index_for,N,G,L);
    [Gamma_int,Gamma_olag,Gamma_cvlag,Gamma_di] = get_grhyp(Gamma,index_dom,index_for,N,G,L);
    [Beta_int,Beta_olag,Beta_cvlag,Beta_di] = get_grhyp(PAI,index_dom,index_for,N,G,L);

    
    % Step 3: Draw lambda with adaptive MH algorithm
    z_loc(1)=z(1)+sqrt(scaleg(1))*randn;
    lambda_loc(1)=0.5*exp(z_loc(1))/(1+exp(z_loc(1)));
    alpha(1)=min(exp(like_a(lambda_loc(1),kappa(1),c(1),Gamma_int',Beta_int',alpha_h,beta_h)-...
        like_a(lambda(1),kappa(1),c(1),Gamma_int',Beta_int',alpha_h,beta_h)),1);
    if alpha(1)>rand
        lambda(1)=lambda_loc(1);
        z(1)=log(lambda(1)/(0.5-lambda(1)));
    end
    
    z_loc(2)=z(2)+sqrt(scaleg(2))*randn;
    lambda_loc(2)=0.5*exp(z_loc(2))/(1+exp(z_loc(2)));
    alpha(2)=min(exp(like_a(lambda_loc(2),kappa(2),c(2),Gamma_olag,Beta_olag,alpha_h,beta_h)-...
        like_a(lambda(2),kappa(2),c(2),Gamma_olag,Beta_olag,alpha_h,beta_h)),1);
    if alpha(2)>rand
        lambda(2)=lambda_loc(2);
        z(2)=log(lambda(2)/(0.5-lambda(2)));
    end
    
    z_loc(3)=z(3)+sqrt(scaleg(3))*randn;
    lambda_loc(3)=0.5*exp(z_loc(3))/(1+exp(z_loc(3)));
    alpha(3)=min(exp(like_a(lambda_loc(3),kappa(3),c(3),Gamma_cvlag,Beta_cvlag,alpha_h,beta_h)-...
        like_a(lambda(3),kappa(3),c(3),Gamma_cvlag,Beta_cvlag,alpha_h,beta_h)),1);
    if alpha(3)>rand
        lambda(3)=lambda_loc(3);
        z(3)=log(lambda(3)/(0.5-lambda(3)));
    end
    
    for jj=1:N
        z_loc(3+jj)=z(3+jj)+sqrt(scaleg(3+jj))*randn;
        lambda_loc(3+jj)=0.5*exp(z_loc(3+jj))/(1+exp(z_loc(3+jj)));
        alpha(3+jj)=min(exp(like_a(lambda_loc(3+jj),kappa(3+jj),c(3+jj),Gamma_di(:,jj),Beta_di(:,jj),alpha_h,beta_h)-...
        like_a(lambda(3+jj),kappa(3+jj),c(3+jj),Gamma_di(:,jj),Beta_di(:,jj),alpha_h,beta_h)),1);
        if alpha(3+jj)>rand
            lambda(3+jj)=lambda_loc(3+jj);
            z(3+jj)=log(lambda(3+jj)/(0.5-lambda(3+jj)));
        end
    end
    
    z_loc(end)=z(end)+sqrt(scaleg(end))*randn;
    lambda_loc(end)=0.5*exp(z_loc(end))/(1+exp(z_loc(end)));
    alpha(end)=min(exp(like_a(lambda_loc(end),kappa(end),c(end),Gamma_a,alp,alpha_h,beta_h)-...
        like_a(lambda(end),kappa(end),c(end),Gamma_a,alp,alpha_h,beta_h)),1);
    if alpha(end)>rand
        lambda(end)=lambda_loc(end);
        z(end)=log(lambda(end)/(0.5-lambda(end)));
    end
    phi=2*c./(lambda.*kappa);
    
    % Step 4: Draw Omega from GIG distribution
    for ie=1:NG
        Omega_int(ie)=gigrnd(lambda(1)-0.5,2,Gamma_int(ie)*Beta_int(ie)^2/phi(1),1);
    end
    for ie=1:nm1
        Omega_olag(ie)=gigrnd(lambda(2)-0.5,2,Gamma_olag(ie)*Beta_olag(ie)^2/phi(2),1);
    end
    for ie=1:nm2
        Omega_cvlag(ie)=gigrnd(lambda(3)-0.5,2,Gamma_cvlag(ie)*Beta_cvlag(ie)^2/phi(3),1);
    end
    for jt=1:N
        for ie=1:nd*G
            Omega_di(ie,jt)=gigrnd(lambda(jt+3)-0.5,2,Gamma_di(ie,jt)*Beta_di(ie,jt)^2/phi(jt+3),1);
        end
    end
    for ie=1:m
        Omega_a(ie)=gigrnd(lambda(end)-0.5,2,Gamma_a(ie)*alp(ie)^2/phi(end),1);
    end
    
    % Step 5: Update c with Adaptive MH algorithm
    z_loc2(1)=z2(1)+sqrt(scaleg2(1))*randn;
    c_loc(1)=0.5*exp(z_loc2(1))/(1+exp(z_loc2(1)));
    alpha2(1)=min(exp(like_c(lambda(1),kappa(1),c_loc(1),Omega_int',Beta_int',alpha_h,beta_h)-...
        like_c(lambda(1),kappa(1),c(1),Omega_int',Beta_int',alpha_h,beta_h)),1);
    if alpha2(1)>rand
        c(1)=c_loc(1);
        z2(1)=log(c(1)/(0.5-c(1)));
    end
    
    z_loc2(2)=z2(2)+sqrt(scaleg2(2))*randn;
    c_loc(2)=0.5*exp(z_loc2(2))/(1+exp(z_loc2(2)));
    alpha2(2)=min(exp(like_c(lambda(2),kappa(2),c_loc(2),Omega_olag,Beta_olag,alpha_h,beta_h)-...
        like_c(lambda(2),kappa(2),c(2),Omega_olag,Beta_olag,alpha_h,beta_h)),1);
    if alpha2(2)>rand
        c(2)=c_loc(2);
        z2(2)=log(c(2)/(0.5-c(2)));
    end
    
    z_loc2(3)=z2(3)+sqrt(scaleg2(3))*randn;
    c_loc(3)=0.5*exp(z_loc2(3))/(1+exp(z_loc2(3)));
    alpha2(3)=min(exp(like_c(lambda(3),kappa(3),c_loc(3),Omega_cvlag,Beta_cvlag,alpha_h,beta_h)-...
        like_c(lambda(3),kappa(3),c(3),Omega_cvlag,Beta_cvlag,alpha_h,beta_h)),1);
    if alpha2(3)>rand
        c(3)=c_loc(3);
        z2(3)=log(c(3)/(0.5-c(3)));
    end
    
    for jj=1:N
        z_loc2(3+jj)=z2(3+jj)+sqrt(scaleg2(3+jj))*randn;
        c_loc(3+jj)=0.5*exp(z_loc2(3+jj))/(1+exp(z_loc2(3+jj)));
        alpha2(3+jj)=min(exp(like_c(lambda(3+jj),kappa(3+jj),c_loc(3+jj),Omega_di(:,jj),Beta_di(:,jj),alpha_h,beta_h)-...
            like_c(lambda(3+jj),kappa(3+jj),c(3+jj),Omega_di(:,jj),Beta_di(:,jj),alpha_h,beta_h)),1);
        if alpha2(3+jj)>rand
            c(3+jj)=c_loc(3+jj);
            z2(3+jj)=log(c(3+jj)/(0.5-c(3+jj)));
        end
    end
    
    z_loc2(end)=z2(end)+sqrt(scaleg2(end))*randn;
    c_loc(end)=0.5*exp(z_loc2(end))/(1+exp(z_loc2(end)));
    alpha2(end)=min(exp(like_c(lambda(end),kappa(end),c_loc(end),Omega_a,alp,alpha_h,beta_h)-...
        like_c(lambda(end),kappa(end),c(end),Omega_a,alp,alpha_h,beta_h)),1);
    if alpha2(end)>rand
        c(end)=c_loc(end);
        z2(end)=log(c(end)/(0.5-c(end)));
    end
    
    phi=2*c./(lambda.*kappa);
    
    % Update the variance of the increments
    if irep>100
        for iia=1:N+4
            scaleg(iia)=scaleg(iia)*exp(irep^(-0.55)*(alpha(iia)-0.3));
            scaleg2(iia)=scaleg2(iia)*exp(irep^(-0.55)*(alpha2(iia)-0.3));
        end
    end
        
    % Step 6: Draw Gamma from Ga distribution
    for ie=1:NG
        Gamma_int(ie)=gamrnd(0.5+c(1),1/(1+Beta_int(ie)^2/(2*phi(1)*Omega_int(ie))));
    end
    for ie=1:nm1
        Gamma_olag(ie)=gamrnd(0.5+c(2),1/(1+Beta_olag(ie)^2/(2*phi(2)*Omega_olag(ie))));
    end
    for ie=1:nm2
        Gamma_cvlag(ie)=gamrnd(0.5+c(3),1/(1+Beta_cvlag(ie)^2/(2*phi(3)*Omega_cvlag(ie))));
    end
    for jt=1:N
        for ie=1:nd*G
            Gamma_di(ie,jt)=gamrnd(0.5+c(3+jt),1/(1+Beta_di(ie,jt)^2/(2*phi(3+jt)*Omega_di(ie,jt))));
        end
     end
    for ie=1:m
        Gamma_a(ie)=gamrnd(0.5+c(end),1/(1+alp(ie)^2/(2*phi(end)*Omega_a(ie))));
    end
    
    % Step 7: Draw kappa from hiearchical representation
    d(1)=gamrnd(lambda(1)+c(1),1/(kappa(1)+2*c(1)/lambda(1)));
    kappa(1)=gamrnd(NG/2+lambda(1),1/(lambda(1)/(4*c(1))*sum(Gamma_int.*Beta_int.^2./Omega_int)+d(1)));
    d(2)=gamrnd(lambda(2)+c(2),1/(kappa(2)+2*c(2)/lambda(2)));
    kappa(2)=gamrnd(nm1/2+lambda(2),1/(lambda(2)/(4*c(2))*sum(Gamma_olag.*Beta_olag.^2./Omega_olag)+d(2)));
    d(3)=gamrnd(lambda(3)+c(3),1/(kappa(3)+2*c(3)/lambda(3)));
    kappa(3)=gamrnd(nm2/2+lambda(3),1/(lambda(3)/(4*c(3))*sum(Gamma_cvlag.*Beta_cvlag.^2./Omega_cvlag)+d(3)));
    for jj=1:N
        d(3+jj)=gamrnd(lambda(3+jj)+c(3+jj),1/(kappa(3+jj)+2*c(3+jj)/lambda(3+jj)));
        kappa(3+jj)=gamrnd(nd*G/2+lambda(3+jj),1/(lambda(3+jj)/(4*c(3+jj))*sum(Gamma_di(:,jj).*Beta_di(:,jj).^2./Omega_di(:,jj))+d(3+jj)));
    end
    d(end)=gamrnd(lambda(end)+c(end),1/(kappa(end)+2*c(end)/lambda(end)));
    kappa(end)=gamrnd(m/2+lambda(end),1/(lambda(end)/(4*c(end))*sum(Gamma_a.*alp.^2./Omega_a)+d(end)));

    phi=2*c./(lambda.*kappa);
    % We need to update prior variance, before this, we need to form
    % grouped local variance parameters back to K*NG matrix
    Omega = get_hypall(Omega_int,Omega_olag,Omega_cvlag,Omega_di,index_dom,index_for,N,G,L);
    Gamma = get_hypall(Gamma_int,Gamma_olag,Gamma_cvlag,Gamma_di,index_dom,index_for,N,G,L);
    Sigma_b= get_priorV(Omega,Gamma,phi(1:end-1),index_dom,index_for,N,G,L);
    Sigma_a=phi(end)*Omega_a./Gamma_a; 
    Sigma_b(Sigma_b<1e-10)=1e-10;
   % V=sparse(1:K_eq*NG,1:K_eq*NG,Sigma_b(:)); 
    iV=sparse(1:K_eq*NG,1:K_eq*NG,1./Sigma_b(:)); 
    Sigma_a(Sigma_a<1e-10)=1e-10;
     
    % Step 8: Draw volatility state: KSC algorithm with 10 normal mixtures
    Vol_states  = KSC(log((RESID*A_').^2 + 1e-6),Vol_states,PHI_,Vol_0mean,Vol_0var);
    sqrt_ht  = exp(Vol_states/2); 
    
    % Step 9: Draw volatility variance
    eta  = Vol_states(2:end,:) - Vol_states(1:end-1,:); 
    temp = chol(inv(s_PHI + eta'*eta),'lower')*randn(NG,T + d_PHI);
    PHI_ = (temp*temp')\speye(NG);
      
   % simulate y path: use draws from posterior distribution
    if irep>burnin
        Y_pred(irep-burnin,:,:)=get_paths_VARSV(PAI(:),h,L,NG,PHI_,invA_,Vol_states,Y);
    end
           
end