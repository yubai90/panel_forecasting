function Y_pred =MCVAR_HS(data,reps,burnin,N,G,L,h)
% We estimate the model using the triangular algorithm in Carriero, Clark
% and Marcellino (2019, JoE), together with Fast sampling from Gaussian
% distribution algorithm.
% Prior is Horseshoe with hierarchical form. 
% Independent N-IW prior, no SV

% compute and define some preliminaries 
[Traw,NG]=size(data);
T=Traw-L;
K_eq=NG*L+1;

nm1=NG*L;nm2=(G-1)*L*NG;   % # of parameters related to own-lags and cross-variable lags
%nm=G*L+1;                 % # of parameters in each country (for each equation, CSH)
nd=G*L*(N-1);              % # of parameters in each country (for each equation, DI)
%m = NG*(NG-1)/2;           % # of free elements in A
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
[Gamma_int,Gamma_olag,Gamma_cvlag,Gamma_di] = get_grhyp(Gamma,index_dom,index_for,N,G,L);
Gamma_int=Gamma_int';
lambda=ones(N+3,1);c=ones(N+3,1);
Sigma_b= get_priorV2(Omega,index_dom,index_for,N,G,L);
%V=sparse(1:K_eq*NG,1:K_eq*NG,Sigma_b(:)); 
iV=sparse(1:K_eq*NG,1:K_eq*NG,1./Sigma_b(:)); 
iVb_prior=zeros(K_eq*NG,1);


% initialize sigma
sigma=diag(panelarloop(data,1));
scale_iw=diag(sigma);

% prior for PHI and A
A_= eye(NG);                                     
%invA_= inv(A_); 

comp=[eye(NG*(L-1)),zeros(NG*(L-1),NG)]; 
 

%  ========= | MCMC starts here |=============
for irep=1:reps
    
    % step 1: Draw VAR mean coefficients
    sigma_sqrt=repmat(sqrt(sigma'),T,1);
    if irep==1
        PAI=X\Y;
    end
    if irep<=burnin
        PAI=CTA(Y,X,NG,K_eq,A_,sigma_sqrt,iV,iVb_prior,PAI); 
    else 
        stationary=0;
        while stationary==0
            
            % This is the only new step (triangular algorithm).
            PAI=CTA(Y,X,NG,K_eq,A_,sigma_sqrt,iV,iVb_prior,PAI);
            
            if max(abs(eig([PAI(2:K_eq,:)' ; comp]))) < 1 
                stationary = 1; 
            end
        end 
    end
   % PAI=triang(Y,X,NG,K_eq,T,invA_,sqrt_ht,iV,iVb_prior); 
    RESID = Y - X*PAI;
    
    % step 2: Draw Sigma and decompose A and sigma
    SSE=RESID'*RESID;
    Sigma=iwishrnd(SSE+scale_iw,T+NG+2);
    
    [invA_,sigma]=ldlt(Sigma);
    A_=invA_\speye(NG);
    sigma=diag(sigma);

   
    % Before updating other parameters, let us first get blocked parameters
    [Omega_int,Omega_olag,Omega_cvlag,Omega_di] = get_grhyp(Omega,index_dom,index_for,N,G,L);
    [Beta_int,Beta_olag,Beta_cvlag,Beta_di] = get_grhyp(PAI,index_dom,index_for,N,G,L);

    
    % Step 3: Draw Omega from GIG distribution
    for ie=1:NG
        Omega_int(ie)=gigrnd(0,2*Gamma_int(ie),Beta_int(ie)^2,1);
    end
    for ie=1:nm1
        Omega_olag(ie)=gigrnd(0,2*Gamma_olag(ie),Beta_olag(ie)^2,1);
    end
    for ie=1:nm2
        Omega_cvlag(ie)=gigrnd(0,2*Gamma_cvlag(ie),Beta_cvlag(ie)^2,1);
    end
    for jt=1:N
        for ie=1:nd*G
            Omega_di(ie,jt)=gigrnd(0,2*Gamma_di(ie,jt),Beta_di(ie,jt)^2,1);
        end
    end

   
        
    % Step 4: Draw Gamma from Ga distribution
    for ie=1:NG
        Gamma_int(ie)=gamrnd(1,1/(lambda(1)+Omega_int(ie)));
    end
    for ie=1:nm1
        Gamma_olag(ie)=gamrnd(1,1/(lambda(2)+Omega_olag(ie)));
    end
    for ie=1:nm2
        Gamma_cvlag(ie)=gamrnd(1,1/(lambda(3)+Omega_cvlag(ie)));
    end
    for jt=1:N
        for ie=1:nd*G
            Gamma_di(ie,jt)=gamrnd(1,1/(lambda(3+jt)+Omega_di(ie,jt)));
        end
     end
    
    
    % Step 5 and 6: Draw lambda and c from hiearchical representation
    lambda(1)=gamrnd(1+(NG-1)/2,1/(c(1)+sum(Gamma_int)));
    c(1)=gamrnd(1,1/(lambda(1)+1));
    lambda(2)=gamrnd(1+(nm1-1)/2,1/(c(2)+sum(Gamma_olag)));
    c(2)=gamrnd(1,1/(lambda(2)+1));
    lambda(3)=gamrnd(1+(nm2-1)/2,1/(c(3)+sum(Gamma_cvlag)));
    c(3)=gamrnd(1,1/(lambda(3)+1));
    for jj=1:N
        lambda(3+jj)=gamrnd(1+(nd*G-1)/2,1/(c(3+jj)+sum(Gamma_di(:,jj))));
        c(3+jj)=gamrnd(1,1/(lambda(3+jj)+1));
    end
   
    
    % We need to update prior variance, before this, we need to form
    % grouped local variance parameters back to K*NG matrix
    Omega = get_hypall(Omega_int,Omega_olag,Omega_cvlag,Omega_di,index_dom,index_for,N,G,L);
    Sigma_b= get_priorV2(Omega,index_dom,index_for,N,G,L);
    Sigma_b(Sigma_b<1e-10)=1e-10;
   % V=sparse(1:K_eq*NG,1:K_eq*NG,Sigma_b(:)); 
    iV=sparse(1:K_eq*NG,1:K_eq*NG,1./Sigma_b(:)); 
   
  % simulate y path: use draws from posterior distribution
    if irep>burnin
        Y_pred(irep-burnin,:,:)=get_pathsVAR(NG,h,T,L,Y,PAI(:),Sigma);
    end
    
end