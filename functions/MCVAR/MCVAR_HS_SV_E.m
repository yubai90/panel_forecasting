function Y_pred =MCVAR_HS_SV_E(data,reps,burnin,N,G,L,h)
% We estimate the model using the triangular algorithm in Carriero, Clark
% and Marcellino (2019, JoE), together with Fast sampling from Gaussian
% distribution algorithm.
% Prior is Horseshoe grouped by equation.
% Triangular algorithm has been implemented as in the correct version.

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


% Initialize hyperparameters: here is simple as we group parameters
% equation by equation

Omega=0.0001*ones(K_eq,NG);Gamma=ones(K_eq,NG);
Omega_a=ones(m,1);Gamma_a=ones(m,1);
lambda=ones(NG+1,1);c=ones(NG+1,1);
Sigma_b= Omega;
%V=sparse(1:K_eq*NG,1:K_eq*NG,Sigma_b(:)); 
iV=sparse(1:K_eq*NG,1:K_eq*NG,1./Sigma_b(:)); 
iVb_prior=zeros(K_eq*NG,1);
Sigma_a=Omega_a;     % prior variance for a

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

% prior for PHI and A
d_PHI = NG+2;                
s_PHI = 0.01*eye(NG); 
A_= eye(NG);                                     
invA_= inv(A_); 

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


    
    % Step 3: Update the hyperparameters equation by equation
    % For intercept and auto regressive coefficients
    for eq=1:NG
        
        % Update Omega and Gamma
        for ie=1:K_eq
            Omega(ie,eq)=gigrnd(0,2*Gamma(ie,eq),PAI(ie,eq)^2,1);
            Gamma(ie,eq)=gamrnd(1,1/(lambda(eq)+Omega(ie,eq)));
        end
        
        % Update lambda and c
        lambda(eq)=gamrnd(1+(K_eq-1)/2,1/(c(eq)+sum(Gamma(:,eq))));
        c(eq)=gamrnd(1,1/(lambda(eq)+1));
        
    end
    
    % For the free elements in A
    for ia=1:m
        Omega_a(ia)=gigrnd(0,2*Gamma_a(ia),alp(ia)^2,1);
    end
   
    lambda(end)=gamrnd(1+(m-1)/2,1/(c(end)+sum(Gamma_a)));
    c(end)=gamrnd(1,1/(lambda(end)+1));
    
   
    Sigma_b= Omega;
    Sigma_a=Omega_a; 
    Sigma_b(Sigma_b<1e-10)=1e-10;
   % V=sparse(1:K_eq*NG,1:K_eq*NG,Sigma_b(:)); 
    iV=sparse(1:K_eq*NG,1:K_eq*NG,1./Sigma_b(:)); 
    Sigma_a(Sigma_a<1e-10)=1e-10;
     
    % Step 7: Draw volatility state: KSC algorithm with 10 normal mixtures
    Vol_states  = KSC(log((RESID*A_').^2 + 1e-6),Vol_states,PHI_,Vol_0mean,Vol_0var);
    sqrt_ht  = exp(Vol_states/2); 
    
    % Step 8: Draw volatility variance
    eta  = Vol_states(2:end,:) - Vol_states(1:end-1,:); 
    temp = chol(inv(s_PHI + eta'*eta),'lower')*randn(NG,T + d_PHI);
    PHI_ = (temp*temp')\speye(NG);
      
   % simulate y path: use draws from posterior distribution
    if irep>burnin
        Y_pred(irep-burnin,:,:)=get_paths_VARSV(PAI(:),h,L,NG,PHI_,invA_,Vol_states,Y);
    end
           
end