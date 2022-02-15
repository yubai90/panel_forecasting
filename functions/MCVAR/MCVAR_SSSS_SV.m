function  Y_pred  = MCVAR_SSSS_SV(Yraw,N,G,p,nsave,burnin,ntot,h)
% predictive simulation for SSSS prior in Koop and Korobilis (2016) with SV

% Create VAR data matrices
[Traw, NG] = size(Yraw);
if NG ~= N*G; error('wrong specification of N and G'); end  % Check dimensions
Ylag = mlag2(Yraw,p);
n =(1+p*NG)*NG;         % total number of regression coefficients
k = p*NG+1;             % number of coefficients in each equation
m = NG*(NG-1)/2;        % # of parameters in A
% Correct time-series dimension of Y due to taking lags
Y = Yraw(p+1:Traw,:); 
T=Traw-p;
X = [ones(T,1) Ylag(p+1:Traw,:)];   % VAR data (RHS) matrix on the original specification    


% ====| Examine restrictions
% Here I use a different function to directly otain the indices for C-S and DI restrictions. We do
% not check SI restrictions since we have SV.

n_CS = N*(N-1)/2;          % Number of C-S restrictions
n_DI = N*(N-1)*p;   % Number of DI restrictions
%n_DI = N;   % Number of DI restrictions
[index_CS,index_restr_DI]=get_index_S4(N,G,p);

% The "index_CS" variable indexes the matrices A_{1}^{i} of country i.
% We need to test pairs of restrictions of the form A_{1}^{i} = A_{1}^{j}
% in order to test homogeneity of countries i and j.
% a) First create pairs of countries:
pairs_index = combntns(1:N,2);   % Index of pairs
% b) Second obtain index of CS restrictions. For each pair of countries, we
% are testing the equivalence of the GxG matrices A_{1}^{i}.
index_restr_CS = cell(n_CS,1);
for ir = 1:n_CS
    temp =  index_CS(:,:,pairs_index(ir,:));
    temp1 = temp(:,:,1); temp2 = temp(:,:,2);
    index_restr_CS{ir,1} =  [temp1(:),temp2(:)];
end

% ====| Set priors
% initialize alpha with OLS estimates 
%alpha =0.01*ones(n,1);

% alpha ~ N(0,D^-1), where D^-1 = diag(tau2)
h_i = ones(n,1);
tau2_0=0.01;
tau2_1=4;
ksi2_0=0.01;
ksi2_1=4;
%tau2 = zeros(n_DI,1);
%ksi2 = zeros(n_CS,1);
%c1 = 1e-6;
%c2 = 1e-6;

% tau^-2 ~ Gamma(rho1,rho2);
%rho1 = 2;
%rho2 = 0.1;

% ksi^-2 ~ Gamma(miu1,miu2);
%miu1 = 2;
%miu2 = 0.2;

% gamma_j ~ Bernoulli(1,p_j)
%p_j_DI    = 0.5*ones(n_DI,1);
%p_j_CS    = 0.5*ones(n_CS,1);

% Initialize parameters
gamma_DI = ones(n_DI,1);
gamma_CS = ones(n_CS,1);
GAMMA2 = cell(n_CS,1);
for i = 1:n_CS; GAMMA2{i,1} = speye(n); end
GAMMA = speye(n);
for i = 1:n_CS
    GAMMA = GAMMA*GAMMA2{i,1};
end

% Create storage matrices for posteriors
Y_pred = zeros(nsave,NG,h);

% Initialize stochastic volatility and A matrix
iVa=0.1*eye(m);     % prior precision for A
inv_A=eye(NG);  
ARresid=zeros(Traw-1,NG);
for i=1:NG
    [~,~,ARresid(:,i)]=OLS_PVAR(Yraw(:,i),1,1,1);
end
htemp=mean(ARresid'.^2,2);
sqrt_ht=sqrt(repmat(htemp',T,1));       
Vol_states=2*log(sqrt_ht);Volt= sqrt_ht.^2;                        
PHI_=0.0001*eye(NG);  
d_PHI = NG+2;                
s_PHI = 0.01*eye(NG); 
Vol_0mean = zeros(NG,1);   
Vol_0var  = 100*eye(NG);   

% ===============| GIBBS SAMPLER

for irep = 1:ntot
  
   
    %------------------------------------------------------
    % STEP 1: Update VAR coefficients alpha from Normal
    %------------------------------------------------------
  %  stationary=0;
  %  while stationary==0
        for kk = 1:n_DI
            ind_temp = index_restr_DI(:,:,kk)';
        % sample tau^-2 from Gamma
      %  r1 = rho1 + .5*G*(k-G*p-1);
      %  r2 = rho2 + sum(alpha(ind_temp(:)).^2)/(2*(c1^(1-gamma_DI(kk,1))));
      %  tau2(kk,1) = min(1./gamrnd(r1,1./r2),1000);
            if gamma_DI(kk,1) == 0
                h_i(ind_temp(:)) = tau2_0;
            elseif gamma_DI(kk,1) == 1
                h_i(ind_temp(:)) = tau2_1;
            end      
        end
    
        for kk = 1:n_CS
            ind_temp = index_restr_CS{kk,1};
        % Sample ksi^-2 from Gamma
     %   r1 = miu1 + .5*G*(G*p+1);
     %   r2 = miu2 + sum((alpha(ind_temp(:,1))-alpha(ind_temp(:,2))).^2)/(2*(c2^(1-gamma_CS(kk,1))));
     %   ksi2(kk,1) = min(1./gamrnd(r1,1./r2),1000);        
            if gamma_CS(kk,1) == 0
                h_i(ind_temp(:,1)) =ksi2_0;
            elseif gamma_CS(kk,1) == 1
                h_i(ind_temp(:,1)) = ksi2_1;
            end       
        end
    
        D = diag(1./h_i);
        [mutemp,Sigtemp] = getpsmean(Y,X,Volt,inv_A);
        Delta_alpha=(Sigtemp + D)\speye(n);
        miu_alpha = Delta_alpha*(mutemp);    
        alpha = GAMMA*miu_alpha + chol(Delta_alpha)'*randn(n,1);
        alpha_mat = reshape(alpha,k,NG);
    %    if max(abs(eig(alpha_mat(2:k,:)' ))) < 1.1
    %        stationary = 1;
    %    end
  %  end
    RESID=Y-X*alpha_mat;
    
    %----------------------------------------------------------------------
    % STEP 2: Update DI and CS restriction indexes of alpha from Bernoulli
    %---------------------------------------------------------------------- 
    % 1) Examine dynamic interdependencies (DI) 
    p_j_DI = repmat(betarnd(1 + sum(gamma_DI==1),1 + sum(gamma_DI~=1)),n_DI,1);
    for kk = 1:n_DI       
        ind_temp = index_restr_DI(:,:,kk)';
        u_i1 = mvnpdf(alpha(ind_temp(:)),zeros(G*G,1),tau2_0*eye(G*G))*p_j_DI(kk);
        u_i2 = mvnpdf(alpha(ind_temp(:)),zeros(G*G,1),tau2_1*eye(G*G))*(1- p_j_DI(kk));
        gst = u_i2./(u_i1 + u_i2);
        gamma_DI(kk,1) = bernoullirnd(gst);   
    end    
    
    % 2) Examine cross-sectional (CS) heterogeneities
    % update gamma, one at a time (their prior is independent)
    p_j_CS = repmat(betarnd(1 + sum(gamma_CS==1),1 + sum(gamma_CS~=1)),n_CS,1);
    for kk = 1:n_CS       
        ind_temp = index_restr_CS{kk,1};
        u_i1 = mvnpdf(alpha(ind_temp(:,1)),alpha(ind_temp(:,2)),ksi2_0*eye(G*(G*p+1)))*p_j_CS(kk);
        u_i2 = mvnpdf(alpha(ind_temp(:,1)),alpha(ind_temp(:,2)),ksi2_1*eye(G*(G*p+1)))*(1- p_j_CS(kk));
        gst = u_i2./(u_i1 + u_i2);
        gamma_CS(kk,1) = bernoullirnd(gst);
        if gamma_CS(kk) == 0
            for d_G = 1:G*(G*p+1)  
                GAMMA2{kk,1}(ind_temp(d_G,1),ind_temp(d_G,1)) = 0;
                GAMMA2{kk,1}(ind_temp(d_G,2),ind_temp(d_G,2)) = 1;                   
                GAMMA2{kk,1}(ind_temp(d_G,1),ind_temp(d_G,2)) = 1;
            end
        else
            for d_G = 1:G*(G*p+1)
                GAMMA2{kk,1}(ind_temp(d_G,1),ind_temp(d_G,1)) = 1;
                GAMMA2{kk,1}(ind_temp(d_G,2),ind_temp(d_G,2)) = 1;
                GAMMA2{kk,1}(ind_temp(d_G,1),ind_temp(d_G,2)) = 0;
            end
        end
        GAMMA = speye(n);
        for i = 1:n_CS
            GAMMA = GAMMA*GAMMA2{i,1};
        end
    end
        
%     %------------------------------------------------------------------
%     % STEP 3: Update coefficients in A
%     %------------------------------------------------------------------
     count_E = 0;
     E = zeros(T*NG,m);
     for jj=1:NG-1
         E(jj+1:NG:end,count_E+1:count_E+jj) = -RESID(:,1:jj);        
         count_E = count_E+jj;
     end
     iD = sparse(1:T*NG,1:T*NG,reshape(1./(sqrt_ht.^2)',1,T*NG));
     Ka = iVa + E'*iD*E;
     a_hat = Ka\(E'*iD*reshape(RESID',T*NG,1));   
     a = a_hat + chol(Ka,'lower')'\randn(m,1);  
     A_=chofac(NG,a);
     inv_A=A_\speye(NG);
     
     % --------------------------------------------------------
     %  Step 4: Draw volatility states using KSC algorithm
     % --------------------------------------------------------
     Vol_states  = KSC(log((RESID*A_').^2 + 1e-6),Vol_states,PHI_,Vol_0mean,Vol_0var);
     Volt=exp(Vol_states);
     sqrt_ht  = exp(Vol_states/2); 

     % ----------------------------------------------------------
     %  Step 5: Draw error variance for volatility innovations
     % ----------------------------------------------------------
    eta  = Vol_states(2:end,:) - Vol_states(1:end-1,:); 
    temp = chol(inv(s_PHI + eta'*eta),'lower')*randn(NG,T + d_PHI);
    PHI_ = (temp*temp')\speye(NG);
  
    
     % simulate y path: use draws from posterior distribution
    if irep>burnin
        Y_pred(irep-burnin,:,:)=get_paths_VARSV(alpha,h,p,NG,PHI_,inv_A,Vol_states,Y);
    end
    
end
