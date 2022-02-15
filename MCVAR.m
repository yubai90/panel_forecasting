% ---------------------------------------------------------------------------
% This MAIN m file implements estimation and forecasting with multi-country
% VAR-SV models and hierarchical shrinkage.
% There are 5 options available:
% (1) Minnesota prior;
% (2) Horseshoe prior;
% (3) Normal-Gamma prior;
% (4) Normal-Gamma-Gamma prior
% (5) Stochastic search  specification selection (SSSS) prior
% Reference: Bai, Y., Carriero, A., Clark, T. E., & Marcellino, M. G. (2022). Macroeconomic forecasting 
% in a multi-country context. Journal of Applied Econometrics, forthcoming.
% ---------------------------------------------------------------------------
% Notes:
% (1) There are several m files in the subfolder "MACVAR" which are used
% to "replicate" part of the results in the robustness section. These include alternative prior
% specifications with horseshoe prior and Minnesota-type prior with
% Horseshoe shrinkage. You need to mannualy change the main files to obtain
% those results.
% (2) These are not replication files, but the way to show how our MCMC
% estimation and forecasting works, since due to Licensing restrictions, we
% are not able to provide the GFD data we used. However, to allow researchers to
% run the programs without the omitted series,we have added a random number
% from (0,1) to the original series of interest rate data, using the excel
% command RAND.
% --------------------------------------------------------------------------
% Contact information:  yu.bai [AT] unibocconi.it
% --------------------------------------------------------------------------

clear
clc
rng(1,'twister')

% Add path of random number generators
addpath(genpath('functions'))
addpath('data');

% Inputs
N=7;                          % cross-sectional units
G=3;                          % variables in each country
NG=N*G;                       % total number of variables
h=12;                         % forecast horizon
reps=3000;                   % total reps
burnin=1000;                 % burn in
nsave=reps-burnin;            % draws saved
L=4;                          % number of lags
theta=[0.003 0.0015 100 0.001];     % Shrinkage hyperparameter for MIN prior: This will only be used as an initialization.
model=5;                      % 1 Minnesota prior;
                              % 2 Horseshoe prior;
                              % 3 Normal-Gamma prior;
                              % 4 Normal-Gamma-Gamma prior
                              % 5 Stochastic search  specification selection (SSSS) prior

% Load data 
load dataQ.mat
Yraw=data_panel;

t = size(Yraw,1);
M = size(Yraw,2);
t0 = 83;                  % forecasting evaluation period start from 1995Q1

% Reshape the dataset in order to construct the parallel loop
Yreg=cell(t-h-t0+1,1);
Yt=cell(t-h-t0+1,1);
for i=1:t-h-t0+1
    Yreg{i}=Yraw(i:i+t0-1,:);
    Yt{i}=Yraw(i+t0:i+t0+h-1,:);
end

% Store the forecasting results: both for each country and each model
MSFE_SV = zeros(NG,h,t-h-t0+1);
CRPS_SV=zeros(NG,h,t-h-t0+1);

% Variable index
index1=1:3:21;
index2=2:3:21;
index3=3:3:21;
index=[index1 index2 index3];

tic
% Estimation and forecasting
for nMC=1:t-h-t0+1
    
    disp(nMC+t0-1)
        
    Y = Yreg{nMC,1};
    Y_f = Yt{nMC,1};
    T=nMC+t0-1;

         
    % Estimation and forecasting
    switch model
        case 1
            Y_predcon_sim_SV= MCVAR_MIN_SV(Y,nsave,burnin,N,G,L,h,theta);
        case 2
            Y_predcon_sim_SV= MCVAR_HS_SV(Y,reps,burnin,N,G,L,h);  
        case 3
            Y_predcon_sim_SV= MCVAR_NG_SV(Y,reps,burnin,N,G,L,h);  
        case 4
            Y_predcon_sim_SV= MCVAR_NGG_SV(Y,reps,burnin,N,G,L,h);  
        case 5
            Y_predcon_sim_SV= MCVAR_SSSS_SV(Y,N,G,1,nsave,burnin,reps,h);           
    end
                  
    Y_pred_con_SV=squeeze(mean(Y_predcon_sim_SV,1));
    
    % point forecast evaluation based on MSFE
    MSFE_SV(:,:,nMC)=(Y_pred_con_SV - Y_f').^2;
         
    % density forecast
    [~,CRPS_SV(:,:,nMC)]=densfore(Y_predcon_sim_SV,Y_f,h,NG); 
     
end

% Output growth
MSFE_y= MSFE_SV(index1,:,:);
CRPS_y=CRPS_SV(index1,:,:);

% Inflation
MSFE_cpi= MSFE_SV(index2,:,:);
CRPS_cpi=CRPS_SV(index2,:,:);

% Interest rates
MSFE_r= MSFE_SV(index3,:,:);
CRPS_r=CRPS_SV(index3,:,:);


toc