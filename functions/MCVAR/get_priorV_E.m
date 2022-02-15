function sigmab = get_priorV_E(omega,gamma,phi)
% Get prior variance matrix with dimension K*NG. omega and gamma are
% collections of hyperparameters. Phi are global shrinkage parameters which takes the form:
% 2c/(ab), and we have NG *1 global shrinkage parameters 

local_V=omega./gamma;
[K,NG]=size(local_V);
sigmab=zeros(K,NG);

for eq=1:NG
    sigmab(:,eq)=phi(eq)*local_V(:,eq);
end


