function sigmab = get_priorV2(omega,index_dom,index_di,N,G,L)
% Get prior variance matrix with dimension K*NG. omega are
% collections of hyperparameters. Phi are global shrinkage parameters which takes the form:
% 2c/(ab)(we have N+3, 1 for intercept, 1 for own lags, 1 for cross-variable lags, N for DI).
% This one is used for NG and HS priors.

local_V=omega;
[K,NG]=size(local_V);
sigmab=zeros(K,NG);
sigmab(1,:)=local_V(1,:);

for ii=1:N
    indx_v=(ii-1)*G+1:ii*G;
    indexc=index_dom(ii,:);
    for jj=1:G
        inda=1:G*L;
        ind=jj:G:G*L;
        inda(ind)=[];
        sigmab(indexc(ind), indx_v(jj))= local_V(indexc(ind), indx_v(jj));
        sigmab(indexc(inda), indx_v(jj))=local_V(indexc(inda), indx_v(jj));
    end
end


for ii=1:NG
    inc=ceil(ii/G);           % get country index
   % ind_csh=index_dom(inc,:);
   % csh_temp=local_V(ind_csh,ii);
   % sigmab(ind_csh,ii)=phi(1)*csh_temp;
    ind_di=index_di(inc,:);
    dip=local_V(ind_di,ii);
    sigmab(ind_di,ii)=dip;
end
