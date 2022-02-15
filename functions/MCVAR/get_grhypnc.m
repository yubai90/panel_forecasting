function [olag,cvlag,di] = get_grhypnc(parameters,index_dom,index_di,N,G,L)
% We need to transfer K*NG (hyper)parameters to different blocks. Then, we
% use this in other steps for MCMC.
% We have:
% (1) own lags; (2) cross-variable lages; (3) DI (k2 * N matrix); (no elements for A here)

NG=N*G;                      % number of variables
k2=G*L*(N-1);                % DI in each equation
olag=[];
cvlag=[];
ditemp=[];

for ii=1:N
    indx_v=(ii-1)*G+1:ii*G;
    indexc=index_dom(ii,:);
    for jj=1:G
        inda=1:G*L;
        ind=jj:G:G*L;
        inda(ind)=[];
        olagtemp=parameters(indexc(ind),indx_v(jj));
        olag=[olag;olagtemp];
        cvlatemp=parameters(indexc(inda),indx_v(jj));
        cvlag=[cvlag;cvlatemp];
    end
end

for ii=1:NG
    inc=ceil(ii/G);           % get country index
    ind_di=index_di(inc,:);
    dip=parameters(ind_di,ii);
    ditemp=[ditemp;dip];
end
di=reshape(ditemp,k2*G,N);
    




