function  yhatx = get_paths_VARSV(PAI,h,L,N,PHI,inv_A,hlast,y)
% simulate the y path: VAR & SV

hhat=zeros(h+L,N);
hhat(L,:)=hlast(end,:);
  
yhat=zeros(h+L,N);
yhat(1:L,:)=y(end-L+1:end,:);
Cchol=chol(PHI);

for m=L+1:h+L
    % simulate the volatility path
    hhat(m,:)=hhat(m-1,:)+randn(1,N)*Cchol;
     
    % simulate y path
    xhat=[];
    for i=1:L
        xhat=[xhat yhat(m-i,:)];
    end
    xhat=[1 xhat];
    smat=inv_A*diag(exp(hhat(m,:)))*inv_A';
    yhat(m,:)=xhat*reshape(PAI,N*L+1,N)+randn(1,N)*chol_new(smat); % find nearest cholesky decomposition if smat is nearly singular 
end

yhatx=yhat(L+1:end,:);
yhatx=yhatx';

