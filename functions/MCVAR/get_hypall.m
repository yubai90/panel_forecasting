function parameters = get_hypall(int,olag,cvlag,di,index_dom,index_for,N,G,L)
% From the grouped parameters to K*NG matrix

NG=N*G;
K=1+NG*L;
%k_m=G*L+1;
olag_mat=reshape(olag,L,NG);
cvlag_mat=reshape(cvlag,(G-1)*L,NG);
k_d=G*L*(N-1);
parameters=zeros(K,NG);
parameters(1,:)=int;

for ii=1:N
    indx_v=(ii-1)*G+1:ii*G;
    olag_temp=olag_mat(:,indx_v);
    cvlag_temp=cvlag_mat(:,indx_v);
    indexc=index_dom(ii,:);
    for jj=1:G
        inda=1:G*L;
        ind=jj:G:G*L;
        inda(ind)=[];
        parameters(indexc(ind),indx_v(jj))=olag_temp(:,jj);
        parameters(indexc(inda),indx_v(jj))=cvlag_temp(:,jj);
    end
end
        
        
for ii=1:N
    indx_v=(ii-1)*G+1:ii*G;
  %  csh_temp=reshape(csh((ii-1)*k_m*G+1:ii*k_m*G),k_m,G);
    di_temp=reshape(di(:,ii),k_d,G);
  %  parameters(index_dom(ii,:),indx_v)=csh_temp;
    parameters(index_for(ii,:),indx_v)=di_temp;
end
    


