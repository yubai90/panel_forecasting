function [index_CS,index_restr_DI]=get_index_S4(N,G,L)
% obtain the index of parameters for both C-S and DI restrictions: with
% intercept

NG=N*G;
n=N*G*L+1;
%index_restr_DI = zeros(G,n-G*L-1,N);
index_restr_DI = zeros(G,G,N*(N-1)*L);
index_CS=zeros(G,G*L+1,N);   
index_tempdi=[];
index_c=1:n*G:n*NG;

for i_country = 1:N
    
    index_temp2=[];
    for i_variable=1:G
        index_temp=[];
        index_all=(i_country-1)*n*G+(i_variable-1)*n+1:(i_country-1)*n*G+(i_variable-1)*n+n;
        for p=1:L
            index_tem = G*n*(i_country-1)+2+(i_country-1)*G+(p-1)*NG:G*n*(i_country-1)+G+1 ...
                  +(i_country-1)*G+(p-1)*NG;
            index_temp=[index_temp index_tem];
        end
        index_temp=[index_c(i_country) index_temp];
        index_all(index_temp-(i_country-1)*n*G)=[];
        index_CS(i_variable,:,i_country) = index_temp + (i_variable-1)*n;
        index_temp2=[index_temp2;index_all];
    end
    index_tempdi=[index_tempdi index_temp2];
  %  index_restr_DI(:,:,i_country)=index_temp2;
    
    for i=(i_country-1)*(N-1)*L+1:(N-1)*L*i_country
        index_restr_DI(:,:,i)=index_tempdi(:,(i-1)*G+1:i*G);
    end
    
end