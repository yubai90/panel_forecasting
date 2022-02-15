% This function constructs the second moments of the Minnesota prior
% for multi-country VARs
% 

function [C,idx_kappa1,idx_kappa2,idx_kappa3] = get_C_MC(N,G,p,sig2)
n=N*G;
k_beta = n^2*p+n;
C = zeros(k_beta,1);
idx_kappa1 = [];
idx_kappa2 = [];
idx_kappa3 = [];
count = 1;

for ii = 1:n
    inc=ceil(ii/G);           % country index
    idxg=(inc-1)*G+1:inc*G;
    Ci = zeros(n*p+1,1);
        % construct Ci
    for j=1:n*p+1    
        l = ceil((j-1)/n); % lag length
        idx = mod(j-1,n);  % variable index
        if idx==0
            idx = n;
        end         
        if j==1 % intercept            
            Ci(j) = sig2(ii);
        elseif idx == ii % own lag
            Ci(j) = 1/l^2;
            idx_kappa1 = [idx_kappa1; count];
        elseif idx>=idxg(1) && idx<=idxg(end) % lag of other variables in the country            
            Ci(j) = sig2(ii)/(l^2*sig2(idx));
            idx_kappa2 = [idx_kappa2; count];
        elseif idx>idxg(end) || idx<idxg(1)
            Ci(j) = sig2(ii)/(l^2*sig2(idx));
            idx_kappa3 = [idx_kappa3; count];
        end 
            count = count + 1;
    end
    C((ii-1)*(n*p+1)+1:ii*(n*p+1)) = Ci;
end
end