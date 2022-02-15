function A_chol=chol_new(A)
% Perform cholesky decomposition of the covariance matrix. If matrix is not
% positive definite, find the nearest SPD matrix in Frobenius norm.

[R,p]=chol(A);
if p==0
    A_chol=R;
else
    A_new=nearestSPD(A);
    A_chol=chol(A_new);
end