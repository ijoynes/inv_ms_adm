function [B_inv] = compute_inverse_source_error_covariance_matrix(tri,xy)
%compute_inverse_source_error_covariance_matrix constructs the inverse 
% of the covariance matrix of the estimated background error aka B^-1.
% This matrix is used by the objective function and its gradient to 
% penalize source distributions that deviate from an assume prior 
% spatial distribution (s_b).
%
% This implementation of B^-1 numerically integrates the square of the 
% residual between the discrete representation of the candidate source 
% distribution (s_k) and the discrete representation of the assumed 
% prior source distribution (s_b), such that this relationship is true.
%
% int((s_k-s_b)^2, Omega) = (s_k-s_b)'*B^-1*(s_k-s_b)
%
% INPUTS:
%
% OUTPUTS:
%   B_inv   A sparse symmetric matrix positive definite matrix that 
%           can be used to numerically integrate the discretized spatial
%           distributions of the square of the residual error 

assert( nargin == 2);
assert(nargout == 1);

assert(isnumeric(tri));
assert(ismatrix(tri));
assert(size(tri,2)==3);
assert(all(~isnan(tri)));

assert(isnumeric(xy));
assert(ismatrix(xy));
assert(size(xy,2)==2);
assert(all(~isnan(xy)));


nTris = size(tri,1);
Bv = nan(nTris*3*3,1);
row = nan(nTris*3*3,1);
col = nan(nTris*3*3,1);
index = 0;

for i = 1 : nTris
    Be = [ 2 1 1; 1 2 1; 1 1 2] * det([ones(3,1), xy(tri(i,:),:)]);
    for j = 1 : 3
        for k = 1 : 3
            index = index + 1;
            Bv(index) = Be(j,k);
            row(index) = tri(i,j);
            col(index) = tri(i,k);
        end
    end
end

Bv = Bv/24;
B_inv = sparse(row,col,Bv,nNodes,nNodes);
