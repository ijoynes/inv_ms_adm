function B_inv = compute_inverse_source_error_covariance_matrix(T, X)
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
%   T     M-by-N matrix representing the connectivity of the mesh.
%           Where M is the number of elements and N is the number of
%           nodes per element.
%             N = 2 line elements
%             N = 3 triangular elements
%             N = 4 tetrahedral elements
%
%           T(i,j) contains the global node index of the jth node of the
%           ith element.
%
%   X     P-by-Q matrix representing the coordinates of the mesh nodes.
%           Where P is the number of nodes in the mesh and Q represents 
%           the number of spatial dimensions of the mesh.  Q must fall
%           within the range 1<=Q<=3.
%
% OUTPUTS:
%   B_inv   A sparse symmetric matrix positive definite matrix that can
%             be used to numerically integrate the discretized spatial
%             distributions of the square of the residual error 

assert( nargin == 2, '2 arguments must be supplied');
assert(nargout == 1);

assert(isnumeric(T));
assert(ismatrix(T));
assert(all(isfinite(T(:))));

assert(isnumeric(X));
assert(ismatrix(X));
assert(all(isfinite(X(:))));
assert(size(T,2) == size(X,2)+1);

[nTris N] = size(T);
nNodes    = size(X,1);

assert(all(1<= T(:)) & all(T(:)<=nNodes));

Bv  = nan(nTris*N*N,1);
row = nan(nTris*N*N,1);
col = nan(nTris*N*N,1);
idx = 0;

for i = 1 : nTris
  Be = (ones(N) + eye(N)) * det([ones(N,1), X(T(i,:),:)]);
  for j = 1 : N
    for k = 1 : N
      idx = idx + 1;
      Bv(idx) = Be(j,k);
      row(idx) = T(i,j);
      col(idx) = T(i,k);
    end
  end
end
Bv = Bv/factorial(N+1);

assert(all(isfinite(Bv)));
assert(all(isfinite(row)));
assert(all(isfinite(col)));

B_inv = sparse(row,col,Bv,nNodes,nNodes);