function H = compute_observation_matrix(T, X, x)
%compute_observation_matrix generates the sparse observation matrix H 
% for computing the receptor observations from the discretized 
% concentration field.  The observation matrix performs a spatial
% interpolation on the discretized concentration field to compute
% receptor observations.
%
% Only simplex element types are supported by this function:
%   - line (2 nodes, linear shape functions)
%   - triangular (3 nodes, linear shape functions)
%   - tetrahedral (4 nodes, linear shape functions)
%
% INPUTS:
%   T   M-by-N matrix representing the connectivity of the mesh.
%         Where M is the number of elements and N is the number of
%         nodes per element.
%           N = 2 line elements
%           N = 3 triangular elements
%           N = 4 tetrahedral elements
%
%         T(i,j) contains the global node index of the jth node of the
%         ith element.
%
%   X   P-by-Q matrix representing the coordinates of the mesh nodes.
%         Where P is the number of nodes in the mesh and Q represents 
%         the number of spatial dimensions of the mesh.
%
%   x   R-by-Q matrix representing the coordinates of the receptors.
%         Where R is the number of receptors in the domain and Q 
%         represents the number of spatial dimensions of the mesh.
%
% OUTPUS:
%   H   R-by-P sparse observation matrix for computing the receptor 
%         observations from the discretized concentration field.  The 
%         observation matrix H is intended to be used in the following
%         manner.
%           r_obs = H*c
%         Where: 
%           c     is a column vector of length P and represents the 
%                   discretized concentration field at a given time 
%                   step.
%           r_obs is a column vector of length R and represents the 
%                   receptor observations at a given time step.
%
% NOTES:
%   This function could be extended to support quadrilateral, 
%     hexahedral and triangular prism elements.
%
%   This function currently uses a brute force search to find the
%     elements which contain the receptors.  Other containing simplex 
%     search methods such as 'pointLocation' should be considered if 
%     the number of elements is so large that unreasonable run times 
%     are expected.

% Ensure that the number of spatial dimensions of the mesh nodes 
% coordinates and receptor coordinates match.
if nargin < 3
  assert(nargin == 3, 'Not enough input arguments.');
elseif nargin > 3
  assert(nargin == 3, 'Too many input arguments.');
end
assert(nargout <= 1, 'Too many output arguments.');

assert(ismatrix(T), 'Node tessellation list must be a matrix.');
assert(ismatrix(X), 'Mesh node coordinate list must be a matrix.');
assert(ismatrix(x), 'Receptor coordinate list must be a matrix.');

assert(isnumeric(T), 'Node tessellation list must contain only numeric values.');
assert(isnumeric(X), 'Mesh node coordinate list must contain only numeric values.');
assert(isnumeric(x), 'Receptor coordinate list must contain only numeric values.');

assert(all(isfinite(T(:))), 'Node tessellation list must contain only finite values.');      % there are no NaNs or Infs
assert(all(isfinite(X(:))), 'Mesh node coordinate list must contain only finite values.');
assert(all(isfinite(x(:))), 'Receptor coordinate list must contain only finite values.');

assert(all(rem(T(:),1)==0), 'Node tessellation list must contain integers.');      % the values of T are integers

assert(size(X, 2) == size(x, 2), 'Number of spatial dimensions of the mesh must match that of the receptors.');
assert(size(T,2) == size(X,2)+1, 'Tessellation must represent simplexes.'); % ensure T represents simplexes

% Predetermine constants        % Matrix sizes from description
nReceptors = size(x, 1);        % R
nNodes     = size(X, 1);        % P
[nElements, nNodesPerElement] = size(T);   % [M, N]

assert(all(1 <= T(:)) & all(T(:) <= nNodes), 'Node tessellation list contains invalid values.');

% Find the indices of the containing triangles and the corresponding
% barycentric coordinates (shape function values) with a brute force 
% search.
B = nan(nReceptors, nNodesPerElement);
T_ids = nan(nReceptors,1);
for i = 1 : nElements
  A = [ X(T(i, :), :)'; ones(1, nNodesPerElement) ];
  b = [ x'; ones(1,nReceptors) ];
  bcc = A\b;

  % This loop scales with the number of receptors.  For a small number 
  % of receptors the run-time is trivial, but not for a large number of 
  % receptors.  This loop has been replaced with a vectorized version 
  % that uses boolean logic indexing.
  %
  %for j = 1 : nReceptors
  %  if all( -eps <= bcc(:, j) ) && all(bcc(:, j) < 1 + eps )
  %    T_ids(j) = i;
  %    B(j, :) = bcc(:, j)';
  %  end
  %end
  %
  % Vectorized version of the previous loop that was removed.
  index = all(-eps <= bcc) & all(bcc<= 1 + eps);
  T_ids(index) = i;
  B(index, :) = bcc(:, index)';

  if all(~isnan(T_ids)) % early exit if the containing simplexes of 
    break;              %   all the receptors have been found
  end
end

% Compute the sparse observation matrix H
Hv  = nan(nReceptors * nNodesPerElement, 1);
row = nan(nReceptors * nNodesPerElement, 1);
col = nan(nReceptors * nNodesPerElement, 1);
index = 0;
for i = 1 : nReceptors
  for j = 1 : nNodesPerElement
    index = index + 1;
    Hv(index) = B(i, j);
    row(index) = i;
    col(index) =  T(T_ids(i), j);
  end
end

assert(all(isfinite(Hv)), 'Observation matrix contains invalid values.');
assert(all(isfinite(row)), 'Observation matrix contains invalid row indexes.');
assert(all(isfinite(col)), 'Observation matrix contains invalid col indexes.');

H = sparse(row, col, Hv, nReceptors, nNodes);