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
%         the number of spatial dimensions of the mesh.  Q must fall
%         within the range 1<=Q<=3.
%
%   x   R-by-Q matrix representing the coordinates of the receptors.
%         Where R is the number of receptors in the domain and Q 
%         represents the number of spatial dimensions of the mesh.  Q 
%         must fall within the range 1<=Q<=3.
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
assert(size(X, 2) == size(x, 2));

% Consider testing T for integers

assert(isnumeric(T));
assert(isnumeric(X));
assert(isnumeric(x));
assert(size(X,2) == size(x,2));
assert(size(T,2) == size(X,2)+1); % ensure T represents simplexes

% Predetermine constants        % Matrix sizes from description
nReceptors = size(x, 1);        % R
nNodes     = size(X, 1);        % P
nElements  = size(T, 1);        % M
nNodesPerElement = size(T, 2);  % N

% Find the indices of the containing triangles and the corresponding
% barycentric coordinates (shape function values) with a brute force 
% search.
B = nan(nReceptors, nNodesPerElement);
T_ids = nan(nReceptors,1);
for i = 1 : nElements
    A = [ X(T(i, :), :)'; ones(1, nNodesPerElement) ];
    b = [ x'; ones(1,nReceptors) ];
    bcc = A\b;
    for j = 1 : nReceptors
      if all( -eps <= bcc(:, j) ) && all(bcc(:, j) < 1 + eps )
          T_ids(j) = i;
          B(j, :) = bcc(:, j)';
      end
    end
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
H = sparse(row, col, Hv, nReceptors, nNodes);