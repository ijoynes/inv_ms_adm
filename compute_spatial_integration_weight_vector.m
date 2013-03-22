function sp_int_wgt = compute_spatial_integration_weight_vector(T, X)
% The integration of a scalar over the entire domain can be 
% approximated numerically by summing the integral of the scalar over
% each element.  The integral of the scalar over an element can be 
% approximated summing the integrals of the product of finite element 
% shape functions and the nodal value of the scalar field.  This 
% integration operation can be streamlined by pre-computing a spatial 
% integration weight vector (sp_int_wgt in this case), which is a
% vector of weights based on the area (or volume) of elements that a 
% given mesh node shares, and indicates the influence that a given 
% node has on the spatial integration of the scalar field.  With this 
% precomputed spatial integration weight vector the spatial 
% integration of any scalar field (represented as a vector of nodal 
% values) is computed with the dot product between these to vectors.
% For instance the total emission rate of a source volumetric emission 
% rates would be computed as follows:
%
% total_emission_rate = dot(sp_int_wgt, s);
%
% Construct the spatial integration weight vector

assert( nargin == 2);
assert(nargout == 1);

assert(isnumeric(T),'Node tessellation list must contain only numeric values.');
assert(ismatrix(T), 'Node tessellation list must be a matrix.');
assert(all(isfinite(T(:))), 'Node tessellation list must contain only finite values.');      % there are no NaNs or Infs
assert(all(rem(T(:),1)==0), 'Node tessellation list must contain integers.');      % the values of T are integers

assert(isnumeric(X), 'Mesh node coordinate list must contain only numeric values.');
assert(ismatrix(X), 'Mesh node coordinate list must be a matrix.');
assert(all(isfinite(X(:))), 'Mesh node coordinate list must contain only finite values.');
assert(size(T,2) == size(X,2)+1, 'Tessellation must represent simplexes.');  % the mesh represents simplexes

[nSimplex] = size(T,1);
[nNodes, nDims] = size(X);

assert(all(1 <= T(:)) & all(T(:) <= nNodes), ...
  'Node tessellation list contains invalid values.');

sp_int_wgt = zeros(nNodes,1);
for i = 1 : nSimplex
  sp_int_wgt(T(i,:)) = sp_int_wgt(T(i,:)) + ...
                             det([ones(nDims+1,1) X(T(i,:),:) ] );
end

assert(all(sp_int_wgt > 0), 'Spatial integration weight vector contains non-positive values.');
sp_int_wgt = sp_int_wgt/factorial(nDims+1);