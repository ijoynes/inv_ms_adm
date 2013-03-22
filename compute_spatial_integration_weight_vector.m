function sp_int_wgt = compute_spatial_integration_weight_vector(tri,xy)
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

assert(isnumeric(tri));
assert(ismatrix(tri));
assert(size(tri,2)==3);
assert(all(~isnan(tri)));

assert(isnumeric(xy));
assert(ismatrix(xy));
assert(size(xy,2)==2);
assert(all(~isnan(xy)));

sp_int_wgt = zeros(nNodes,1);
for i = 1 : nTris
  sp_int_wgt(tri(i,:)) = sp_int_wgt(tri(i,:)) + ...
                             det([ones(3,1) xy(tri(i,:),:) ] );
end
sp_int_wgt = sp_int_wgt/6;