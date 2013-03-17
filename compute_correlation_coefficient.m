function r = compute_correlation_coefficient(x, y)
% x   model response
% y   observations

assert( isnumeric(x) );
assert( isnumeric(y) );
assert( all(size(x) == size(y)) );

x  = reshape(x, 1, [])';
y  = reshape(y, 1, [])';

Sx  = sum(x);
Sy  = sum(y);
Sxx = sum(x.^2);
Syy = sum(y.^2);
Sxy = sum(x.*y);
n   = numel(y);
den = sqrt(n*Sxx-Sx^2)*sqrt(n*Syy-Sy^2);

if den == 0 % if the denominator is zero then the correlation 
  r = 0;    %   coefficient is assumed to be 0
else
  r = (n*Sxy-Sx*Sy)/den;
end