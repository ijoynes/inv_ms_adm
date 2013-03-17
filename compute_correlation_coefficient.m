function r = compute_correlation_coefficient(c_k, c_star)
% consider changing the variable names to a more general format, 
% (i.e. f for the model and y for the observations).

assert( isnumeric(c_k) );
assert( isnumeric(c_star) );
assert( size(c_k,1) == size(c_star,1));
assert( size(c_k,2) == size(c_star,2));

Sy = sum(sum(c_k));

if Sy == 0
  r = 0;
else
  Sx  = sum(sum(c_star));
  Sxx = sum(sum(c_star.^2));
  Syy = sum(sum(c_k.^2));
  Sxy = sum(sum(c_star.*c_k));
  n = numel(c_star);
  r = (n*Sxy-Sx*Sy)/(sqrt(n*Sxx-Sx^2)*sqrt(n*Syy-Sy^2));
end