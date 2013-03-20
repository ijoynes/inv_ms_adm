function r2 = compute_coefficient_of_determination(f, y)
% f   model response
% y   observations

assert( isnumeric(f) );
assert( isnumeric(y) );
assert( all(size(f) == size(y)) );

f  = reshape(f, 1, [])';
y  = reshape(y, 1, [])';

y_mean = mean(y);
den = sum((y-y_mean).^2);

if den == 0
  r2 = 0;
else
  r2 = 1 - sum((y-f).^2)/sum((y-y_mean).^2);
end