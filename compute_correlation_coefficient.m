function r = compute_correlation_coefficient(c_k, c_star)

Sx = sum(sum(signal_o));
Sy = sum(sum(signal_c));
Sxx = sum(sum(signal_o.^2));
Syy = sum(sum(signal_c.^2));
Sxy = sum(sum(signal_o.*signal_c));
n = numel(signal_o);

if Sy == 0
  r = 0;
else
  r = (n*Sxy-Sx*Sy)/(sqrt(n*Sxx-Sx^2)*sqrt(n*Syy-Sy^2));
end