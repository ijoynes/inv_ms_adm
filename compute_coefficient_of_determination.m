function r2 = compute_coefficient_of_determination(c_k, c_star)
  ybar = mean(mean(signal_o));
  r2 = 1 - sum(sum((signal_o-signal_c).^2))/sum(sum((signal_o-ybar).^2));