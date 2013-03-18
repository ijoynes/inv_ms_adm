function [f,df]=objective_function(x_k)

%global simDir operDir domainPath paramPath sensorPath sourcePath ...
%        passPath tMax noise reg_par signal_o signal_c ...
%       spaceIntWeight m_star nPass c_star B_inv


global sim_dir oper_dir iter_dir t reg_par c_star B_inv H noise

assert(nargout == 2);

nt = length(t);

s_k = volumetric_emission_rate(x_k);

c_k = solve_advection_diffusion_equation(s_k, t, H, false, conc_dir, oper_dir);
c_k = add_obs_noise(c_k, noise);

f_obs = 0.5*sum(sum((c_k-c_star).^2));  % f_obs = 0.5*(Hc(x)-c*)'R^-1(Hc(x)-c*), R^-1 is assumed to be the identity matrix


file_num = generate_file_num(0, nt);
initial_source_path = fullfile(iter_dir [iter_label file_num '.mat']);
load(initial_source_path,'x');
s_b = volumetric_emission_rate(x);  % initial candidate source distribution


p = reg_par;
f = f_obs + 0.5*p*((s_k-s_b)'*(B_inv*(s_k-s_b)));



% solve for the sensitivity of the objective function to the receptor
% observations 
df_obs = solve_adjoint_advection_diffusion(c_k, c_star, t, H, false, adj_dir, oper_dir);
%df_dx = CostFncGrad2(t, nNodes);

% add sensitivity of objective function to the regularization term
df = df_obs + p*(B_inv*(s_k-s_b));


% This is a hack to ensure that the results for the initial starting 
% conditions of the minimization routine are saved.  The first iteration
% of L-BFGS-B typically makes 2 or 3 objective function and gradient 
% evaluations while performing a line-search to ensure that the 
% objective is reduced in the next iteration.
if nPass == 0
  s = s_k;
  c = c_k;
  x = x_k;
  r  = compute_correlation_coefficient(c_k, c_star);
  r2 = compute_coefficient_of_determination(c_k, c_star);

  m = dot(spaceIntWeight, E_approx);
  m_norm = m/m_star;
  fprintf('| %s | %4d | %8.6e |  %8.6e | %9.6e | %8.6e | %8.6e |\n', datestr(now), nPass, f, max(abs(df)),r2,r, m_norm);

  file_num = generate_file_num(nPass, nt);
  results_path = fullfile(iter_dir, [iter_label, file_num, '.mat']);
  save(results_path,'f','f_obs', 'df', 'df_obs', 's', 'x' 'c','r','r2','m','m_norm');
end

nPass = nPass + 1;