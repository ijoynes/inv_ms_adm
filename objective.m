function [f,g]=objective(x)

global_vars

assert(nargin == 1);
assert(nargout <= 2);


%nt = length(t);

s = volumetric_emission_rate(x);
assert( length(s) == length(s_b), 's must have the same length as s_b');

c = solve_advection_diffusion_equation(s, t, H, save_flag, conc_dir, oper_dir,'m',time_offset,nt_max,0.5);

assert( all(size(c) == size(c_star)), 'c must have the same size as c_star');
f_obs = 0.5*sum(sum((c-c_star).^2));  % f_obs = 0.5*(Hc(x)-c*)'R^-1(Hc(x)-c*), R^-1 is assumed to be the identity matrix

f = f_obs + 0.5*reg_par*((s-s_b)'*(B_inv*(s-s_b)));

% solve for the sensitivity of the objective function to the receptor
% observations 
if nargout == 1
  g_obs = nan(length(x),1);
  g = nan(length(x),1);
else
  g_obs = solve_adjoint_advection_diffusion_equation(c, c_star, t, H, save_flag, adj_dir, oper_dir,time_offset,nt_max, 0.5);
  %df_dx = CostFncGrad2(t, nNodes);
  % add sensitivity of objective function to the regularization term
  g = g_obs + reg_par*(B_inv*(s-s_b));
end

nfevals = nfevals + 1;

% This is a hack to ensure that the results for the initial starting 
% conditions of the minimization routine are saved.  The first iteration
% of L-BFGS-B typically makes 2 or 3 objective function and gradient 
% evaluations while performing a line-search to ensure that the 
% objective is reduced in the next iteration.
if nPass == 0
  
  g_proj = projgr(x, g, lb, ub, nbd);
  r  = compute_correlation_coefficient(c, c_star);
  r2 = compute_coefficient_of_determination(c, c_star);
  
  m = dot(space_int_wgt, s);
  m_norm = m/m_star;

  x_0 = x;
  s_0 = s;
  f_0 = f;
  f_obs_0 = f_obs;
  g_0 = g;
  g_proj_0 = g_proj;
  g_obs_0 = g_obs;

  r_0 = r;
  r2_0 = r2;
  m_0 = m;
  c_0 = c;
  
  fprintf('   f0 = %e [ kg^2 * m^-6 ]\n', f_0);
  fprintf('|pg0| = %e [ kg * s * m^-3 ]\n' , norm(g_proj_0,Inf));
  fprintf('   m* = %e [ kg * s^-1 ]\n', m_star);
%  fprintf('\n');
%  fprintf('-------------------------------------------------------------------------------------------------------\n');
%  fprintf('|       Date/Time      | Iter |    fk/f0    | |pgk|/|pg0|  |      r^2     |      r      |    m/m*     |\n');
%  fprintf('-------------------------------------------------------------------------------------------------------\n');
%                                           nfevals |
  fprintf('\n');
  fprintf('----------------------------------------------------------------------------------------------------------------\n');
  fprintf('|       Date/Time      | Iter | nfevals |    fk/f0    | |pgk|/|pg0| |      r^2     |      r      |    m/m*     |\n');
  fprintf('----------------------------------------------------------------------------------------------------------------\n');

  if r2 < 0
    fprintf('| %s | %4d | %7d | %8.5e | %8.5e | %9.5e | %8.5e | %8.5e |\n', datestr(now), nPass, nfevals, f/f_0, norm(g_proj,Inf)/norm(g_proj_0,Inf),r2,r, m_norm);
  else
    fprintf('| %s | %4d | %7d | %8.5e | %8.5e |  %8.5e | %8.5e | %8.5e |\n', datestr(now), nPass, nfevals, f/f_0, norm(g_proj,Inf)/norm(g_proj_0,Inf),r2,r, m_norm);
  end
  file_num = generate_file_num(nPass, max_iters);
  results_path = fullfile(iter_dir, [iter_label, file_num, '.mat']);
  save(results_path, 'f', 'f_obs', 'g', 'g_proj', 'g_obs', 's', 'x', 'c', 'r', 'r2', 'm', 'm_norm','nfevals');
  nfevals = 0;
end

nPass = nPass + 1;

