function [f,df_dx]=objective_function(x_k)

%global simDir operDir domainPath paramPath sensorPath sourcePath ...
%        passPath tMax noise reg_par signal_o signal_c ...
%       spaceIntWeight m_star nPass c_star B_inv


global sim_dir oper_dir t reg_par c_star B_inv H

assert(nargout == 2);

s_k = volumetric_emission_rate(x_k);
%load(paramPath)


working_dir = fullfile( sim_dir, 'Concentration');
c_k = solve_advection_diffusion_equation(s_k, t, H, working_dir, oper_dir);
%SolveConcentrationTransport(t, s_k, c_0, 'c', 0);

%signal_c = nan(nt, length(sensorIndex));  % this should be a local variable

% this step should not be necessary if c_k is returned from the 
% solve_advection_diffusion_equation.m
%for n = 1 : nt
%  concentrationPath = fullfile(simDir,'Concentration',['Concentration_' num2str(n-1) '.mat']);
%  load(concentrationPath)
%  signal_c(n,:) = c(sensorIndex);
%end

f_obs = 0.5*sum(sum((c_k-c_star).^2));  % f_obs = 0.5*(Hc(x)-c*)'R^-1(Hc(x)-c*), R^-1 is assumed to be the identity matrix


initialSourcePath = fullfile(simDir, 'Source', 'Source_0.mat');
load(initialSourcePath,'x');
s_b = volumetric_emission_rate(x);  % initial candidate source distribution


p = reg_par;
f = f_obs + 0.5*p*((s_k-s_b)'*(B_inv*(s_k-s_b)));



% solve for the sensitivity of the objective function to the receptor
% observations 
working_dir = fullfile( sim_dir, 'Adjoint');
df_dx = solve_adjoint_advection_diffusion(c_k, t, H, working_dir, oper_dir);
%df_dx = CostFncGrad2(t, nNodes);

% add sensitivity of objective function to the regularization term
df_dx = dJ_dE + p*(B_inv*(s_k-s_b));  


if nPass == 0
  resultsPath = fullfile(simDir,'Source',['Source_', num2str(nPass), '.mat']);
  %load(fullfile(simDir,'Source',['Source_Correct.mat']),'signal_o')
  f = J;
 
  g = df_dx;

  s = E_approx;

  r  = compute_correlation_coefficient(c_k, c_star);
  r2 = compute_coefficient_of_determination(c_k, c_star);

	m = dot(spaceIntWeight, E_approx);
	m_norm = m/m_star;
  fprintf('| %s | %4d | %8.6e |  %8.6e | %9.6e | %8.6e | %8.6e |\n', datestr(now), nPass, f, max(abs(g)),r2,r, m_norm);

  save(resultsPath,'f','g','s','signal_c','r','r2','m','m_norm');
end

nPass = nPass + 1;