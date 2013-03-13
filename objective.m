function [f,df_dx]=objective(x_k)

global simDir operDir domainPath paramPath sensorPath sourcePath ...
        passPath tMax noise reg_par signal_o signal_c ...
        spaceIntWeight m_star nPass c_star B_inv


s_k = volumetric_emission_rate(x_k);
%load(paramPath)

c_k = solve_advection_diffusion_equation(s_k)
%SolveConcentrationTransport(t, s_k, c_0, 'c', 0);

%signal_c = nan(nt, length(sensorIndex));  % this should be a local variable

% this step should not be necessary if c_k is returned from the 
% solve_advection_diffusion_equation.m
%for n = 1 : nt
%  concentrationPath = fullfile(simDir,'Concentration',['Concentration_' num2str(n-1) '.mat']);
%  load(concentrationPath)
%  signal_c(n,:) = c(sensorIndex);
%end

f_obs = 0.5*sum(sum((c_k-c_star).^2));


p = reg_par;

% should B_inv be finite element capacitance matrix instead of a diagonal matrix?
%B_inv = spdiags(spaceIntWeight,0,length(E_approx),length(E_approx));


initialSourcePath = fullfile(simDir, 'Source', 'Source_0.mat');
load(initialSourcePath,'x');
s_b = volumetric_emission_rate(x);



f = f_obs + 0.5*p*((s_k-s_b)'*(B_inv*(s_k-s_b)));

% solve for the sensitivity of the objective function to the receptor
% observations 
% df_dx = solve_adjoint_advection_diffusion;
df_dx = CostFncGrad2(t, nNodes);

% add sensitivity of objective function to the regularization term
df_dx = dJ_dE + p*(B_inv*(s_k-s_b));  

if nPass == 0
  resultsPath = fullfile(simDir,'Source',['Source_', num2str(nPass), '.mat']);
  %load(fullfile(simDir,'Source',['Source_Correct.mat']),'signal_o')
  f = J;
  g = dJ_dE;
  s = E_approx;
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

   ybar = mean(mean(signal_o));
   r2 = 1 - sum(sum((signal_o-signal_c).^2))/sum(sum((signal_o-ybar).^2));
	m = dot(spaceIntWeight, E_approx);
	m_norm = m/m_star;
  fprintf('| %s | %4d | %8.6e |  %8.6e | %9.6e | %8.6e | %8.6e |\n', datestr(now), nPass, f, max(abs(g)),r2,r, m_norm);

  save(resultsPath,'f','g','s','signal_c','r','r2','m','m_norm');
end

nPass = nPass + 1;