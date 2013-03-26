global_vars
controls_file_path = './controls_run_dummy.ini';


idx = 10000;

[sim_dir, oper_dir, domainPath, paramPath, sensorPath, sourcePath, ... 
  passPath, tMax, dt, noise, max_iters, reg_par ...
  factr, pgtol, m, iprint, save_flag] = readControls(controls_file_path)

load(domainPath);
load(sensorPath);
load(sourcePath);


x = zeros(nNodes,1);

save_flag = false;
obs_dir = '.';
fprintf('\nComputing inverse of source error covariance matrix...\n\n');
B_inv = compute_inverse_source_error_covariance_matrix(tri,xy);

%-----------------------------------------------------------------------
% Move sensors to the closest mesh nodes                               |
% This step could be ignored if an observation matrix H was            |
% constructed to perform a weighted average of the concentration at    |
% the nodes of the element which contains the receptor.  The weights   |
% would be the values of the finite element shape functions at the     |
% location of the receptor.                                            |
fprintf('\nComputing observation matrix...\n\n');
H = compute_observation_matrix(tri, xy, receptor_xy); %                |
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% Determine discrete representation of the emission source
% Again, this function moves the sources to the nearest mesh nodes.
% The volumetric emission rate of each source is scaled by the inverse
% of the size of the element such that a cluster of large or small 
% elements would have the same total emission rate.
fprintf('\nComputing discretized volumetric emission rate of known source...\n\n');
s_star = place_sources(tri, xy, source_xy, source_m);
%-----------------------------------------------------------------------

% Create the time step vector
t = (0:dt:tMax)';
nt = length(t);


fprintf('\nComputing synthetic receptor observations...\n\n');
c_star_without_noise = solve_advection_diffusion_equation(s_star, t, H, save_flag, obs_dir, oper_dir);

[f,g] = objective(x);

%dx = (eps)^(1/3);
%xp1 = x;
%xp1(idx) = x(idx) + dx;
%[fp1] = objective(xp1);

%xn1 = x;
%xn1(idx) = x(idx) - dx;
%[fn1] = objective(xn1);

%g(idx)
%(fp1-fn1)/(2*dx)

dx = sqrt(eps);
xcp1 = x;
xcp1 = x(idx)-i*dx;
[fcp1] = objective(xcp1);

imag(fcp1)/dx