function [exit_flag, x, f, user_data] = driver(controls_file_path)
%DRIVER setups up and runs the Inverse Micro-Scale Atmospheric 
% Dispersion Model to find the optimum volumetric source to profile 
% to recreate synthetic receptor observations from a known source
%
% Last Modified:
%   March 04, 2013 by Ian Joynes
%
% PROCEDURE
%   This function reads a control file which gives the values of the 
%   Inverse Micro-Scale Atmospheric Dispersion Model parameters which 
%   include: sim_dir, operDir, tMax, noise, max_iters, reg_par, factr, 
%   pgtol, m and iprint.  The function then loads the 2D unstructured 
%   triangular mesh, the location of the receptors and the location 
%   and emission rate of the known sources.  The concentration of the 
%   methane released from the known sources is determined by 
%   numerically solving the advection-diffusion equation.  The 
%   synthetic receptor observation are extracted from the 
%   concentration field.  With the synthetic receptor observations, 
%   source reconstruction begins by calling the general lbfgs driver.
%   This function iteratively improves the fitness of candidate 
%   volumetric emission source profiles by minimizing an objective 
%   function based on the misfit between the synthetic receptor 
%   observations and the modelled observations from a candidate 
%   volumetric emission source profile.
%
% REQUIRED FILES
%
% INPUTS
%   controls_file_path    file path to a text file containing the values
%                           of the parameters
%
% OUTPUS:
%   x           final value
%   f           fn(x)
%   user_data   user data from callback function
%   exit_flag   0: algorithm converged or max iterations reached
%               1: abnormal termination
%               2: error
%               3: unknown cause of termination
%
% VARIABLES
%   sim_dir   Directory path the the simulation results.
%   operDir   Directory path the stiffness and capacitance matrices.
%   tMax      The max time of the simulation.
%   noise     The slope of the line which relates the observation 
%               concentration to the standard deviation of the 
%               Gaussian distributed noise (0.1 = 10% noise).
%   max_iters   Maximum number of iterations.
%   reg_par   Regularization parameter (reg_par >= 0).
%   signal_o  Matrix representing the synthetic receptor observations
%               each row represents a sample time, and each column 
%               represents an individual receptor.  Measurements are 
%               reported in the units kg/m^3.
%   signal_c  Matrix representing the receptor observations for a 
%               candidate volumetric source profile each row 
%               represents a sample time, and each column represents 
%               an individual receptor.  Measurements are reported in 
%                the units kg/m^3.
%   space_int_wgt  A weight matrix for the regularization term of the 
%                     objective function.  This weight matrix 
%                     penalizes emission rates from larger elements 
%                     more than smaller elements.
%   m_star    The total emission rate of the known source in kg/s.


% Declare the global variables. 
%   The LBFGS function does not allow for additional variables to be 
%   passed to the objective function or callback function.  Critical 
%   variables which are required to pass between the driver, objective 
%   function and callback function are declared as 
%   global variables.
global_vars
fprintf('Start of Inverse Micro-Scale Atmospheric Dispersion Model: %s\n', datestr(now) );
assert(exist(controls_file_path,'file')==2, 'The control file "%s" was not found.', controls_file_path);
exit_flag = 2;  % initialize to the error condition

% Write the simulation start date and time to the log file.


% Read the run parameters from the  control file and write their 
% values to the log file.
fprintf('\n');
fprintf('----------------------------------------------------------\n');
fprintf('|               Control parameter values:                |\n');
fprintf('----------------------------------------------------------\n');
fprintf('controls_file_path = %s\n', controls_file_path);
[sim_dir, oper_dir, domainPath, paramPath, sensorPath, sourcePath, ... 
  passPath, tMax, dt, noise, max_iters, reg_par ...
  factr, pgtol, m, iprint, save_flag] = readControls(controls_file_path)
fprintf('----------------------------------------------------------\n');
fprintf('\n');

obs_dir = fullfile(sim_dir, 'Observation');
conc_dir = fullfile(sim_dir, 'Concentration');
adj_dir = fullfile(sim_dir, 'Adjoint');
iter_dir = fullfile(sim_dir, 'Source');

iter_label = 'Source_';
conc_label = 'Concentration_';
adj_label  = 'Adjoint_';
grad_label = 'Gradient_';

% Declare the file paths the critical simulation files
%domainPath = fullfile(sim_dir,'Domain.mat');
%paramPath = fullfile(sim_dir,'Parameters.mat');
%sensorPath = fullfile(sim_dir,'Sensors.mat');
%sourcePath = fullfile(sim_dir,'Source.mat');
%passPath = fullfile(sim_dir,'pass.mat');

% Save these file paths for later reference
%save(fullfile(sim_dir,'dirSettings.mat'), 'domainPath','paramPath', ...
%  'sensorPath','sourcePath','passPath');

% Make directories to store data generated by the atmospheric 
% transport model, its adjoint and the LBFGS.
mkdir(iter_dir);
mkdir(conc_dir);
mkdir(obs_dir);
mkdir(fullfile(sim_dir,'Concentration_Initial_Guess'));
mkdir(adj_dir);
mkdir(fullfile(sim_dir,'Gradient'));
mkdir(fullfile(sim_dir,'Noise'));

% Load the mesh, sensor placement and emission source details
load(domainPath);
load(sensorPath);
load(sourcePath);

gd_idx = unique(boundary(:));
gd_val = zeros(length(gd_idx),1);

%-----------------------------------------------------------------------
% The integration of a scalar over the entire domain can be 
% approximated numerically by summing the integral of the scalar over
% each element.  The integral of the scalar over an element can be 
% approximated summing the integrals of the product of finite element 
% shape functions and the nodal value of the scalar field.  This 
% integration operation can be streamlined by pre-computing a spatial 
% integration weight vector (space_int_wgt in this case), which is a
% vector of weights based on the area (or volume) of elements that a 
% given mesh node shares, and indicates the influence that a given 
% node has on the spatial integration of the scalar field.  With this 
% precomputed spatial integration weight vector the spatial 
% integration of any scalar field (represented as a vector of nodal 
% values) is computed with the dot product between these to vectors.
% For instance the total emission rate of a source volumetric emission 
% rates would be computed as follows:
%
% total_emission_rate = dot(space_int_wgt, s);
%
% Construct the spatial integration weight vector
fprintf('\nComputing spatial integration weight vector...\n\n');
space_int_wgt = compute_spatial_integration_weight_vector(tri,xy);


%-----------------------------------------------------------------------
% Construct the inverse of the covariance matrix of the estimated 
% background error.  This matrix will numerically integrate the square
% of the residual between the discrete representation of the candidate 
% source distribution (s_k) and the discrete representation of the 
% assumed prior source distribution (s_b), such that this relationship
% is true.
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

%-----------------------------------------------------------------------
% Compute the total emission rate in kg/s of the known source by 
% integrating the volumetric source emission rate over the spatial 
% domain.  This integration is approximated numerically with the 
% spatial integration weight vector.  
m_star = dot(space_int_wgt, s_star);
%-----------------------------------------------------------------------

% Reorient the domain such that the lower left corner has the 
% coordinate (0,0)
xy = move_domain_origin(xy);

% Create the time step vector
t = (0:dt:tMax)';
nt = length(t);


% ----------------------------------------------------------------------
% If an initial source guess is present in the file path
%   sim_dir/Source/Source_0.mat
% then load that initial guess.  However if no initial guess is present
% assume an all zero initial guess and save it.
file_num = generate_file_num(0, nt);
if exist(fullfile(iter_dir,[iter_label file_num '.mat'])) == 2
  load(fullfile(iter_dir,[iter_label file_num '.mat']), 'x');
  s = volumetric_emission_rate(x);
  save(fullfile(iter_dir,[iter_label file_num '.mat']), 's', '-append');
else
  x = zeros(nNodes,1);
  s = volumetric_emission_rate(x);
  save(fullfile(iter_dir,[iter_label file_num '.mat']), 'x', 's');
end

% Save the volumetric source emission rate profile of the known source 
% with its total emission rate.
save(fullfile(iter_dir, [iter_label 'Correct.mat']), 's_star','m_star');


% Compute the synthetic receptor observation from the known source.
fprintf('\nComputing synthetic receptor observations...\n\n');
c_star_without_noise = solve_advection_diffusion_equation(s_star, t, H, save_flag, obs_dir, oper_dir);

% Add Gaussian distributed noise to the synthetic receptor observations.
c_star = add_observation_noise(c_star_without_noise, noise);

% Append the synthetic receptor observations to the known source file.
save(fullfile(iter_dir,[iter_label 'Correct.mat']),'c_star', 'c_star_without_noise', '-append');

% These should be global variables
%%%%save(paramPath,'t', 'c_0', 'boundary', 'sensorIndex','nNodes', 'nt','tri','xy','E')
% This should also be a global variable
nPass = 0;

% Set the parameters for the operation of the L-BFGS-B routines
fn = @objective;          % objective function handle

x0 = x;                   % initial parameter estimate
                          %   In this case the source parameters 
                          %   represent the nodal values of the 
                          %   volumetric source emission rate.
                          
lb = zeros(nNodes,1);     % lower bound on parameter search space
                          %   Setting this lower bound to zero ensures
                          %   that unphysical negative sources (sinks) 
                          %   are not predicted. This will need to be 
                          %   adjusted in future implementations if 
                          %   the source parameters x represent 
                          %   quantities other than volumetric 
                          %   emission rate (e.g. source location). 

ub = Inf*ones(nNodes,1);  % upper bound on parameter search space
                          %   Since the is no upper bound on the 
                          %   volumetric emission rate, this vector is 
                          %   not used by the L-BFGS-B routines and 
                          %   can assume any value.  A value of +ve 
                          %   infinity is assigned for clarity on the 
                          %   intent of the use of the upper bound

nbd = ones(nNodes,1);     % the search space bound types
                          % nbd(i)=0 if x(i) is unbounded,
                          %        1 if x(i) has only a lower bound
                          %        2 if x(i) has lower & upper bounds
                          %        3 if x(i) has only an upper bound 

% Setup the L-BFGS-B options
%   iprint  controls the frequency and type of output
%   maxits  maximum number of iterations
%   factr   exit criteria for small changes to the objective function
%   pgtol   exit criteria for small changes to the gradient
%   m       number of histories to use to approximate the Hessian
%   cb      callback function handle
opts = lbfgs_options('iprint', iprint, ...
                     'maxits', max_iters, ...
                     'factr', factr, ...
                     'pgtol', pgtol, ...
                     'm', m, ...
                     'cb', @callback);
           

% Set the assumed prior source distribution to the same as the initial 
% candidate source distribution.
s_b = volumetric_emission_rate(x);
[x, f, exit_flag, user_data] = lbfgs(fn, x0, lb, ub, nbd, opts);

save(fullfile(sim_dir,'Source','Optimum_Source.mat'),'x','f');
save(fullfile(sim_dir,'exit_dump.mat'));
fprintf('Done: %s\n', datestr(now));
