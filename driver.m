function exitflag = driver(controlsFilePath)
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
%   include: simDir, operDir, tMax, noise, maxIter, reg_par, factr, 
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
%
%
% VARIABLES
%   simDir    Directory path the the simulation results.
%   operDir   Directory path the stiffness and capacitance matrices.
%   tMax      The max time of the simulation.
%   noise     The slope of the line which relates the observation 
%               concentration to the standard deviation of the 
%               Gaussian distributed noise (0.1 = 10% noise).
%   maxIter   Maximum number of iterations.
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
%   spaceIntWeight  A weight matrix for the regularization term of the 
%                     objective function.  This weight matrix 
%                     penalizes emission rates from larger elements 
%                     more than smaller elements.
%   m_star    The total emission rate of the known source in kg/s.


% Declare the global variables. 
% The LBFGS function does not allow for additional variables to be 
% passed to the objective function (or any other function the 
% objective function may call).  Critical variables which are 
% determined before the objective function is called are declared as 
% global variables.
%%%global simDir operDir domainPath paramPath sensorPath sourcePath ...
%%%        passPath tMax noise reg_par signal_o signal_c ...
%%%        spaceIntWeight m_star nPass B_inv

global sim_dir oper_dir iter_dir conc_dir adj_dir iter_label conc_label adj_label grad_label H c_star save_flag

obs_dir = fullfile(sim_dir, 'Observation');

iter_label = 'Source_';
conc_label = 'Concentration_';
adj_label  = 'Adjoint_';
grad_label = 'Gradient_';



% Write the simulation start date and time to the log file.
fprintf('%s\n', datestr(now) );

% Set the control file to read only.
fileattrib(controlsFilePath,'-w'); 

% Read the run parameters from the  control file and write their 
% values to the log file.
[sim_dir, oper_dir, domainPath, paramPath, sensorPath, sourcePath, ... 
  passPath, tMax, dt, noise, maxIter, reg_par ...
  factr, pgtol, m, iprint, save_flag] = readControls(controlsFilePath)

% Declare the file paths the critical simulation files
%domainPath = fullfile(simDir,'Domain.mat');
%paramPath = fullfile(simDir,'Parameters.mat');
%sensorPath = fullfile(simDir,'Sensors.mat');
%sourcePath = fullfile(simDir,'Source.mat');
%passPath = fullfile(simDir,'pass.mat');

% Save these file paths for later reference
%save(fullfile(simDir,'dirSettings.mat'), 'domainPath','paramPath', ...
%  'sensorPath','sourcePath','passPath');

% Make directories to store data generated by the atmospheric 
% transport model, its adjoint and the LBFGS.
mkdir(fullfile(sim_dir,'Source'));
mkdir(fullfile(sim_dir,'Concentration'));
mkdir(obs_dir);
mkdir(fullfile(sim_dir,'Concentration_Initial_Guess'));
mkdir(fullfile(sim_dir,'Adjoint'));
mkdir(fullfile(sim_dir,'Gradient'));
mkdir(fullfile(sim_dir,'Noise'));

% Load the mesh, sensor placement and emission source details
load(domainPath);
load(sensorPath);
load(sourcePath);

%---------------------------------------------------------------------
% The integration of a scalar over the entire domain can be 
% approximated numerically by summing the integral of the scalar over
% each element.  The integral of the scalar over an element can be 
% approximated summing the integrals of the product of finite element 
% shape functions and the nodal value of the scalar field.  This 
% integration operation can be streamlined by pre-computing a spatial 
% integration weight vector (spaceIntWeight in this case), which is a
% vector of weights based on the area (or volume) of elements that a 
% given mesh node shares, and indicates the influence that a given 
% node has on the spatial integration of the scalar field.  With this 
% precomputed spatial integration weight vector the spatial 
% integration of any scalar field (represented as a vector of nodal 
% values) is computed with the dot product between these to vectors.
% For instance the total emission rate of a source volumetric emission 
% rates would be computed as follows:
%
% total_emission_rate = dot(spaceIntWeight, s);
%
% Construct the spatial integration weight vector
spaceIntWeight = zeros(nNodes,1);
for i = 1 : nTris
  spaceIntWeight(tri(i,:)) = spaceIntWeight(tri(i,:)) + ...
                             det([ones(3,1) xy(tri(i,:),:) ] );
end
spaceIntWeight = spaceIntWeight/6;
%---------------------------------------------------------------------
% Construct the inverse of the covariance matrix of the estimated 
% background error.  This matrix will numerically integrate the square
% of the residual between the discrete representation of the candidate 
% source distribution (s_k) and the discrete representation of the 
% assumed prior source distribution (s_b), such that this relationship
% is true.
%
% int((s_k-s_b)^2, Omega) = (s_k-s_b)'*B^-1*(s_k-s_b)
Bv = nan(nTris*3*3,1);
row = nan(nTris*3*3,1);
col = nan(nTris*3*3,1);
index = 0;
for i = 1 : nTris
    Be = [ 2 1 1; 1 2 1; 1 1 2] * det([ones(3,1), xy(tri(i,:),:)]);
    for j = 1 : 3
        for k = 1 : 3
            index = index + 1;
            Bv(index) = Be(j,k);
            row(index) = tri(i,j);
            col(index) = tri(i,k);
        end
    end
end
Bv = Bv/24;
B_inv = sparse(row,col,Bv,nNodes,nNodes);
clear Bv Be row col index i j k

% Move sensors to the closest mesh nodes
% This step could be ignored if an observation matrix H was 
% constructed to perform a weighted average of the concentration at 
% the nodes of the element which contains the receptor.  The weights 
% would be the values of the finite element shape functions at the 
% location of the receptor. 
H = compute_observation_matrix(tri, xy, receptor_xy);

% Determine discrete representation of the emission source
% Again, this function moves the sources to the nearest mesh nodes.
% The volumetric emission rate of each source is scaled by the inverse
% of the size of the element such that a cluster of large or small 
% elements would have the same total emission rate.
s_star = place_sources(tri, xy, source_xy, source_m);

% Compute the total emission rate in kg/s of the known source by 
% integrating the volumetric source emission rate over the spatial 
% domain.  This integration is approximated numerically with the 
% spatial integration weight vector.  
m_star = dot(spaceIntWeight,s_star);

% Reorient the domain such that the lower left corner has the 
% coordinate (0,0)
xy = xy - ones(nNodes,1)*min(xy);

% Create the time step vector
t = (0:dt:tMax)';
nt = length(t);

c_0 = zeros(nNodes,1);

% --------------------------------------------------------------------
% If an initial source guess is present in the file path
%   simDir/Source/Source_0.mat
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

% Compute the evolution of the plume from the known source
%SolveConcentrationTransport(t,E,c_0,'o',noise);
c_star = solve_advection_diffusion_equation(s_star, t, H, true, obs_dir, oper_dir);
c_star = add_observation_noise(c_star, noise);

% Pre-allocate space for the synthetic and candidate receptor 
% observations
%signal_o=nan(nt,length(sensorIndex));
%signal_c=nan(nt,length(sensorIndex));

% Extract the synthetic receptor observations from the known source
%for i = 1 : nt
%  load(fullfile(simDir, 'Observation', ...
%                  ['Observation_' int2str(i-1) '.mat']),'o')
%  if noise > 0 || noise == -1 % add noise to the receptor observations
%    load(fullfile(simDir, 'Noise', ...
%                    ['Noise_' int2str(i-1) '.mat']), 'o_error')
%    signal_o(i,:) = o(sensorIndex) + o_error(sensorIndex);
%  else
%    signal_o(i,:) = o(sensorIndex);
%  end
%end

% Append the synthetic receptor observations to the known source file
save(fullfile(sim_dir,'Source','Source_Correct.mat'), ...
  'E','signal_o', 'm_star');

% These should be global variables
%%%%save(paramPath,'t', 'c_0', 'boundary', 'sensorIndex','nNodes', 'nt','tri','xy','E')
% This should also be a global variable
nPass = 0;

% Set the parameters for the operation of the L-BFGS-B routines
fn = @objective_function; % objective function handle

x0 = E_0;                 % initial parameter estimate
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
                     'maxits', maxIter, ...
                     'factr', factr, ...
                     'pgtol', pgtol, ...
                     'm', m, ...
                     'cb', @callback);
           

[x,fx,exitflag,userdata] = lbfgs(fn,x0, lb, ub, nbd, opts);


save(fullfile(simDir,'Source','Optimum_Source.mat'),'x','fx');
save(fullfile(simDir,'exit_dump.mat'));
fprintf('Done: %s\n', datestr(now));
