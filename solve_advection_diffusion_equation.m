function r_obs = solve_advection_diffusion_equation( s, t, H, ...
                    working_dir, operator_dir, theta, c_0 )
%solve_advection_diffusion_equation numerically models the transport 
% of a tracer from volumetric emission source and returns matrix of 
% the receptor observations.
%
% INPUTS:
%   s     vector of length P representing discretized volumetric 
%           source emission rate function with the units kg/m^3/s and
%           P is the number of nodes in the discretized domain.  s is 
%           assumed to be a constant emission source in time.
%
%   t     vector of length nt that represents the time steps in 
%           seconds. t can contain nonuniform time steps but is 
%           subject to the conditions:
%               t(1) = 0 and t(n) < t(n+1) for 1 <= n < nt.
%
%   H     R-by-P sparse observation matrix for computing the receptor 
%           observations from the discretized concentration field.  
%           Where R is the number of receptors in the domain.  The
%           observation matrix H is pre-computed with the function 
%           'compute_observation_matrix'.
%
%   working_dir   character string of the system path to the directory
%                   where the discretized concentration distributions
%                   will be saved.
%
%   operator_dir  character string of the system path to the directory
%                   where the finite element capacitance and stiffness
%                   matrices are saved.
% 
%   theta   optional time step integration weight factor for solving 
%             the transient evolution of the tracer concentration.
%
%               theta = 0     Forward difference method (Forward Euler
%                               method), 1st order accurate in time
%                               and conditionally stable.
%
%               theta = 0.5   Crank-Nicolson method, 2nd order 
%                               accurate in time and unconditionally 
%                               stable.
%
%               theta = 2/3   Galerkin method, unconditionally stable.
%                                 
%
%               theta = 1     Backward difference method (Backward 
%                               Euler method), 1st order accurate in
%                               time and unconditionally stable.
%
%               theta >= 0.5  unconditionally stable.
%
%             If no value is supplied for theta then the default value
%             of theta = 0.5 is used (Crank-Nicolson method, 2nd order 
%             accurate in time and unconditionally stable).
%
%   c_0     optional initial tracer concentration condition in kg/m^3.
%
%             1) If this parameter is not supplied then the initial  
%                 tracer concentration is computed from the steady  
%                 state of the tracer for the initial flow conditions.
%
%             2) If this parameter is a scalar then a uniform initial 
%                 tracer concentration of c_0 is assumed.
%
%             3) If this parameter is a vector of length P then c_0 
%                 represents the 
%
% OUTPUTS:
%   r_obs   nt-by-R matrix of rector observation in kg/m^3. 
%             r_obs(i, j) represents the tracer concentration measured
%             by the jth receptor at time t(i).
%
% EXAMPLES:
%   theta and c_0 are optional arguments.  Below are a few examples of
%   how 'solve_advection_diffusion_equation' can be called.
%
%   r_obs = solve_advection_diffusion_equation( s, t, H, ...
%                working_dir, operator_dir );
%
%     defaults invoked: theta = 0.5 and c_0 is the steady state tracer
%                       concentration distribution for the initial 
%                       flow conditions.
%
%
%   r_obs = solve_advection_diffusion_equation( s, t, H, ...
%                working_dir, operator_dir, theta );
%
%     defaults invoked: c_0 is the steady state tracer concentration
%                       distribution for the initial flow conditions.
%
%
%   r_obs = solve_advection_diffusion_equation( s, t, H, ...
%                working_dir, operator_dir, 7.06E-7 );
%
%     defaults invoked: c_0 is a uniform initial tracer concentration 
%                       of 7.06E-7 kg/m^3.
%
%
%   r_obs = solve_advection_diffusion_equation( s, t, H, ...
%                working_dir, operator_dir, 0.5, c_0 );

%-----------------------------------------------------------------------
% CHECK FOR VALID INPUT                                                |
%-----------------------------------------------------------------------
% Ensure the require number of arguments are supplied.
assert( 5 <= nargin && nargin <= 7 ); 

% Check s
assert( isvector(s) );                % ensure s is a vector
assert( isnumeric(s) );               % ensure s is numeric array
if size(s, 1) == 1 && size(s, 2) > 1  % ensure s is a column vector
  s = s';
end
nNodes = length(s);                   % number of nodes in the mesh

% Check t
assert( isvector(t) );            % ensure t is a vector
assert( isnumeric(t) );           % ensure t is numeric array
assert( t(1) = 0 );               % ensure t starts at 0
nt = length(t);
for n = 2 : nt                    % ensure t contains sequential time 
  assert( t(n-1) < t(n) );        %   steps
end

% Check H
assert( issparse(H) );            % ensure H is a sparse matrix
assert( isnumeric(H) );           % ensure H is a numeric array
assert( size(H,2) == size(s,1) ); % ensure dimensions of H and s agree

% Check working_dir
assert( ischar(working_dir) );  % ensure working_dir is a string
assert( exist( working_dir, 'dir') == 7 );

% Check operator_dir
assert( ischar(operator_dir) ); % ensure operator_dir is a string
assert( exist( operator_dir, 'dir') == 7 );

% If theta is supplied make sure it contains valid data.
if nargin >= 6                        % if theta is supplied
  assert( isscalar(theta) );          % ensure theta is a scalar
  assert( isnumeric(theta) );         % ensure theta is a number
  assert( 0 <= theta && theta <= 1 ); % ensure 0  <= theta <= 1
else
  theta = 0.5 % default value (Crank-Nicholson method)
end

% If c_0 is supplied make sure it contains valid data.  If c_0 is not 
% supplied then the initial condition is assumed to be the steady 
% state concentration distribution for the initial flow conditions.
if nargin == 7              % if c_0 is supplied
  assert( isnumeric(c_0) ); % ensure c_0 is numeric array
  assert( isvector(c_0) );  % ensure c_0 is a vector(or scalar 1-by-1)
  if isscalar(c_0)          % constant initial condition
    c_0 = c_0 * ones(nNodes, 1);
  else  % c_0 must be a vector, isvector(c_0) = true
    if size(c_0, 1) == 1 && size(c_0, 2) > 1  % ensure c_0 is a column 
      c_0 = c_0';                             %   vector
    end
    assert( length(c_0) == size(H,2));  % ensure dimensions of H and 
  end                                   %   c_0 agree
end
%---------------------------------------------------------------------

conc_label = 'Concentration_';
oper_label = 'Operators_';
nReceptors = size(H, 1);

% pre-allocate the receptor observation matrix
r_obs = nan(nt, nReceptors);


% load the initial operators 
file_num = generate_file_num(0, nt);
operator_path = fullfile(operator_dir, [oper_label file_num '.mat']);
load( operator_path, 'C', 'K' );
uPath = fullfile(working_dir, '..', 'U.mat');
load(uPath);
Cn = C;
Kn = K;

% Compute initial concentration distribution from a steady state
% solution of the source with in initial flow conditions.
if nargin < 7
  c_0 = U'*((U*Kn*U')\(U*Cn*s));  % compute the initial concentration
end
c = c_0;
clear c_0   % c_0 is no longer required
r_obs(:, 1) = (H*c)';   % compute the initial receptor observations
working_path = fullfile( working_dir, ...
                [conc_label, file_num, '.mat'] );
save(working_path, 'c');

% Compute the transient evolution of tracer concentration within the
% domain and the associated receptor observations
for n = 1 : nt - 1
  file_num = generate_file_num(n, nt);
  operator_path = fullfile(operator_dir, ...
                    [oper_label, file_num, '.mat'] );
  working_path = fullfile(working_dir, ...
                    [conc_label, file_num, '.mat'] );
  load(operator_path,'C','K');
  dt = t(n+1) - t(n);
  A = U*(C/dt+theta*K)*U';
  b = U*C*(c/dt-(1-theta)*(Cn\(Kn*c))+s);
  c = U'*(A\b);
  save(working_path, 'c');
  r_obs(:,n+1) = (H*c)';
  Cn = C;
  Kn = K;
end

%function receptor_observations = solve_advection_diffusion_equation(t,E,c_0,type,noise)
%global simDir operDir, H
%B = ComputeBoundaryMatrix(boundary,nNodes);
%if type == 'o'
%  o = c_0;
%  o = U'*((U*Kn*U')\(U*Cn*E));  % compute the initial concentration

%  concentrationPath = fullfile(simDir,'Observation','Observation_0.mat');
%  noisePath = fullfile(simDir,'Noise','Noise_0.mat');

%  o_error = CalcMeasurementNoise(o, noise); % compute the noise
%  save(concentrationPath,'o')
%  save(noisePath,'o_error')

%  receptor_observations(:,1) = (H*(o+o_error))';

%elseif type == 'c'
%  c = c_0;
%  c = U'*((U*Kn*U')\(U*Cn*E));
%  concentrationPath = fullfile(simDir,'Concentration','Concentration_0.mat');
%  save(concentrationPath,'c')
%  if nPass == 0
%    concentrationPath = fullfile(simDir,'Concentration_Initial_Guess','Concentration_0.mat');
%    save(concentrationPath,'c')
%  end

%else  % if type ~= 'o' or 'c' then display an error message
%  fprintf('ERROR:');
%end

%for n = 1 : nt - 1

%  operatorPath = fullfile(operDir, ['Operators_' num2str(n) '.mat']);
%  load(operatorPath,'C','K');
%  dt = t(n+1) - t(n);
%  A = U*(C/dt+theta*K)*U';

%  if type == 'o'
%    b = U*C*(o/dt-(1-theta)*(Cn\(Kn*o))+E);
%    o = U'*(A\b);

%    concentrationPath = fullfile(simDir,'Observation',['Observation_' num2str(n) '.mat']);
%    noisePath = fullfile(simDir,'Noise', ['Noise_' num2str(n) '.mat']);
%    o_error = CalcMeasurementNoise(o,noise);
%    save(concentrationPath,'o');
%    save(noisePath, 'o_error');

%  elseif type == 'c'
%    b = U*C*(c/dt-(1-theta)*(Cn\(Kn*c))+E);
%    c = U'*(A\b);

%    concentrationPath = fullfile(simDir,'Concentration',['Concentration_' num2str(n) '.mat']);

%    save(concentrationPath,'c');

%    if nPass == 0
%      concentrationPath = fullfile(simDir,'Concentration_Initial_Guess',['Concentration_' num2str(n) '.mat']);
%      save(concentrationPath,'c')
%    end
%  else
%    fprintf('ERROR');    
%  end

%  Cn = C;
%  Kn = K;
%end