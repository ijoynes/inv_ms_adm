function r_obs = solve_advection_diffusion_equation( s, t, H, ...
  working_dir, operator_dir, theta, c_0 )
%function receptor_observations = solve_advection_diffusion_equation(t,E,c_0,type,noise)
%global simDir operDir, H

%solve_advection_diffusion_equation numerically models the transport 
% of a tracer from volumetric emission source and returns matrix of 
% the receptor observations.
%
% INPUTS:
%   s     vector of discretized volumetric source emission rate 
%           function.
%
%   t     vector of 
%

%
%   H
%
%   working_dir
%
%   operator_dir
% 
%   theta
%
%   c_0   vector of 
%
% OUTPUTS:
%   

%---------------------------------------------------------------------
% CHECK FOR VALID INPUT
%---------------------------------------------------------------------
% Ensure the require number of arguments are supplied.
assert( 5 <= nargin && nargin <= 7 ); 

% Check s
assert( isvector(s) );                % ensure s is a vector
assert( isnumeric(s) );               % ensure s is numeric array
if size(s, 1) == 1 && size(s, 2) > 1  % ensure s is a column vector
  s = s';
end

% Check t
assert( isvector(t) );          % ensure t is a vector
assert( isnumeric(t) );         % ensure t is numeric array

% Check H
assert( issparse(H) );          % ensure H is a sparse matrix
assert( isnumeric(H) );         % ensure H is a numeric array

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
% supplied then the inital condition is assumed to be the steady state
% concentration distribution for the initial flow conditions.
if nargin == 7              % if c_0 is supplied
  assert( isvector(c_0) );  % ensure c_0 is a vector
  assert( isnumeric(c_0) ); % ensure c_0 is numeric array
  if size(c_0, 1) == 1 && size(c_0, 2) > 1  % ensure c_0 is a column 
    c_0 = c_0';                             %   vector
  end
end
%---------------------------------------------------------------------





%load(fullfile(simDir,'dirSettings.mat'));
%if type == 'c'
%load(passPath);
%end

%theta = 0.5;
nt = length(t);

for n = 2 : nt
  assert( t(n-1) < t(n) );
end
r_obs = nan(nt,nReceptors)


%B = ComputeBoundaryMatrix(boundary,nNodes);
operatorPath = fullfile(operDir, ['Operators_0.mat']);
load(operatorPath,'C','K');
uPath = fullfile(simDir,'U.mat');
load(uPath);
Cn = C;
Kn = K;


% Compute initial concentration distribution from a steady state
% solution of the source with in initial flow conditions.
if type == 'o'
  o = c_0;
  o = U'*((U*Kn*U')\(U*Cn*E));  % compute the initial concentration

  concentrationPath = fullfile(simDir,'Observation','Observation_0.mat');
  noisePath = fullfile(simDir,'Noise','Noise_0.mat');

  o_error = CalcMeasurementNoise(o, noise); % compute the noise
  save(concentrationPath,'o')
  save(noisePath,'o_error')

  receptor_observations(:,1) = (H*(o+o_error))';

elseif type == 'c'
  c = c_0;
  c = U'*((U*Kn*U')\(U*Cn*E));
  concentrationPath = fullfile(simDir,'Concentration','Concentration_0.mat');
  save(concentrationPath,'c')
  if nPass == 0
    concentrationPath = fullfile(simDir,'Concentration_Initial_Guess','Concentration_0.mat');
    save(concentrationPath,'c')
  end

else  % if type ~= 'o' or 'c' then display an error message
  fprintf('ERROR:');
end

for n = 1 : nt - 1

  operatorPath = fullfile(operDir, ['Operators_' num2str(n) '.mat']);
  load(operatorPath,'C','K');
  dt = t(n+1) - t(n);
  A = U*(C/dt+theta*K)*U';

  if type == 'o'
    b = U*C*(o/dt-(1-theta)*(Cn\(Kn*o))+E);
    o = U'*(A\b);

    concentrationPath = fullfile(simDir,'Observation',['Observation_' num2str(n) '.mat']);
    noisePath = fullfile(simDir,'Noise', ['Noise_' num2str(n) '.mat']);
    o_error = CalcMeasurementNoise(o,noise);
    save(concentrationPath,'o');
    save(noisePath, 'o_error');

  elseif type == 'c'
    b = U*C*(c/dt-(1-theta)*(Cn\(Kn*c))+E);
    c = U'*(A\b);

    concentrationPath = fullfile(simDir,'Concentration',['Concentration_' num2str(n) '.mat']);

    save(concentrationPath,'c');

    if nPass == 0
      concentrationPath = fullfile(simDir,'Concentration_Initial_Guess',['Concentration_' num2str(n) '.mat']);
      save(concentrationPath,'c')
    end
  else
    fprintf('ERROR');    
  end

  Cn = C;
  Kn = K;
end

