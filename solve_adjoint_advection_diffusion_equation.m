function df_dx = solve_adjoint_advection_diffusion_equation(c_k, c_star, t, H, working_dir, oper_dir, theta, c_0)
%solve_adjoint_advection_diffusion_equation computes the gradient of 
% the observational component of objective function with respect to
% the volumetric source emission rate parameters x with the discrete 
% adjoint method.

% Ensure the require number of arguments are supplied.
assert( 6 <= nargin && nargin <= 8 ); 

% Check c_k
assert( isvector(c_k) );                % ensure s is a vector
assert( isnumeric(c_k) );               % ensure s is numeric array

% Check c_star
assert( isvector(c_star) );                % ensure s is a vector
assert( isnumeric(c_star) );               % ensure s is numeric array


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
assert( size(H, 1) == size(c_k, 2) );  % ensure dimensions of H and c_k agree
assert( size(H, 1) == size(c_star, 2) );  % ensure dimensions of H and c_star agree

% Check working_dir
assert( ischar(working_dir) );  % ensure working_dir is a string
assert( exist( working_dir, 'dir') == 7 );

% Check operator_dir
assert( ischar(operator_dir) ); % ensure operator_dir is a string
assert( exist( operator_dir, 'dir') == 7 );

% If theta is supplied make sure it contains valid data.
if nargin >= 7                        % if theta is supplied
  assert( isscalar(theta) );          % ensure theta is a scalar
  assert( isnumeric(theta) );         % ensure theta is a number
  assert( 0 <= theta && theta <= 1 ); % ensure 0  <= theta <= 1
else
  theta = 0.5 % default value (Crank-Nicholson method)
end

% If c_0 is supplied make sure it contains valid data.  If c_0 is not 
% supplied then the initial condition is assumed to be the steady 
% state concentration distribution for the initial flow conditions.
if nargin == 8              % if c_0 is supplied
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

nNodes = size(H, 2);
nx = nNodes;

df_dx = zeros(nx,1);
y = H'*(c_k(nt,:)-c_star(nt,:))';


uPath = fullfile(simDir,'U.mat');
load(uPath);

% load operators for the last time step
operatorPath = fullfile(operDir, ['Operators_' num2str(nt-1) '.mat']);
load(operatorPath,'C','K');

% copy operators so that they can still be accessed in the next time step.
Cn = C;
Kn = K;

% compute the backwards in time evolution of the adjoint variable and its contribution to the gradient of the objective function.
for n = nt-1 :-1: 1

  operatorPath = fullfile(operDir, ['Operators_' num2str(n-1) '.mat']);
  load(operatorPath,'C','K');

  dt = t(n+1) - t(n);

  temp = Cn'*(U'*((U*(Cn'/dt+theta*Kn')*U')\(U*y)));
  df_dx = df_dx + temp;
  if nPass == 0
    gradPath = fullfile(simDir, 'Gradient', ['Gradient_' num2str(n) '.mat']);
    save(gradPath,'df_dx');
  end

  y = temp/dt -(1-theta)*K'*(C'\temp);
  y = y + H'*(c_k(n,:) - c_star(n,:))';
  if nPass == 0
    adjointPath = fullfile(simDir,'Adjoint', ['Adjoint_' num2str(n-1) '.mat']);
    save(adjointPath,'y');
  end

  % copy operators so that they can still be accessed in the next time step.
  Cn = C;
  Kn = K;

end

% if no initial conditions for the tracer concentration field are supplied then assume that
% the initial conditions of the tracer concentration field are dependent on the
% volumetric source emission rate then compute the contribution of the initial conditions 
% to the objective function gradient.
if nargin < 8
  temp = Cn'*U'*((U*Kn'*U')\(U*y));
  df_dx = df_dx + temp;
end

if nPass == 0
  gradPath = fullfile(simDir, 'Gradient', 'Gradient_0.mat');
  save(gradPath,'df_dx');
end

%###################################################################################################################
%###################################################################################################################
%###################################################################################################################

%function dJ_dE = CostFncGrad2(t, nNodes)

global simDir operDir tMax signal_o signal_c

load(fullfile(simDir,'dirSettings.mat'));
load(passPath);
nt = length(t);
theta = 0.5;

   %%% observationPath = fullfile(simDir,'Observation',['Observation_' num2str(nt-1) '.mat']);
   %%% concentrationPath = fullfile(simDir,'Concentration',['Concentration_' num2str(nt-1) '.mat']);
   %%% noisePath = fullfile(simDir,'Noise',['Noise_' num2str(nt-1) '.mat']);
   %%% load(observationPath)
   %%% load(concentrationPath)
   %%% load(noisePath)
   load(sensorPath)
y = zeros(nNodes,1);   
%%% y(sensorIndex) = c(sensorIndex)-(o(sensorIndex)+o_error(sensorIndex));
y(sensorIndex) = signal_c(nt,:)-signal_o(nt,:);
if nPass == 0
adjointPath = fullfile(simDir,'Adjoint', ['Adjoint_' num2str(nt-1) '.mat']);
save(adjointPath,'y');
end
dJ_dE = zeros(nNodes,1);

uPath = fullfile(simDir,'U.mat');
load(uPath);
operatorPath = fullfile(operDir, ['Operators_' num2str(nt-1) '.mat']);
   load(operatorPath,'C','K');
   Cn = C;
   Kn = K;
for n = nt-1 :-1: 1   
   
	%fprintf('%s %d\n','y @ t =', n);
   	operatorPath = fullfile(operDir, ['Operators_' num2str(n-1) '.mat']);
   	load(operatorPath,'C','K');

   
   dt = t(n+1) - t(n);

   temp = Cn'*(U'*((U*(Cn'/dt+theta*Kn')*U')\(U*y)));
   
   dJ_dE = dJ_dE + temp;
if nPass == 0
 gradPath = fullfile(simDir, 'Gradient', ['Gradient_' num2str(n) '.mat']);
save(gradPath,'dJ_dE');
end  
   %if n > 1
		   
   %%% observationPath = fullfile(simDir,'Observation',['Observation_' num2str(n-1) '.mat']);
   %%% concentrationPath = fullfile(simDir,'Concentration',['Concentration_' num2str(n-1) '.mat']);
   %%% noisePath = fullfile(simDir,'Noise',['Noise_' num2str(nt-1) '.mat']);

   %%% load(observationPath)
   %%% load(concentrationPath)
   %%% load(noisePath)
   y = temp/dt -(1-theta)*K'*(C'\temp);
   %%% y(sensorIndex) = y(sensorIndex) +c(sensorIndex)-(o(sensorIndex)+o_error(sensorIndex));
   y(sensorIndex) = y(sensorIndex) + (signal_c(n,:) - signal_o(n,:))';
   if nPass == 0
adjointPath = fullfile(simDir,'Adjoint', ['Adjoint_' num2str(n-1) '.mat']);
save(adjointPath,'y');
end
   Cn = C;
   Kn = K;
   %end

end

temp = Cn'*U'*((U*Kn'*U')\(U*y));
dJ_dE = dJ_dE + temp;
if nPass == 0
gradPath = fullfile(simDir, 'Gradient', 'Gradient_0.mat');
save(gradPath,'dJ_dE');
end
