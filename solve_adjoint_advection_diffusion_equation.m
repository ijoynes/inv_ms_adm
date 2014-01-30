function df_dx = solve_adjoint_advection_diffusion_equation(c_k, c_star, t, H, save_flag, working_dir, oper_dir, time_offset, nt_max, theta, c_0)
%solve_adjoint_advection_diffusion_equation computes the gradient of 
% the observational component of objective function with respect to
% the volumetric source emission rate parameters x with the discrete 
% adjoint method.


%-----------------------------------------------------------------------
% ########################## VARIABLE NOTES ########################## |
%-----------------------------------------------------------------------
% adj_label   string for the file label of the saved discrete adjoint  |
%               variable states                                        |
% adj_path    string for the file path to save a .mat file containing  |
%               the values of the discrete adjoint variable            |
% C           sparse finite element capacitance matrix saved in the    |
%               operator file                                          |
% Cn          sparse finite element capacitance matrix for the nth     |
%                time step                                             |
% dt          time step size                                           |
% file_num    string of digits containing a file number                |
% grad_label  string for the file label of the saved objective         |
%               function gradient states                               |
% grad_path   string for the file path to save a .mat file containing  |
%               the current state of the objective function gradient   |
% K           sparse finite element stiffness matrix saved in the      |
%               operator file                                          |
% Kn          sparse finite element stiffness matrix for the nth time  |
%               step                                                   |
% n           time step index                                          |
% nt          number of time steps                                     |
% nNodes      number nodes in the discretized domain                   |
% nx          number of volumetric source emission rate parameters     |
% oper_label  string for the file label of the precomputed finite      |
%               element operators                                      |
% oper_path   string for the file path to load a .mat file containing  |
%               finite element operators                               |
% temp        intermediate result from the computation of the adjoint  |
%               variable at a previous time step (backwards marching   |
%               in time)                                               |
% y           discrete adjoint variable (was chosen for its similarity |
%               to an upside down lower case lambda)                   |
%-----------------------------------------------------------------------

% Ensure the require number of arguments are supplied.
assert( 7 <= nargin && nargin <= 11 ); 

% Check c_k
assert( isnumeric(c_k) );               % ensure s is numeric array

% Check c_star
assert( isnumeric(c_star) );               % ensure s is numeric array


% Check t
assert( isvector(t) );            % ensure t is a vector
assert( isnumeric(t) );           % ensure t is numeric array
assert( t(1) == 0 );               % ensure t starts at 0
nt = length(t);
for n = 2 : nt                    % ensure t contains sequential time 
  assert( t(n-1) < t(n) );        %   steps
end

% Check H
assert( issparse(H) );            % ensure H is a sparse matrix
assert( isnumeric(H) );           % ensure H is a numeric array
assert( size(H, 1) == size(c_k, 2) );  % ensure dimensions of H and c_k agree
assert( size(H, 1) == size(c_star, 2) );  % ensure dimensions of H and c_star agree

% Check save_flag
assert( islogical(save_flag) );
assert( isscalar(save_flag) );

% Check working_dir
%   if save_flag is false then working_dir can contain any value.
if save_flag
  assert( ischar(working_dir) );  % ensure working_dir is a string
  assert( exist( working_dir, 'dir') == 7 );
end

% Check operator_dir
assert( ischar(oper_dir) ); % ensure operator_dir is a string
assert( exist( oper_dir, 'dir') == 7 );

% If time_offset is supplied make sure it contains valid data.
if nargin >= 8                        % if time_offset is supplied
  assert( isscalar(time_offset) );          % ensure time_offset is a scalar
  assert( isnumeric(time_offset) );         % ensure time_offset is a number
  assert( time_offset >= 0 ); % ensure time_offset is positive
else
  time_offset = 0;
end

% If nt_max is supplied make sure it contains valid data.
if nargin >= 9                         % if nt_max is supplied
  assert( isscalar(nt_max) );          % ensure nt_max is a scalar
  assert( isnumeric(nt_max) );         % ensure nt_max is a number
  assert( nt_max >= 1 ); % ensure nt_max is greater than 0
else
  nt_max = nt;
end

% If theta is supplied make sure it contains valid data.
if nargin >= 10                        % if theta is supplied
  assert( isscalar(theta) );          % ensure theta is a scalar
  assert( isnumeric(theta) );         % ensure theta is a number
  assert( 0 <= theta && theta <= 1 ); % ensure 0  <= theta <= 1
else
  theta = 0.5; % default value (Crank-Nicholson method)
end

% If c_0 is supplied make sure it contains valid data.  If c_0 is not 
% supplied then the initial condition is assumed to be the steady 
% state concentration distribution for the initial flow conditions.
if nargin == 11              % if c_0 is supplied
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


nNodes = size(H, 2);
nx = nNodes;
oper_label = 'Operators_';

if save_flag
  adj_label = 'Adjoint_';
  grad_label = 'Gradient_';
end

df_dx = zeros(nx, 1);
y = H'*(c_k(nt,:)-c_star(nt,:))';

% Save the current state of the discretized adjoint field if it is desired.
%   Make sure that the final value of the adjoint field is supposed to be the receptor residual error and not zero instead!
if save_flag
  file_num = generate_file_num(nt-1, nt);
  adjoint_path = fullfile(working_dir, [adj_label file_num '.mat']);
  %fprintf([datestr(now) ' - Saving:  ' adjoint_path '\n']);
  save(adjoint_path,'y');
end

%uPath = fullfile(oper_dir, '..', 'U.mat');
uPath = fullfile(oper_dir, 'U.mat');
load(uPath);

% load operators for the last time step
file_num = generate_file_num(nt-1+time_offset, nt_max);
oper_path = fullfile(oper_dir, [oper_label file_num '.mat']);
%fprintf([datestr(now) ' - Loading: ' oper_path '\n']);
load(oper_path,'C','K');

% copy operators so that they can still be accessed in the next time step.
Cn = C;
Kn = K;

% compute the backwards in time evolution of the adjoint variable and its contribution to the gradient of the objective function.
for n = nt-1 :-1: 1
  file_num = generate_file_num(n-1+time_offset, nt_max);
  oper_path = fullfile(oper_dir, [oper_label file_num '.mat']);
  %fprintf([datestr(now) ' - Loading: ' oper_path '\n']);
  load(oper_path,'C','K');

  dt = t(n+1) - t(n);

  temp = Cn'*(U'*((U*(Cn'/dt+theta*Kn')*U')\(U*y)));
  df_dx = df_dx + temp;
  
  % Save the current state of the objective function gradient if it is desired.
  if save_flag
    file_num = generate_file_num(n, nt);
    grad_path = fullfile(working_dir, [grad_label file_num '.mat']);
    %fprintf([datestr(now) ' - Saving:  ' grad_path '\n']);
    save(grad_path,'df_dx');
  end

  y = temp/dt -(1-theta)*K'*(C'\temp);
  y = y + H'*(c_k(n,:) - c_star(n,:))';

  % Save the current state of the discretized adjoint field if it is desired.
  if save_flag
    file_num = generate_file_num(n-1, nt);
    adjoint_path = fullfile(working_dir, [adj_label file_num '.mat']);
    %fprintf([datestr(now) ' - Saving:  ' adjoint_path '\n']);
    save(adjoint_path,'y');
  end

  % copy operators so that they can still be accessed in the next time step.
  Cn = C;
  Kn = K;

end

% if no initial conditions for the tracer concentration field are supplied then assume that
% the initial conditions of the tracer concentration field are dependent on the
% volumetric source emission rate then compute the contribution of the initial conditions 
% to the objective function gradient.
if nargin < 11
  temp = Cn'*(U'*((U*Kn'*U')\(U*y)));
  df_dx = df_dx + temp;
end

% Save the current state of the objective function gradient if it is desired.
if save_flag
  file_num = generate_file_num(0,nt);
  grad_path = fullfile(working_dir, [grad_label file_num '.mat']);
  save(grad_path,'df_dx');
end

%###################################################################################################################
%###################################################################################################################
%###################################################################################################################

%function dJ_dE = CostFncGrad2(t, nNodes)

%global simDir operDir tMax signal_o signal_c

%load(fullfile(simDir,'dirSettings.mat'));
%load(passPath);
%nt = length(t);
%theta = 0.5;


%   load(sensorPath)
%y = zeros(nNodes,1);   
%y(sensorIndex) = signal_c(nt,:)-signal_o(nt,:);
%if nPass == 0
%adjointPath = fullfile(simDir,'Adjoint', ['Adjoint_' num2str(nt-1) '.mat']);
%save(adjointPath,'y');
%end
%dJ_dE = zeros(nNodes,1);

%uPath = fullfile(simDir,'U.mat');
%load(uPath);
%operatorPath = fullfile(operDir, ['Operators_' num2str(nt-1) '.mat']);
%   load(operatorPath,'C','K');
%   Cn = C;
%   Kn = K;
%for n = nt-1 :-1: 1   
   
%   	operatorPath = fullfile(operDir, ['Operators_' num2str(n-1) '.mat']);
%   	load(operatorPath,'C','K');

%   
%   dt = t(n+1) - t(n);

%   temp = Cn'*(U'*((U*(Cn'/dt+theta*Kn')*U')\(U*y)));
   
%   dJ_dE = dJ_dE + temp;
%if nPass == 0
% gradPath = fullfile(simDir, 'Gradient', ['Gradient_' num2str(n) '.mat']);
%save(gradPath,'dJ_dE');
%end  

%   y = temp/dt -(1-theta)*K'*(C'\temp);
%   y(sensorIndex) = y(sensorIndex) + (signal_c(n,:) - signal_o(n,:))';
%   if nPass == 0
%adjointPath = fullfile(simDir,'Adjoint', ['Adjoint_' num2str(n-1) '.mat']);
%save(adjointPath,'y');
%end
%   Cn = C;
%   Kn = K;


%end

%temp = Cn'*U'*((U*Kn'*U')\(U*y));
%dJ_dE = dJ_dE + temp;
%if nPass == 0
%gradPath = fullfile(simDir, 'Gradient', 'Gradient_0.mat');
%save(gradPath,'dJ_dE');
%end
