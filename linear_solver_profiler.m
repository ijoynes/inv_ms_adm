sim_dir = '/home/ijoynes/run_dummy/';


oper_dir = '/files/j_group/ijoynes/balzac_2d_simulations/transient/FlatExportTrans775/Operators/';

operator_path = fullfile(oper_dir, 'Operators_0.mat');

load( operator_path, 'C', 'K' );
Cn = C;
Kn = K;

operator_path = fullfile(oper_dir, 'Operators_1.mat');
load(operator_path,'C','K');
load(fullfile(sim_dir, 'U.mat'));
dt = 1;
theta = 0.5;
n = size(C,1);
s = zeros(n,1);
s(10000) = 1;
A = U*(C/dt+theta*K)*U';
t = nan(10,1);
r = nan(10,1);

% Backslash
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = A\b;
c = U'*x;
toc;
t(1) = toc;
r(1) = norm(A*x-b)/norm(b);

% lsqr
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = lsqr(A,b,1E-15,100);
c = U'*x;
toc;
t(2) = toc;
r(2) = norm(A*x-b)/norm(b);

% cgs
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = cgs(A,b,1E-15,100);
c = U'*x;
toc;
t(3) = toc;
r(3) = norm(A*x-b)/norm(b);

% cgs
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = gmres(A,b,[],1E-15,100);
c = U'*x;
toc;
t(4) = toc;
r(4) = norm(A*x-b)/norm(b);

% bicg
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = bicg(A,b,1E-15,100);
c = U'*x;
toc;
t(5) = toc;
r(5) = norm(A*x-b)/norm(b);

% bicgstab
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = bicgstab(A,b,1E-15,100);
c = U'*x;
toc;
t(6) = toc;
r(6) = norm(A*x-b)/norm(b);

% bicgstabl
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = bicgstabl(A,b,1E-15,100);
c = U'*x;
toc;
t(7) = toc;
r(7) = norm(A*x-b)/norm(b);

% qmr
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = qmr(A,b,1E-15,100);
c = U'*x;
toc;
t(8) = toc;
r(8) = norm(A*x-b)/norm(b);

% pcg
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = pcg(A,b,1E-15,100);
c = U'*x;
toc;
t(9) = toc;
r(9) = norm(A*x-b)/norm(b);

% minres
c = zeros(n,1);
tic;
temp = Cn\(Kn*c);
b = U*(C*(c/dt-(1-theta)*(temp+s)));
x = minres(A,b,1E-15,100);
c = U'*x;
toc;
t(10) = toc;
r(10) = norm(A*x-b)/norm(b);

fprintf('')


% 1    pcg         - Preconditioned Conjugate Gradients Method.
% 2    bicg        - BiConjugate Gradients Method.
% 3    bicgstab    - BiConjugate Gradients Stabilized Method.
% 4    bicgstabl   - BiCGStab(l) Method.
% 5    cgs         - Conjugate Gradients Squared Method.
% 6    gmres       - Generalized Minimum Residual Method.
% 7    lsqr        - LSQR Method.
% 8    minres      - Minimum Residual Method.
% 9    qmr         - Quasi-Minimal Residual Method.
% 10   symmlq      - Symmetric LQ Method.
% 11   tfqmr       - Transpose-Free QMR Method.


t
r