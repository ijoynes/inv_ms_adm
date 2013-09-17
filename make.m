% Run this script to build an executable file that can be run on cluster84

fprintf('Compiling L-BFGS-B mex files...\n');

cd lbfgs

% lbfgsb 2.1
mex lbfgs_mex.c list.c utils.c routines.f

% lbfgsb 3.0
% For an unknown reason Matlab is unable to compile lbfgsb 3.0 into a 
% mex file.  For now this updated version of L-BFGS-B is not used by the
% Inverse Micro-Scale Atmospheric Dispersion Model
%
%   mex lbfgs_mex.c list.c utils.c blas.f linpack.f lbfgsb.f

movefile(['lbfgs_mex.' mexext], '..');
cd ..

fprintf('Compiling Inverse Micro-Scale Atmospheric Dispersion Model...\n');
mcc -m driver -O all -o inv_ms_adm
fprintf('Done.\n');