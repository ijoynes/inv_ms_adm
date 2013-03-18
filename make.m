% Run this script to build an executable file that can be run on cluster84

fprintf('Compiling L-BFGS-B mex files...\n');

cd lbfgs
mex lbfgs_mex.c list.c utils.c routines.f
movefile(['lbfgs_mex.' mexext], '..');
cd ..

% fprintf('Compiling Inverse Micro-Scale Atmospheric Dispersion Model...\n');
% mcc -m driver.m -o imsadm.out
fprintf('Done.\n');