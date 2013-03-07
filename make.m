% Run this script to build an executable file that can be run on cluster84

fprintf('Compiling lbfgs mex files...\n');
fprintf('Change directory to lbfgs for this to work.\n');
cd lbfgs\
mex lbfgs_mex.c list.c utils.c routines.f

fprintf('Done.\n');

mcc -m driver.m -o imsadm.out