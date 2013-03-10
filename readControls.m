function [simDir, operDir, tMax, noise, maxIter, reg_par, factr, ... 
    pgtol, m, iprint] = readControls(filePath)
%readControls reads an ASCII text file containing the values of the 
% parameters that control the execution of the Inverse Mirco-Scale
% Atmospheric Dispersion Model.
%
% INPUTS:
%
%   filePath  character string containing the complete file path to
%               the text file containing the values of the control
%               parameters.
%
%
% OUTPUTS:
%
%   simDir    character string containing the complete path to the  
%               directory containing the simulation inputs and 
%               results.  If no value is supplied by the control 
%               file, then the default value of an empty string is 
%               used.
%
%   operDir   character string containing the complete path to the  
%               directory containing the simulation operators
%               (sparse finite element capacitance and stiffness 
%               matrices).  If no value is supplied by the control  
%               file, then the default value of an empty string is
%               used.
%
%   noise     The slope of the line which relates the observation 
%               concentration to the standard deviation of the 
%               Gaussian distributed noise (e.g. 0.1 = 10% noise).
%               A concentration measurement of 0 mg/m^3 or less is
%               associated will a minimum Gaussian noise standard 
%               deviation of ~19.24 mg/m^3.
%

% set default simulation control parameters
simDir = '';
operDir = '';
noise = 0;
maxIter = 200;  % The default for lbfgs_options.m is actually 100
reg_par = 0;
pgtol = 1E-5;   % default from lbfgs_options.m
factr = 1E7;    % default from lbfgs_options.m
m = 5;      % default from lbfgs_options.m
iprint = 0;     % default from lbfgs_options.m

% open file control file and read the simulation control parameters
fid = fopen(filePath,'r');
assert(fid>0);
while ~feof(fid)
    varName = fscanf(fid,'%s',1);

    switch varName
        case 'simDir'
            simDir = fscanf(fid,'%s',1);
        case 'operDir'
            operDir = fscanf(fid, '%s', 1);
        case 'tMax'
            tMax = fscanf(fid, '%d',1);
            assert(tMax>0);
        case 'noise'
            noise = fscanf(fid, '%e', 1);
            %assert(noise>=0);
        case 'maxIter'
            maxIter = fscanf(fid,'%d',1);
            assert( maxIter > 0 );
        case 'reg_par'
            reg_par = fscanf(fid,'%e',1);
            assert(reg_par>=0);
        case 'factr'
            factr = fscanf(fid,'%e',1);
            assert(factr>=0);
        case 'pgtol'
            pgtol = fscanf(fid,'%e',1);
            assert(pgtol>=0);
        case 'm'
            m = fscanf(fid,'%d',1);
            assert(m>0);
        case 'iprint'
            iprint = fscanf(fid,'%d',1);
    end

    if ~feof(fid)
        garbage = fgetl(fid);
    end
end
fclose(fid);

assert( exist( simDir, 'dir') == 7 );
assert( exist(operDir, 'dir') == 7 );
%assert(   noise >= 0 );
assert( maxIter >  0 );
assert(    tMax >  0 );
assert( reg_par >= 0 );
assert(   factr >= 0 );
assert(   pgtol >= 0 );
assert(       m >  0 );
