function [simDir, operDir, domainPath, paramPath, sensorPath, ... 
  sourcePath, passPath, tMax, dt, noise, maxIter, reg_par, factr, ...
  pgtol, m, iprint, save_flag, time_offset, nt_max, est_min] = readControls(filePath)
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
%   tMax      The final time of the simulation in seconds.  The 
%               default value is tMax = 900.
%
%   dt        The time step size in seconds.  The default value is 
%               dt = 1.
%
%   noise     The slope of the line which relates the observation 
%               concentration to the standard deviation of the 
%               Gaussian distributed noise (e.g. 0.1 = 10% noise).
%               A concentration measurement of 0 mg/m^3 or less is
%               associated will a minimum Gaussian noise standard 
%               deviation of ~19.24 mg/m^3.
%
%   maxIter   The maximum number of iterations attempted by the
%               Inverse Micro-Scale Dispersion Model.  The default
%               value is 200 iterations.
%
%   reg_par   The regularization parameter which determines the 
%               importance of having a sparse source distribution vs. 
%               source distribution that explains the receptor 
%               observations.  The default value is 0, which adds no
%               penalty to complex source distributions.
%
%   factr     The exit criteria for small changes to the objective 
%               function.  The default value is factr = 0, which 
%               practically disables this minimization exit test.  
%               More details on this parameter is given in section 3
%               of:
% 
%               Zhu, Ciyou; Byrd, Richard H.; Lu, Peihuang and 
%               Nocedal, Jorge., (1997) "Algorithm 778: L-BFGS-B: 
%               Fortran subroutines for large-scale bound-constrained 
%               optimization", ACM Transactions on Mathematical 
%               Software 23(4), 550-560.
%               
% 
%   pgtol     The exit criteria for small changes to the gradient.  
%               The default value is pgtol = 0, which practically
%               disables this minimization exit test.  More details on 
%               this parameter is given in section 3 of:
% 
%               Zhu, Ciyou; Byrd, Richard H.; Lu, Peihuang and 
%               Nocedal, Jorge., (1997) "Algorithm 778: L-BFGS-B: 
%               Fortran subroutines for large-scale bound-constrained 
%               optimization", ACM Transactions on Mathematical 
%               Software 23(4), 550-560.
%
%   m         The number of histories to use to approximate the 
%               Hessian matrix.  The default value is m = 5.
%
%  iprint     An integer variable that controls the frequency and type
%               of output generated by the L-BFGS-B routine:
%               
%               iprint<0    no output is generated;
%               iprint=0    print only one line at the last iteration;
%               0<iprint<99 print also f and |proj g| every iprint 
%                             iterations;
%               iprint=99   print details of every iteration except 
%                             n-vectors;
%               iprint=100  print also the changes of active set and 
%                             final x;
%               iprint>100  print details of every iteration including
%                             x and g;
%               
%               When iprint > 0, the file iterate.dat will be created 
%               to summarize the iteration.
%               
%               When iprint is not specified the default value of
%               iprint = -1 is used to suppress all output
%
%  save_flag
%
%  time_offset  A positive integer variable that controls the number of 
%               time-steps the flow conditions are offset from the 
%               observations.  The default value is time_offset = 0
%
%-----------------------------------------------------------------------


% set default simulation control parameters
simDir  = '';
operDir = '';
domainName = 'Domain.mat';
sensorName = 'Sensors.mat';
sourceName = 'Source.mat';
tMax    = 900;
dt      = 1;
noise   = 0;
maxIter = 200;  % The default for lbfgs_options.m is actually 100
reg_par = 0;
pgtol   = 0;      % The default from lbfgs_options.m is actually 1E-5
factr   = 0;      % The default from lbfgs_options.m is actually 1E7
m       = 5;      % default from lbfgs_options.m
iprint  = -1;     % The default from lbfgs_options.m is actually 0
save_flag = false;
time_offset = 0;
nt_max = 901;
est_min = false;

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
        case 'domainName'
            domainName = fscanf(fid, '%s', 1);
        case 'sensorName'
            sensorName = fscanf(fid, '%s', 1);
        case 'sourceName'
            sourceName = fscanf(fid, '%s', 1);
        case 'tMax'
            tMax = fscanf(fid, '%d',1);
            assert(tMax>0);
        case 'dt'
            dt = fscanf(fid, '%d',1);
            assert(dt>0);
        case 'noise'
            noise = fscanf(fid, '%e', 1);
            assert(noise>=0 || noise == -1);
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
        case 'save_flag'
            temp = fscanf(fid, '%s', 1);
            if strcmpi(temp, 'true')
                save_flag = true;
            else
                save_flag = false;
            end
        case 'time_offset'
            time_offset = fscanf(fid,'%d',1);
        case 'nt_max'
            nt_max = fscanf(fid, '%d',1);
        case 'est_min'
            temp = fscanf(fid, '%s', 1);
            if strcmpi(temp, 'true')
                est_min = true;
            else
                est_min = false;
            end
    end

    if ~feof(fid)
        garbage = fgetl(fid);
    end
end
fclose(fid);

domainPath = fullfile(simDir, domainName);
paramPath  = fullfile(simDir,'Parameters.mat');
sensorPath = fullfile(simDir, sensorName);
sourcePath = fullfile(simDir, sourceName);
passPath   = fullfile(simDir,'pass.mat');

% ensure that the parameters contain valid values

assert( exist( simDir, 'dir') == 7  );
assert( exist(operDir, 'dir') == 7  );
assert(   noise >= 0 || noise == -1 );
assert( maxIter >  0 );
assert(    tMax >  0 );
assert(      dt >  0 );
assert( reg_par >= 0 );
assert(   factr >= 0 );
assert(   pgtol >= 0 );
assert(       m >  0 );
assert( time_offset >= 0);
assert( nt_max >= 1 );