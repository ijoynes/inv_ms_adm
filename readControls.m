function [simDir, operDir, tMax, noise, maxIter, reg_par, factr, pgtol, m, iprint] = readControls(filePath)

% set default simulation control parmeters
simDir = '';
operDir = '';
noise = 0;
maxIter = 200;	% The default for lbfgs_options.m is actually 100
reg_par = 0;
pgtol = 1E-5; 	% default from lbfgs_options.m
factr = 1E7;	% default from lbfgs_options.m
m = 5;		% default from lbfgs_options.m
iprint = 0; 	% default from lbfgs_options.m

% open file control file and read the simulation control parmeters
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
