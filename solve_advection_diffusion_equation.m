function receptor_observations = solve_advection_diffusion_equation(t,E,c_0,type,noise)
global simDir operDir, H
load(fullfile(simDir,'dirSettings.mat'));
if type == 'c'
load(passPath);
end
theta = 0.5;
nt = length(t);





%B = ComputeBoundaryMatrix(boundary,nNodes);
operatorPath = fullfile(operDir, ['Operators_0.mat']);
load(operatorPath,'C','K');
uPath = fullfile(simDir,'U.mat');
load(uPath);
Cn = C;
Kn = K;



if type == 'o'
    % Refactored the add measurmement noise code
    % Aug 1, 2012 by Ian Joynes
    %o = c_0;
    %o = U'*((U*Kn*U')\(U*Cn*E));
    %concentrationPath = fullfile(simDir,'Observation','Observation_0.mat');
    %o_temp = o;
    %o=AddMeasurementNoise(o,noise);
    %save(concentrationPath,'o')
    %o = o_temp;

    o = c_0;
    o = U'*((U*Kn*U')\(U*Cn*E));
%    o_temp = o;
%    o_temp(o_temp<0)=0;
%    o = (o+o_temp)/2; % this is an attempt at stabilizing the solution
%    fprintf('||o(%3d)||_inf = %e\n', 0, max(abs(o)));
    concentrationPath = fullfile(simDir,'Observation','Observation_0.mat');
    noisePath = fullfile(simDir,'Noise','Noise_0.mat');

    o_error = CalcMeasurementNoise(o,noise);
    save(concentrationPath,'o')
    save(noisePath,'o_error')

elseif type == 'c'
    c = c_0;
    c = U'*((U*Kn*U')\(U*Cn*E));
%    c_temp = c;
%    c_temp(c_temp<0)=0;
%    c=(c+c_temp)/2;

 %   fprintf('||c(%d)||_inf = %e\n', 0, max(abs(c)));
    concentrationPath = fullfile(simDir,'Concentration','Concentration_0.mat');
    save(concentrationPath,'c')
	if nPass == 0
	 concentrationPath = fullfile(simDir,'Concentration_Initial_Guess','Concentration_0.mat');
    save(concentrationPath,'c')
	
	end


else
        fprintf('ERROR');
end

for n = 1 : nt - 1
%if type == 'o'
%fprintf('%s %d\n','Observation t =', n);
%elseif type == 'c'
%fprintf('%s %d\n','Concentration t =', n);    
%end
    %load(['E:\Adjoint TR 2D [2011_1_28]\Simulation_2\Operators\Operators_' num2str(5*n) '.mat']);
    operatorPath = fullfile(operDir, ['Operators_' num2str(n) '.mat']);
    load(operatorPath,'C','K');
    %U = B;
    
    dt = t(n+1) - t(n);
    A = U*(C/dt+theta*K)*U';
    %b = U*C*(c(:,n)+dt*E);
    %c(:,n+1) = U'*(A\b);
    
    
    
    
    
    if type == 'o'
        % Refactored the add measurmement noise code
        % Aug 1, 2012 by Ian Joynes
         
        %b = U*C*(o/dt-(1-theta)*(Cn\(Kn*o))+E);
        %o = U'*(A\b);
        %concentrationPath = fullfile(simDir,'Observation',['Observation_' num2str(n) '.mat']);
        %save(concentrationPath,'o')
        %o_temp = o;
        %o=AddMeasurementNoise(o,noise);
        %save(concentrationPath,'o')
        %o = o_temp;

        b = U*C*(o/dt-(1-theta)*(Cn\(Kn*o))+E);
        o = U'*(A\b);
        %o(o<0) = 0;
%	o_temp = o;
%	o_temp(o_temp<0)=0;
%	o=(o+o_temp)/2;
  %      fprintf('||o(%3d)||_inf = %e\n', n, max(abs(o)));
        concentrationPath = fullfile(simDir,'Observation',['Observation_' num2str(n) '.mat']);
        noisePath = fullfile(simDir,'Noise', ['Noise_' num2str(n) '.mat']);
        o_error = CalcMeasurementNoise(o,noise);
        save(concentrationPath,'o');
        save(noisePath, 'o_error');

    elseif type == 'c'
        b = U*C*(c/dt-(1-theta)*(Cn\(Kn*c))+E);
        c = U'*(A\b);
        %c(c<0) = 0;
%	c_temp = c;
%	c_temp(c_temp<0)=0;
%	c = (c + c_temp)/2;
   %     fprintf('||c(%3d)||_inf = %e\n', n, max(abs(c)));
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
%fprintf('END OF SolveConcentrationTransport.m\n');

