function df_dx = solve_adjoint_advection_diffusion_equation(c_k, t, H, working_dir, oper_dir, theta, c_0)









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
