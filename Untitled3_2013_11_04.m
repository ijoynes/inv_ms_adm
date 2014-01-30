clc
close
clear


cdf_threshold = [.999 .99 .95 .9];

beta = [0,1,10,100,1000,0,1,10,100,1000]';

% run_dir = {fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_058');
%     fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_059');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_060');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_061');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_062');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_063');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_064');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_065');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_066');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_067'); };

run_dir = {fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_031');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_068');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_069');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_070');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_071');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_046');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_055');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_047');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_045');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_048') };


beta = [0, 1, 10, 100, 1000, 0, 1, 10, 100, 1000]';

domainDir = run_dir{1};
load( fullfile(domainDir, 'Domain.mat') );
xy =xy - ones(nNodes,1)*min(xy);
s2m = zeros(nNodes,1); 
for i = 1 : nTris
    s2m(tri(i,:)) = s2m(tri(i,:)) + det([ones(3,1) xy(tri(i,:),:)]);
end
s2m = s2m/6;  

f_min = nan(10,1);
f_v = nan(201,10);
for i = 1 : 10
i
load(fullfile(run_dir{i}, 'Source','Source_Correct.mat'),'signal_o')
c_o = signal_o;
c_o_mean = mean(mean(c_o));

load(fullfile(run_dir{1}, 'Source','Source_Correct.mat'),'signal_o','E')
c = signal_o;

f_min(i) = sum(sum((c-c_o).^2))/2+beta(i)*dot(s2m,E.^2)/2;

load(fullfile(run_dir{1},'Source','Source_0.mat'))
f0 = f;



%for j = 0 : 200
%    load(fullfile(run_dir{i},'Source',['Source_' num2str(j) '.mat']))

%   f_v(j+1,i) = f;
%end



end
figure(1)
semilogy(0:200, (f_v(:,1:5)-ones(201,1)*f_min(1:5)')/f0)


figure(2)
semilogy(0:200, (f_v(:,6:10)-ones(201,1)*f_min(6:10)')/f_v(1,6))

load(fullfile(run_dir{10},'Source','Source_Correct.mat'))
fprintf('f/f0 = %e\n', f/f0);
fprintf('f_obs/f0 = %e\n', (f-beta(i)/2*dot(s2m,s.^2))/f0);