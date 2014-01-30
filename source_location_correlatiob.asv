clc
close
clear

nRuns = 10;

data_dir = fullfile('E:', 'ijoynes', 'thesis_data_backup');

run_dir = {fullfile(data_dir, 'run_031');
   fullfile(data_dir, 'run_068');
   fullfile(data_dir, 'run_069');
   fullfile(data_dir, 'run_070');
   fullfile(data_dir, 'run_071');
   fullfile(data_dir, 'run_046');
   fullfile(data_dir, 'run_055');
   fullfile(data_dir, 'run_047');
   fullfile(data_dir, 'run_045');
   fullfile(data_dir, 'run_048') };

% run_dir = {fullfile(data_dir, 'run_058');
%     fullfile(data_dir, 'run_059');
%     fullfile(data_dir, 'run_060');
%     fullfile(data_dir, 'run_061');
%     fullfile(data_dir, 'run_062');
%     fullfile(data_dir, 'run_063');
%     fullfile(data_dir, 'run_064');
%     fullfile(data_dir, 'run_065');
%     fullfile(data_dir, 'run_066');
%     fullfile(data_dir, 'run_067'); };

run_name = {'No noise, No reg';
    'No noise, theta = 1';
    'No noise, theta = 10';
    'No noise, theta = 100';
    'No noise, theta = 1000';
    '10% noise, No reg';
    '10% noise, theta = 1'
    '10% noise, theta = 10';
    '10% noise, theta = 100';
    '10% noise, theta = 1000'};

theta = [0 1 10 100 1000 0 1 10 100 1000];

domainDir = run_dir{1};
load( fullfile(domainDir, 'Domain.mat') );

s2m = zeros(nNodes,1); 
for i = 1 : nTris
    s2m(tri(i,:)) = s2m(tri(i,:)) + det([ones(3,1) xy(tri(i,:),:)]);
end
s2m = s2m/6;  

iRun = 1;
nIters = zeros(nRuns,1);
for i = 1 : nRuns
    while 1
        srcPath = fullfile(run_dir{i}, 'Source', ['Source_' int2str(nIters(i)) '.mat']);
        if exist(srcPath )~= 2
            break;
        end
        nIters(i) = nIters(i) + 1;
    end
end

load(fullfile(run_dir{1}, 'Source', 'Source_Correct.mat'));
src_idx = find(E>0);

wgt = zeros(nNodes,1);
for i = 1 : length(src_idx)
    wgt = wgt + 1./sum((xy-ones(nNodes,1)*xy(src_idx(i),:)).^2,2);
end
wgt(src_idx) = 1E6;
wgt=wgt/length(src_idx);

B_inv = compute_inverse_source_error_covariance_matrix(tri, xy);

%E_norm = E/max(E);
E_norm = E;
r_s = nan(max(nIters),nRuns);

for iRun = 1 : nRuns
iter = 0 : (nIters(iRun)-1);
for i = 1 : nIters(iRun)
    i
    load(fullfile(run_dir{iRun}, 'Source', ['Source_' int2str(iter(i)) '.mat']));
    %s_norm = s/max(s);
    s_norm = s;
    %r_s(i,iRun) = (wgt'*B_inv*s_norm)/(wgt'*B_inv*E_norm);
    r_s(i,iRun) = (E'*B_inv*s)/(E'*B_inv*E);
end
end

figure(1)
plot(iter,r_s(:,1:5))
axis([0,200,0,0.4])
grid
legend(run_name(1:5),'Location','NorthWest')
title('Convergence History of Single Source Location Correlation (No Noise)')
xlabel('Iteration #')
ylabel('Emission Source Location Correlation')

figure(2)
plot(iter,r_s(:,6:10))
axis([0,200,0,0.4])
grid
legend(run_name(6:10),'Location','NorthWest')
title('Convergence History of Single Source Location Correlation (10% Gaussian Noise)')
xlabel('Iteration #')
ylabel('Emission Source Location Correlation')

xy = xy - ones(nNodes,1)*min(xy);
trisurf(tri,xy(:,1),xy(:,2),log(wgt),'edgecolor','interp','facecolor','interp')
view(2)
axis image

