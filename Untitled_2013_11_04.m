clc
close
clear

nRuns = 10;

data_dir = fullfile('E:', 'ijoynes', 'thesis_data_backup');

% run_dir = {fullfile(data_dir, 'run_031');
%    fullfile(data_dir, 'run_068');
%    fullfile(data_dir, 'run_069');
%    fullfile(data_dir, 'run_070');
%    fullfile(data_dir, 'run_071');
%    fullfile(data_dir, 'run_046');
%    fullfile(data_dir, 'run_055');
%    fullfile(data_dir, 'run_047');
%    fullfile(data_dir, 'run_045');
%    fullfile(data_dir, 'run_048') };

run_dir = {fullfile(data_dir, 'run_058');
    fullfile(data_dir, 'run_059');
    fullfile(data_dir, 'run_060');
    fullfile(data_dir, 'run_061');
    fullfile(data_dir, 'run_062');
    fullfile(data_dir, 'run_063');
    fullfile(data_dir, 'run_064');
    fullfile(data_dir, 'run_065');
    fullfile(data_dir, 'run_066');
    fullfile(data_dir, 'run_067'); };

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
xy =xy - ones(nNodes,1)*min(xy);
s2m = zeros(nNodes,1); 
for i = 1 : nTris
    s2m(tri(i,:)) = s2m(tri(i,:)) + det([ones(3,1) xy(tri(i,:),:)]);
end
s2m = s2m/6;  

load(fullfile(run_dir{10},'Source','Source_200.mat'))

[max_s_v, max_s_i]=max(s);
s(s<0.01*max_s_v) = 1E-100;
trisurf(tri,xy(:,1),xy(:,2),log(s),'edgecolor','interp','facecolor','interp')
view(2)
axis image
colorbar
