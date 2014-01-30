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



domainDir = run_dir{1};
load( fullfile(domainDir, 'Domain.mat') );
xy =xy - ones(nNodes,1)*min(xy);
s2m = zeros(nNodes,1); 
for i = 1 : nTris
    s2m(tri(i,:)) = s2m(tri(i,:)) + det([ones(3,1) xy(tri(i,:),:)]);
end
s2m = s2m/6;  

iter = [25,50,100,200]';
for j = 1:4
    
load(fullfile(run_dir{6},'Source',['Source_' num2str(iter(j)) '.mat']))
s_max = max(s);
s(s<0.5*s_max) = 0;

area = 0;
num_elem = 0;
for i = 1 : nTris
    if any(s(tri(i,:)) > 0)
        area = area + det([ones(3,1), xy(tri(i,:),:)])/2;
        num_elem = num_elem + 1;
    end
end
fprintf('%.0f (%d)\n', area, num_elem)
end