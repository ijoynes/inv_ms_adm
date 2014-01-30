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

load(fullfile(run_dir{1}, 'Source', 'Source_Correct.mat'))

[E_max, E_max_i] = max(E);

loc_error = nan(201,1);

for i = 0 : 200
    i
    load(fullfile(run_dir{6}, 'Source', ['Source_' num2str(i) '.mat']))
    [s_max, s_max_i] = max(s);
    loc_error(i+1) =  sqrt(sum((xy(s_max_i,:)-xy(E_max_i,:)).^2,2));
end
find(loc_error == 0, 1, 'first')-1
plot(0:200,loc_error)


loc_error([26,51,101,201])