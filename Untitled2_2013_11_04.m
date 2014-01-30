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


load(fullfile(run_dir{6}, 'Source','Source_Correct.mat'),'signal_o')
c_o = signal_o;
c_o_mean = mean(mean(c_o));

load(fullfile(run_dir{1}, 'Source','Source_Correct.mat'),'signal_o')
c = signal_o;

r_max = 1-sum(sum((c-c_o).^2))/sum(sum((c_o_mean-c_o).^2))
