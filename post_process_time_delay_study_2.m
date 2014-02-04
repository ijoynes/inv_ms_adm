% VARIABLES
%  1) x [m]             <times><i>x</i> [m]</times>
%  2) y [m]             <times><i>y</i> [m]</times>
%  3) s [mg*s^-1*m^-3]  <times><i>s</i>(<b>x</b><sub><i>k</i></sub>) [mg</times><math>W</math><times>s<sup>-1</sup></times><math>W</math><times>m<sup>-3</sup>]</times>
%  4) grad(f)           <math>Q&</math><times>(<b>x</b><sub><i>k</i></sub>) [g</times><math>W</math><times>s</times><math>W</math><times>m<sup>-3</sup>]</times>
%     grad(f)           <math>Q&</math><times>(<b>x</b><sub><i>k</i></sub>) [g</times><math>W</math><times>s</times><math>W</math><times>m<sup>-3</sup>]</times>
%  5) f                 <math>&</math><times>(<b>x</b><sub><i>k</i></sub>) [kg<sup>2</sup></times><math>W</math><times>m<sup>-3</sup>]</times>
%  6) ||grad(f)||       <times>||</times><math>Q&</math><times>(<b>x</b><sub><i>k</i></sub>)||</times><sub><math>%</math></sub><times> [kg</times><math>W</math><times>s</times><math>W</math><times>m<sup>-3</sup>]</times>
%  7) ||proj-grad(f)||              <times>||</times><math>Q</math><times><sub>proj</sub></times><math>&</math><times>(<b>x</b><sub><i>k</i></sub>)||</times><sub><math>%</math></sub><times> [kg</times><math>W</math><times>s</times><math>W</math><times>m<sup>-3</sup>]</times>
%  8) Q                             <times><i>Q</i>(<b>x</b><sub><i>k</i></sub>) [g</times><math>W</math><times>s<sup>-1</sup>]</times>
%  9) f/f_0                         <math>&</math><times>(<b>x</b><sub><i>k</i></sub>) /</times><math>&</math><times>(<b>x</b><sub><i>k</i></sub>)</times>
% 10) ||grad(f)||/||grad(f_0)||     <times>||</times><math>Q&</math><times>(<b>x</b><sub><i>k</i></sub>)||</times><sub><math>%</math></sub><times> / ||</times><math>Q&</math><times>(<b>x</b><sub>0</sub>)||</times><sub><math>%</math></sub>
% 11) <times>||</times><math>Q</math><times><sub>proj</sub></times><math>&</math><times>(<b>x</b><sub><i>k</i></sub>)||</times><sub><math>%</math></sub><times> / ||</times><math>Q</math><times><sub>proj</sub></times><math>&</math><times>(<b>x</b><sub>0</sub>)||</times><sub><math>%</math></sub>
% 12) <times><i>Q</i>(<b>x</b><sub><i>k</i></sub>) / <i>Q</i>*</times>
% 13) <times><i>r</i><sup>2</sup>(<b>x</b><sub><i>k</i></sub>)</times>
% 14) <times><i>t</i> [min]</times>
% 15) <times><i>c</i>(<b>x</b><sub><i>k</i></sub>) [mg</times><math>W</math><times>m<sup>-3</sup>]</times>

close
clear
clc

dir_name = {'E:\ijoynes\thesis_data_backup\run_072', ...
            'E:\ijoynes\thesis_data_backup\run_073', ...
            'E:\ijoynes\thesis_data_backup\run_074', ...
            'E:\ijoynes\thesis_data_backup\run_075', ...
            'E:\ijoynes\thesis_data_backup\run_076', ...
            'E:\ijoynes\thesis_data_backup\run_077', ...
            'E:\ijoynes\thesis_data_backup\run_078', ...
            'E:\ijoynes\thesis_data_backup\run_079'};


nRuns = length(dir_name);
time_delay = [0, 1, 2, 5, 10, 20, 30, 60]';
nIters = [36, 39, 38, 39, 23, 23, 42, 42]';
f_all = nan(max(nIters)+1,nRuns);
g_norm_all = nan(max(nIters)+1,nRuns);
tecplot_dat_file_path = 'E:\ijoynes\thesis_data_backup\time_delay_study.dat';
tecplot_plt_file_path = 'E:\ijoynes\thesis_data_backup\time_delay_study.plt';

fid = fopen(tecplot_dat_file_path,'w');
fprintf(fid, 'TITLE = "Time Delay Study"\n');
fprintf(fid, 'FILETYPE = FULL\n');
fprintf(fid, 'VARIABLES = ');
fprintf(fid, '"<times><i>x</i> [m]</times>",\n'); % x
fprintf(fid, '"<times><i>y</i> [m]</times>",\n'); % y
fprintf(fid, '"<times><i>s</i> [mg</times><math>W</math><times>s<sup>-1</sup></times><math>W</math><times>m<sup>-3</sup>]</times>",\n'); % s
fprintf(fid, '"<math>Q&</math><times>(<b>x</b><sub><i>k</i></sub>) [g</times><math>W</math><times>s</times><math>W</math><times>m<sup>-3</sup>]</times>",\n'); % grad(f)
fprintf(fid, '"<math>Q</math><times><sub>proj</sub></times><math>&</math><times>(<b>x</b><sub><i>k</i></sub>) [g</times><math>W</math><times>s</times><math>W</math><times>m<sup>-3</sup>]</times>",\n'); % grad(f)
fprintf(fid, '"<times><i>k</i></times>" \n');
fprintf(fid, '"<math>&</math><times>(<b>x</b><sub><i>k</i></sub>) [kg<sup>2</sup></times><math>W</math><times>m<sup>-6</sup>]</times>",\n'); % f
fprintf(fid, '"<times>||</times><math>Q&</math><times>(<b>x</b><sub><i>k</i></sub>)||</times><sub><math>%%</sub></math><times> [g</times><math>W</math><times>s</times><math>W</math><times>m<sup>-3</sup>]</times>",\n'); % ||grad(f)||
%fprintf(fid, '"<times>||</times><math>Q</math><times><sub>proj</sub></times><math>&</math><times>(<b>x</b><sub><i>k</i></sub>)||</times><sub><math>%%</sub></math><times> [g</times><math>W</math><times>s</times><math>W</math><times>m<sup>-3</sup>]</times>",\n'); % ||proj-grad(f)||
%fprintf(fid, '"<times><i>Q</i>(<b>x</b><sub><i>k</i></sub>) [g</times><math>W</math><times>s<sup>-1</sup>]</times>",\n'); % Q
%fprintf(fid, '"<math>&</math><times>(<b>x</b><sub><i>k</i></sub>) /</times><math>&</math><times>(<b>x</b><sub><i>k</i></sub>)</times>",\n'); % f/f_0
%fprintf(fid, '"<times>||</times><math>Q&</math><times>(<b>x</b><sub><i>k</i></sub>)||</times><sub><math>%%</math></sub><times> / ||</times><math>Q&</math><times>(<b>x</b><sub>0</sub>)||</times><sub><math>%%</math></sub>",\n'); % ||grad(f)||/||grad(f_0)||
%fprintf(fid, '"<times>||</times><math>Q</math><times><sub>proj</sub></times><math>&</math><times>(<b>x</b><sub><i>k</i></sub>)||</times><sub><math>%%</math></sub><times> / ||</times><math>Q</math><times><sub>proj</sub></times><math>&</math><times>(<b>x</b><sub>0</sub>)||</times><sub><math>%%</math></sub>",\n'); % ||proj-grad(f)||/||proj-grad(f_0)||
%fprintf(fid, '"<times><i>Q</i>(<b>x</b><sub><i>k</i></sub>) / <i>Q</i>*</times>",\n');
%fprintf(fid, '"<times><i>r</i><sup>2</sup>(<b>x</b><sub><i>k</i></sub>)</times>",\n');
%fprintf(fid, '"<times><i>t</i> [min]</times>",\n');
%fprintf(fid, '"<times><i>c</i>(<b>x</b><sub><i>k</i></sub>) [mg</times><math>W</math><times>m<sup>-3</sup>]</times>"\n');

% Known Source
fprintf('Writing Known Source ...\n');
load(fullfile(dir_name{1}, 'Domain.mat'), 'xy', 'tri', 'nNodes', 'nTris')
load(fullfile(dir_name{1}, 'Source', 'Source_Correct.mat'), 's_star')
xy = xy - ones(nNodes,1)*min(xy);

fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Known Source"\n');
fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
fprintf(fid, 'NODES = %d\n', nNodes);
fprintf(fid, 'ELEMENTS = %d\n', nTris);
%fprintf(fid, 'PASSIVEVARLIST = [4-16]\n');
fprintf(fid, 'PASSIVEVARLIST = [4-8]\n');
fprintf(fid, '%e\n', xy(:,1));
fprintf(fid, '%e\n', xy(:,2));
fprintf(fid, '%e\n', s_star*1E6);
fprintf(fid, '%d %d %d\n', tri');

for i = 1:nRuns
    fprintf('Writing Run #%d ...\n',i+71);
    for j = 0:nIters(i)
        file_num = generate_file_num(j, 200);
        load(fullfile(dir_name{i}, 'Source', ['Source_' file_num '.mat']),'s','g','g_proj','f')
        f_all(j+1,i) = f;
        g_norm_all(j+1,i) = norm(g,Inf);
fprintf(fid, '\nZONE\n');

fprintf(fid, 'T = "TD = %d, It = %d"\n', time_delay(i), j);
fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
fprintf(fid, 'NODES = %d\n', nNodes);
fprintf(fid, 'ELEMENTS = %d\n', nTris);
fprintf(fid, 'VARSHARELIST = ([1,2]=1)\n');
%fprintf(fid, 'PASSIVEVARLIST = [6-16]\n');
fprintf(fid, 'PASSIVEVARLIST = [6-8]\n');
fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
fprintf(fid, 'STRANDID = %d\n', i);
fprintf(fid, 'SOLUTIONTIME = %d\n', j);
fprintf(fid, '%e\n', s*1E6);
fprintf(fid, '%e\n', g*1E3);
fprintf(fid, '%e\n', g_proj*1E3);
%fprintf(fid, 'AUXDATA time_delay = "%d"\n', time_delay(i));
    end
end

for i = 1 : nRuns
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Convergence: TD = %d"\n', i);
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', nIters(i)+1);
fprintf(fid, 'PASSIVEVARLIST = [1-5]\n');
fprintf(fid, '%d\n', 0:nIters(i));
fprintf(fid, '%e\n', f_all(1:nIters(i)+1,i));
fprintf(fid, '%e\n', g_norm_all(1:nIters(i)+1,i));
end

fclose(fid);
system(['preplot ' tecplot_dat_file_path ' ' tecplot_plt_file_path]);
