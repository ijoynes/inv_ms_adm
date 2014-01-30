close
clear
clc

src_dir = 'E:\ijoynes\thesis_data_backup\run_023';
%flow_src_dir = 'E:\ijoynes\thesis_data_backup\run_023';
dst_dir = 'E:\ijoynes\thesis_data_backup';
dst_name = 'Single_Source_Projected_Gradient_2013-12-11';

C_mu = 0.085;
Sc_t = 0.9;
nt = 900;

load(fullfile(src_dir,'Domain.mat'),'xy','tri','nNodes','nTris');
load(fullfile(src_dir,'Source.mat'));
load(fullfile(src_dir,'Sensors.mat'));



nSources = size(source_xy,1);
nSensors = size(receptor_xy,1);

source_xy = source_xy - ones(nSources,1)*min(xy);
receptor_xy = receptor_xy - ones(nSensors,1)*min(xy);

xy = xy - ones(nNodes,1)*min(xy);
sourceIndex = placeSensors(xy,source_xy);

fid = fopen([fullfile(dst_dir,dst_name) '.dat'],'w');


fprintf(fid, 'TITLE = "Single Source Gradient"\n');
fprintf(fid, [ 'VARIABLES = "x [m]" ', ...  %1
                           '"y [m]" ', ...  %2
                           '"<math>Q&</math> [g<math>W</math>s<math>W</math>m<sup>-3</sup>]"\n', ...
                           '"<math>Q</math><sub>proj</sub><math>&</math> [g<math>W</math>s<math>W</math>m<sup>-3</sup>]"\n'] ); %3

n = 0;
fprintf(fid, '\nZONE\n');
load(fullfile(src_dir,'Source', ['Source_' num2str(n) '.mat']),'g','s');
fprintf(fid, 'T = "Gradient: %d"\n', n);
fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
fprintf(fid, 'NODES = %d\n', nNodes);
fprintf(fid, 'ELEMENTS = %d\n', nTris);
fprintf(fid, 'STRANDID = 1\n');
fprintf(fid, 'SOLUTIONTIME = %d\n', n);
fprintf(fid, '%e\n', xy(:,1));
fprintf(fid, '%e\n', xy(:,2));
fprintf(fid, '%e\n', g*1E3);
idx = g>0;
g(idx) = min(s(idx),g(idx));
fprintf(fid, '%e\n', g*1E3);
fprintf(fid, '%d %d %d\n', tri');

for i = 1 : 200

fprintf(fid, '\nZONE\n');
load(fullfile(src_dir,'Source', ['Source_' num2str(n) '.mat']),'g','s');
fprintf(fid, 'T = "Gradient: %d"\n', n);
fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
fprintf(fid, 'NODES = %d\n', nNodes);
fprintf(fid, 'ELEMENTS = %d\n', nTris);
fprintf(fid, 'STRANDID = 1\n');
fprintf(fid, 'SOLUTIONTIME = %d\n', n);
fprintf(fid, 'VARSHARELIST = ([1,2]=1)\n');
fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
fprintf(fid, '%e\n', g*1E3);
idx = g>0;
g(idx) = min(s(idx),g(idx));
fprintf(fid, '%e\n', g*1E3);

end 

%%% Write the correct receptor location
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Receptor Location"\n');
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', size(receptor_xy,1));
fprintf(fid, 'PASSIVEVARLIST = [3,4]\n');
fprintf(fid, '%e\n', xy(sensorIndex,1));
fprintf(fid, '%e\n', xy(sensorIndex,2));

%%% Write the correct source location
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Source Location"\n');
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', size(source_xy,1));
fprintf(fid, 'PASSIVEVARLIST = [3,4]\n');
fprintf(fid, '%e\n', xy(sourceIndex,1));
fprintf(fid, '%e\n', xy(sourceIndex,2));


fclose(fid);
srcPath = fullfile(dst_dir, [dst_name '.dat']);
dstPath = fullfile(dst_dir, [dst_name '.plt']);
system(['preplot ' srcPath ' ' dstPath]);

