close
clear
clc
%%%% This may need to be changed to run_031
%src_dir = 'E:\ijoynes\thesis_data_backup\run_023';
src_dir = 'E:\ijoynes\thesis_data_backup\run_031';
%flow_src_dir = 'E:\ijoynes\thesis_data_backup\run_023';
dst_dir = 'E:\ijoynes\thesis_data_backup';
dst_name = 'Single_Source_Adjoint_Simulation_2013-12-11_full';

C_mu = 0.085;
Sc_t = 0.9;
nt = 900;

load(fullfile(src_dir,'Domain.mat'),'xy','tri','nNodes','nTris');
load(fullfile(src_dir,'Source.mat'));
load(fullfile(src_dir,'Sensors.mat'));
load(fullfile(src_dir,'Source', 'Source_Correct.mat'),'E','signal_o');


nSources = size(source_xy,1);
nSensors = size(receptor_xy,1);

source_xy = source_xy - ones(nSources,1)*min(xy);
receptor_xy = receptor_xy - ones(nSensors,1)*min(xy);

xy = xy - ones(nNodes,1)*min(xy);
sourceIndex = placeSensors(xy,source_xy);

fid = fopen([fullfile(dst_dir,dst_name) '.dat'],'w');

fprintf(fid, 'TITLE = "Single Source Forward Simulation"\n');
fprintf(fid, [ 'VARIABLES = "x [m]" ', ...  %1
                           '"y [m]" ', ...  %2
                           '"<math>6&$6</math><times>c</times> [mg<math>W</math>m<sup>-3</sup>]" ', ... %3
                           '"<greek>l</greek> [mg<math>W</math>s<math>W</math>m<sup>-3</sup>]" ', ... %4
                           '"t [min]"\n'] ); %5


n = nt;
load(fullfile(src_dir,'Adjoint', ['Adjoint_', num2str(n) '.mat']),'y');
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Adjoint: %d"\n', n);
mins = floor(n/60);
secs = mod(n, 60);

if mins < 10 && secs < 10 
  fprintf(fid, 'AUXDATA time_in_mins="0%d:0%d"\n', mins, secs);

elseif mins < 10 && secs >= 10 
  fprintf(fid, 'AUXDATA time_in_mins="0%d:%d"\n', mins, secs);

elseif mins >= 10 && secs < 10
  fprintf(fid, 'AUXDATA time_in_mins="%d:0%d"\n', mins, secs);

elseif mins >= 10 && secs >= 10
  fprintf(fid, 'AUXDATA time_in_mins="%d:%d"\n', mins, secs);
end

fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
fprintf(fid, 'NODES = %d\n', nNodes);
fprintf(fid, 'ELEMENTS = %d\n', nTris);
fprintf(fid, 'STRANDID = 1\n');
fprintf(fid, 'SOLUTIONTIME = %d\n', nt-n);
fprintf(fid, 'PASSIVEVARLIST = [3,5]\n');
fprintf(fid, '%e\n', xy(:,1));
fprintf(fid, '%e\n', xy(:,2));
fprintf(fid, '%e\n', y*1E6);
fprintf(fid, '%d %d %d\n', tri');


for n = nt-1 : -1 : 0
  n
  load(fullfile(src_dir,'Adjoint', ['Adjoint_', num2str(n) '.mat']),'y');
  fprintf(fid, '\nZONE\n');
  fprintf(fid, 'T = "Adjoint: %d"\n', n);
  mins = floor(n/60);
  secs = mod(n, 60);

  if mins < 10 && secs < 10 
    fprintf(fid, 'AUXDATA time_in_mins="0%d:0%d"\n', mins, secs);

  elseif mins < 10 && secs >= 10 
    fprintf(fid, 'AUXDATA time_in_mins="0%d:%d"\n', mins, secs);

  elseif mins >= 10 && secs < 10
    fprintf(fid, 'AUXDATA time_in_mins="%d:0%d"\n', mins, secs);

  elseif mins >= 10 && secs >= 10
    fprintf(fid, 'AUXDATA time_in_mins="%d:%d"\n', mins, secs);
  end

  fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
  fprintf(fid, 'NODES = %d\n', nNodes);
  fprintf(fid, 'ELEMENTS = %d\n', nTris);
  fprintf(fid, 'STRANDID = 1\n');
  fprintf(fid, 'SOLUTIONTIME = %d\n', nt-n);
  fprintf(fid, 'VARSHARELIST = ([1,2]=1)\n');
  fprintf(fid, 'PASSIVEVARLIST = [3,5]\n');
  fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
  fprintf(fid, '%e\n', y*1E6);

end

%%% Write the correct receptor location
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Receptor Location"\n');
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', size(receptor_xy,1));
fprintf(fid, 'PASSIVEVARLIST = [3-5]\n');
fprintf(fid, '%e\n', xy(sensorIndex,1));
fprintf(fid, '%e\n', xy(sensorIndex,2));

%%% Write the correct source location
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Source Location"\n');
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', size(source_xy,1));
fprintf(fid, 'PASSIVEVARLIST = [3-5]\n');
fprintf(fid, '%e\n', xy(sourceIndex,1));
fprintf(fid, '%e\n', xy(sourceIndex,2));


for i = 1 : nSensors
    fprintf(fid, '\nZONE\n');
    fprintf(fid, 'T = "R# %d"\n', i);
    fprintf(fid, 'ZONETYPE = ORDERED\n');
    fprintf(fid, 'I = %d\n', nt+1);
    fprintf(fid, 'PASSIVEVARLIST = [4]\n');
    fprintf(fid, '%e\n', ones(nt+1,1)*xy(sensorIndex(i),1));
    fprintf(fid, '%e\n', ones(nt+1,1)*xy(sensorIndex(i),2));
    fprintf(fid, '%e\n', -signal_o(:,i)*1E6);
    fprintf(fid, '%e\n', (0:nt)/60);
end

fclose(fid);
srcPath = fullfile(dst_dir, [dst_name '.dat']);
dstPath = fullfile(dst_dir, [dst_name '.plt']);
system(['preplot ' srcPath ' ' dstPath]);

