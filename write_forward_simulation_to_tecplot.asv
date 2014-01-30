close
clear
clc
%%%% This may need to be changed to run_031
src_dir = 'E:\ijoynes\thesis_data_backup\run_031';
flow_src_dir = 'E:\ijoynes\thesis_data_backup\run_023';
dst_dir = 'E:\ijoynes\thesis_data_backup';
dst_name = 'Single_Source_Forward_Simulation_2013-12-11';

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
                           '"u [m<math>W</math>s<sup>-1</sup>]" ', ... %3
                           '"v [m<math>W</math>s<sup>-1</sup>]" ', ... %4
                           '"k [m<sup>2</sup><math>W</math>s<sup>-2</sup>]" ', ... %5
                           '"<greek>e</greek> [m<sup>2</sup><math>W</math>s<sup>-2</sup>]" ', ... %6
                           '"<greek>n</greek><sub>t</sub> [m<sup>2</sup><math>W</math>s<sup>-1</sup>]" ', ... %7
                           '"D [m<sup>2</sup><math>W</math>s<sup>-1</sup>]" ', ...  %8
                           '"s [mg<math>W</math>s<sup>-1</sup><math>W</math>m<sup>-3</sup>]" ', ... %9
                           '"c [mg<math>W</math>m<sup>-3</sup>]" ', ... %10
                           '"t [min]"\n'] ); %11


n = 0;
load(fullfile(flow_src_dir,'Flow', ['Flow_', num2str(n) '.mat']),'uv','k','epsilion');
load(fullfile(src_dir,'Observation', ['Observation_', num2str(n) '.mat']),'o');
fprintf(fid, '\nZONE\n');
fprintf(fid, '\n');
fprintf(fid, 'T = "Flow: %d"\n', n);
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
fprintf(fid, 'SOLUTIONTIME = %d\n', n);
fprintf(fid, 'PASSIVEVARLIST = [11]\n');
fprintf(fid, '%e\n', xy(:,1));
fprintf(fid, '%e\n', xy(:,2));
fprintf(fid, '%e\n', uv(:,1));
fprintf(fid, '%e\n', uv(:,2));
fprintf(fid, '%e\n', k);
fprintf(fid, '%e\n', epsilion);
fprintf(fid, '%e\n', (C_mu/Sc_t)*(k.^2)./epsilion);
fprintf(fid, '%e\n', (C_mu/Sc_t)*(k.^2)./epsilion + 2.2E-5);
fprintf(fid, '%e\n', E*1E6);
fprintf(fid, '%e\n', o*1E6);
fprintf(fid, '%d %d %d\n', tri');


for n = 1 : nt
  n
  load(fullfile(flow_src_dir,'Flow', ['Flow_', num2str(n) '.mat']),'uv','k','epsilion');
  load(fullfile(src_dir,'Observation', ['Observation_', num2str(n) '.mat']),'o');
  fprintf(fid, '\nZONE\n');
  fprintf(fid, '\n');
  fprintf(fid, 'T = "Flow: %d"\n', n);
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
  fprintf(fid, 'SOLUTIONTIME = %d\n', n);
  fprintf(fid, 'VARSHARELIST = ([1,2,9]=1)\n');
  fprintf(fid, 'PASSIVEVARLIST = [11]\n');
  fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
  fprintf(fid, '%e\n', uv(:,1));
  fprintf(fid, '%e\n', uv(:,2));
  fprintf(fid, '%e\n', k);
  fprintf(fid, '%e\n', epsilion);
  fprintf(fid, '%e\n', (C_mu/Sc_t)*(k.^2)./epsilion);
  fprintf(fid, '%e\n', (C_mu/Sc_t)*(k.^2)./epsilion + 2.2E-5);
  fprintf(fid, '%e\n', o*1E6);

end

%%% Write the correct receptor location
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Receptor Location"\n');
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', size(receptor_xy,1));
fprintf(fid, 'PASSIVEVARLIST = [3-11]\n');
fprintf(fid, '%e\n', xy(sensorIndex,1));
fprintf(fid, '%e\n', xy(sensorIndex,2));

%%% Write the correct source location
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Source Location"\n');
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', size(source_xy,1));
fprintf(fid, 'PASSIVEVARLIST = [3-11]\n');
fprintf(fid, '%e\n', xy(sourceIndex,1));
fprintf(fid, '%e\n', xy(sourceIndex,2));


for i = 1 : nSensors
  i
  for n = 0 : nt
    fprintf(fid, '\nZONE\n');
    fprintf(fid, 'T = "R# %d: %d"\n', i, n);
    fprintf(fid, 'ZONETYPE = ORDERED\n');
    fprintf(fid, 'I = %d\n', n+1);
    fprintf(fid, 'PASSIVEVARLIST = [3-9]\n');
    fprintf(fid, '%e\n', ones(n+1,1)*xy(sensorIndex(i),1));
    fprintf(fid, '%e\n', ones(n+1,1)*xy(sensorIndex(i),2));
    fprintf(fid, '%e\n', signal_o(1:(n+1),i)*1E6);
    fprintf(fid, '%e\n', (0:n)/60);
  end
end

fclose(fid);
srcPath = fullfile(dst_dir, [dst_name '.dat']);
dstPath = fullfile(dst_dir, [dst_name '.plt']);
system(['preplot ' srcPath ' ' dstPath]);

