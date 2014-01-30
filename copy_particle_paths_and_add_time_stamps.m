close
clear
clc

src_dir = 'C:\Users\ijoynes\Documents\run_023';
src_file = 'particle_paths_2013-06-21.dat';

dst_dir = 'E:\ijoynes\thesis_data_backup';
dst_file = 'particle_paths_2013-12-04.dat';

fid_in = fopen(fullfile(src_dir, src_file),'r');
fid_out = fopen(fullfile(dst_dir, dst_file),'w');

temp = fgetl(fid_in);
fprintf(fid_out, '%s\n', temp);
temp = fgetl(fid_in);
fprintf(fid_out, 'VARIABLES = "x [m]" "y [m]" "u [m<math>W</math>s<sup>-1</sup>]" "v [m<math>W</math>s<sup>-1</sup>]" "D [m<sup>2</sup><math>W</math>s<sup>-1</sup>]"\n');
i = 0
while ~feof(fid_in)
    if i > 900
        i = 0
    end
    temp = fgetl(fid_in);
    if strcmp(temp, 'T = "Flow"') || strcmp(temp, 'T = "Particle"')
        
          mins = floor(i/60);
    secs = mod(i,60);
    
    if mins < 10
        time_in_mins = ['0', num2str(mins), ':'];
    else
        time_in_mins = [num2str(mins), ':'];
    end
    
    
    if secs < 10
        time_in_mins = [time_in_mins, '0', num2str(secs)];
    else
        time_in_mins = [time_in_mins, num2str(secs)];
    end
    fprintf(fid_out, 'AUXDATA time_in_mins="%s"\n', time_in_mins);
    i = i+1
    end
    fprintf(fid_out, '%s\n', temp);
end


fclose(fid_in);
fclose(fid_out);