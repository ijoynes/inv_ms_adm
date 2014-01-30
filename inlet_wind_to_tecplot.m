%M = CSVREAD('C:\Users\ijoynes\Documents\wind_tower_data\20050521_775_815_0.csv');

%index = find(M(:,3)==2);
fid = fopen('tower_data_for_wind_rose3.txt','w');
%N = length(index);
N = 901;
nodeIndex = 9321;

for i = 1 : N
    load(['C:\Users\ijoynes\Documents\run_023\Flow\Flow_' num2str(i-1) '.mat'],'uv');
   
    u = uv(nodeIndex,1);
    v = uv(nodeIndex,2);
    wdir(i) = 270-atan2(v,u)*180/pi;
    if wdir(i) < 0
        wdir(i) = wdir(i) + 360;
    elseif wdir(i) >= 360
        wdir(i) = wdir(i) - 360;
    end
    wmag(i) = norm([u,v]);
end


fid = fopen('inlet_wind_dir_and_mag.dat','w');
fprintf(fid, 'TITLE = "Inlet Wind Direction and Magnitude"\n');
fprintf(fid, 'VARIABLES = "t [min]" "mag [m/s]" "dir [deg]"\n');
fprintf(fid, '\n');
fprintf(fid, 'ZONE\n');
fprintf(fid, 'I = %d\n', N);
fprintf(fid, 'ZONETYPE = Ordered\n');
fprintf(fid, '%e\n', (0:900)/60);
fprintf(fid, '%e\n', wmag);
fprintf(fid, '%e\n', wdir);
fclose(fid);


