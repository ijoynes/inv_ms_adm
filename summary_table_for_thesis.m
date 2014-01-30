clc
close
clear


cdf_threshold = [.999 .99 .95 .9];

beta = [0,1,10,100,1000,0,1,10,100,1000]';

run_dir = {fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_058');
    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_059');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_060');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_061');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_062');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_063');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_064');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_065');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_066');
   fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_067'); };

% run_dir = {fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_031');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_068');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_069');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_070');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_071');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_046');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_055');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_047');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_045');
%    fullfile('C:', 'Users', 'ijoynes', 'Documents', 'run_048') };

load(fullfile(run_dir{1},'Domain.mat'))
xy = xy - ones(nNodes,1)*min(xy);
lb = zeros(nNodes,1);
ub = Inf*ones(nNodes,1);
nbd = ones(nNodes,1);

% precomute spatial integration vector s2m
s2m = zeros(nNodes,1); 
for i = 1 : nTris
    s2m(tri(i,:)) = s2m(tri(i,:)) + det([ones(3,1) xy(tri(i,:),:)]);
end
s2m = s2m/6;    

ySum = zeros(nNodes,1);
tf = 900;
for i = tf :-1: 0
    load(fullfile('C:\Users\ijoynes\Documents\field_of_view_study\field_of_view', ['field_of_view_' int2str(i) '.mat']));
    y(y<0) = 0;
    y(xy(:,1) < 200) = 0;
    y(xy(:,1) > 1200) = 0;
    y(xy(:,2) < 200) = 0;
    y(xy(:,2) > 1200) = 0;
    
    %fprintf(fid, '\nZONE\n');
    %fprintf(fid, 'T = "Retroplume: t = %d"\n', tf-i);
    %fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
    %fprintf(fid, 'NODES = %d\n', nNodes);
    %fprintf(fid, 'ELEMENTS = %d\n', nTris);
    %fprintf(fid, 'PASSIVEVARLIST = [3-8,10,11]\n');
    
    
    % fprintf(fid, 'VARSHARELIST = ([1,2]=1)\n');
    % fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
    
    
    %fprintf(fid, 'STRANDID = 1\n');
    %fprintf(fid, 'SOLUTIONTIME = %d\n', tf-i);
    
    
    
    %fprintf(fid, '%e\n', y);
    
    
    
    if i == 900 || i == 0
        ySum = ySum + 7*y;
    elseif mod(i,4) == 1 || mod(i,4) == 3
        ySum = ySum +32*y;
    elseif mod(i,4) == 2
        ySum = ySum +12*y;
    elseif mod(i,4) == 0
        ySum = ySum +14*y;
    end
end

yInt4 = ySum*4/90;

P = yInt4/dot(s2m,yInt4);

range = unique(sort(P));
n = length(range);
percentile = nan(n,1);
for i = 1 : n
    P_temp = P;
    P_temp(P_temp<range(i)) = 0;
    percentile(i) = dot(s2m,P_temp);
end
percentile_field = interp1(range,percentile,P);   




load(fullfile(run_dir{1}, 'Source', 'Source_Correct.mat'))
Q_star = dot(s2m,E);

load(fullfile(run_dir{1}, 'Source', 'Source_0.mat'))
f0 = f;
g0 = g;
s0 = s;


for i = 1 :10
    i

load(fullfile(run_dir{i}, 'Source','Source_Correct.mat'),'signal_o')
c_o = signal_o;
c_o_mean = mean(mean(c_o));
load(fullfile(run_dir{1}, 'Source','Source_Correct.mat'),'signal_o')
c = signal_o;

r_max = 1-sum(sum((c-c_o).^2))/sum(sum((c_o_mean-c_o).^2));

    
load(fullfile(run_dir{i},'Source','Source_200.mat'))
load(fullfile(run_dir{1},'Source','Source_Correct.mat'),'E')

fprintf('f/f0 = %e\n', f/f0);
fprintf('f_obs/f0 = %e\n', (f-beta(i)/2*dot(s2m,s.^2))/f0);
fprintf('g_norm = %e\n', norm(g,Inf)/norm(g0,Inf));
fprintf('g_norm = %e\n', norm(projgr(g,s,lb,ub,nbd),Inf)/norm(projgr(g0,s0,lb,ub,nbd),Inf));
fprintf('r2 = %e\n', r2);
fprintf('r2/r_max = %e\n', r2/r_max);
fprintf('Q/Q* = %e\n', dot(s2m,s)/dot(s2m,E))
for j = 1 : 4
load(fullfile(run_dir{1},'Source','Source_Correct.mat'),'E')
load(fullfile(run_dir{i},'Source','Source_200.mat'),'s')
E(percentile_field>=cdf_threshold(j)) = 0;
s(percentile_field>=cdf_threshold(j)) = 0;
%m_star_filter = dot(s2m,E);
%m_star_filter*1E3
fprintf('Q_filter/Q*_filter = %e (Q*_filter = %e)\n', dot(s2m,s)/dot(s2m,E), dot(s2m,E));

end
end


