clc
close
clear


cdf_threshold = [.999 .99 .95 .9];

beta = [0,1,10,100,1000,0,1,10,100,1000]';

%run_dir = {fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_058');
%    fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_059');
%   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_060');
%   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_061');
%   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_062');
%   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_063');
%   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_064');
%   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_065');
%   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_066');
%   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_067'); };

run_dir = {fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_031');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_068');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_069');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_070');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_071');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_046');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_055');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_047');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_045');
   fullfile('E:', 'ijoynes', 'thesis_data_backup', 'run_048') };

run_name = {'No noise, No reg';
    'No noise, theta = 1';
    'No noise, theta = 10';
    'No noise, theta = 100';
    'No noise, theta = 1000';
    '10% noise, No reg';
    '10% noise, theta = 1'
    '10% noise, theta = 10';
    '10% noise, theta = 100';
    '10% noise, theta = 1000'};

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

Q = nan(201,10,5);

for i = 1 :10
  i
  
  for iter = 0 : 200
    load(fullfile(run_dir{i},'Source',['Source_' num2str(iter) '.mat']))
    load(fullfile(run_dir{1},'Source','Source_Correct.mat'),'E')
    Q(iter+1,i,1) = dot(s2m,s)/dot(s2m,E);

    
    %fprintf('Q/Q* = %e\n', dot(s2m,s)/dot(s2m,E))

    for j = 1 : 4
      load(fullfile(run_dir{1},'Source','Source_Correct.mat'),'E')
      load(fullfile(run_dir{i},'Source',['Source_' num2str(iter) '.mat']),'s')
      E(percentile_field>=cdf_threshold(j)) = 0;
      s(percentile_field>=cdf_threshold(j)) = 0;
      %m_star_filter = dot(s2m,E);
      %m_star_filter*1E3
      %fprintf('Q_filter/Q*_filter = %e (Q*_filter = %e)\n', dot(s2m,s)/dot(s2m,E), dot(s2m,E));
      Q(iter+1,i,j+1) = dot(s2m,s)/dot(s2m,E);
    end

  end

end

figure(1)
plot(0:200,Q(:,6,1))
grid

figure(2)
semilogy(1:200,abs(diff(Q(:,6,1))./Q(2:201,6,1)))
grid

temp = nan(201,10,5);
iDelta = 25;
conv_thresh = 5E-3;
for i = iDelta+1:201
    for iRun = 1 : 10
        for iFilter = 1 : 5
            %temp(i,iRun,iFilter) = (Q(i,iRun,iFilter)-Q(i-iDelta,iRun,iFilter))/min(Q(i,iRun,iFilter),Q(i-iDelta,iRun,iFilter));
            %temp(i,iRun,iFilter) = (Q(i,iRun,iFilter)-Q(i-iDelta,iRun,iFilter));
            temp(i,iRun,iFilter) = std(Q(i-iDelta:i,iRun,iFilter))/Q(i,iRun,iFilter);
        end
    end
end

for iFilter = 1 : 5
    
    for iRun = 1 : 5
        if abs(temp(200,iRun,iFilter)) <= conv_thresh
            fprintf(' %.2f ', Q(200,iRun,iFilter));
        else
            fprintf(' %.2f*', Q(200,iRun,iFilter));
        end
    end
    fprintf('\n');
end
fprintf('-----------------------------\n')
         
for iFilter = 1 : 5
    
    for iRun = 6 : 10
        if abs(temp(200,iRun,iFilter)) <= conv_thresh
            fprintf(' %.2f ', Q(200,iRun,iFilter));
        else
            fprintf(' %.2f*', Q(200,iRun,iFilter));
        end
    end
    fprintf('\n');
end


iRun = 2
iFilter = 3;
figure(1)
%plot(0:200,Q(:,6:10,iFilter))
plot(0:200,Q(:,6:10,iFilter))
grid
iter = 0:200;
figure(2)
semilogy(iter,abs(temp(:,6:10,iFilter)))
grid



figure(1)
semilogy(iters,temp(:,6:10,1))
legend(run_name{1:5})
grid

