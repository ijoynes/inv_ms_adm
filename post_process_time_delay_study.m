% This script is intended to read in the source predictions for a regularization
% study and write the results to a tecplot ascii file that will get converted to
% a tecplot binary file with the 'preplot' program
%
% Notes: the results will only contain iteration data.  It will not contain
%		 transient data such as velocity or turbulence field variables.
%
% Writen by:
% 	Ian Joynes
% 	Nov. 23, 2012
%
% Updates:
%

% reset the environment
clc
close
clear

outName = 'single_source_reg_study_2013-07-30'
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

% run_dir = {fullfile('home', 'ijoynes', 'run_031');
%     fullfile('home', 'ijoynes', 'run_068');
%     fullfile('home', 'ijoynes', 'run_069');
%     fullfile('home', 'ijoynes', 'run_070');
%     fullfile('home', 'ijoynes', 'run_071');
%     fullfile('home', 'ijoynes', 'run_046');
%     fullfile('home', 'ijoynes', 'run_055');
%     fullfile('home', 'ijoynes', 'run_047');
%     fullfile('home', 'ijoynes', 'run_045');
%     fullfile('home', 'ijoynes', 'run_048') };

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

theta = [0 1 10 100 1000 0 1 10 100 1000];

domainDir = run_dir{1};
outPath = data_dir
thresh_hold_factor = 0.1;
load( fullfile(domainDir, 'Parameters.mat'));
load( fullfile(domainDir, 'Domain.mat') );
load( fullfile(domainDir, 'Source.mat') );
load( fullfile(domainDir, 'Sensors.mat') );


nSources = size(source_xy,1);
nSensors = size(receptor_xy,1);

source_xy = source_xy - ones(nSources,1)*min(xy);
% AUX DATA
% density =
%    m* =
%   reg_par =
%   factr =
%  pgtol =
% noise =
% m =
% maxIter =
% tMax =

% adjust axis origin
receptor_xy = receptor_xy - ones(nSensors,1)*min(xy);

xy = xy - ones(nNodes,1)*min(xy);
sourceIndex = placeSensors(xy,source_xy);

receptor_obs   = nan(nSensors,nt);
receptor_mod_i = nan(nSensors,nt);
receptor_mod_f = nan(nSensors,nt);

lb = zeros(nNodes,1);
ub = Inf*ones(nNodes,1);
nbd = ones(nNodes,1);

% precomute spatial integration vector s2m
s2m = zeros(nNodes,1); 
for i = 1 : nTris
    s2m(tri(i,:)) = s2m(tri(i,:)) + det([ones(3,1) xy(tri(i,:),:)]);
end
s2m = s2m/6;    

%%% Write the correct source
fid = fopen(fullfile(outPath, [outName, '.dat']),'w');

fprintf(fid, 'TITLE = "regularization_study"\n');
fprintf(fid, ['VARIABLES = "x [m]" ', ...
    '"y [m]" ', ...
    '"s [mg<math>W</math>s<sup>-1</sup><math>W</math>m<sup>-3</sup>]" ', ...
    '"Iteration #" ', ...
    '"f_i/f_0" ', ... %<math>&</math><sub>i</sub>/<math>&</math><sub>0</sub>" ', ...   % objective function
    '"||g_i||_inf/||g_0||_inf" ', ... %"<math>wwQ&</math><sub>i</sub><math>ww<sub>%</sub></math>/<math>wwQ&</math><sub>0</sub><math>ww<sub>%</sub></math>" ', ... % gradient
    '"m<sub>i</sub>/m*" ', ...
    '"r<sup>2</sup>" ',...
    '"c" ', ...
    '"pdf" ', ...
    '"cdf" ', ...
    '"t [min]" ', ...
    '"c [mg<math>W</math>m<sup>-3</sup>]"\n']);

load(fullfile(domainDir,'Source','Source_Correct.mat'),'E')
m_star = dot(s2m,E);
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Source Correct"\n');
fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
fprintf(fid, 'NODES = %d\n', nNodes);
fprintf(fid, 'ELEMENTS = %d\n', nTris);
fprintf(fid, 'PASSIVEVARLIST = [4-13]\n');

fprintf(fid, '%e\n', xy(:,1));
fprintf(fid, '%e\n', xy(:,2));
fprintf(fid, '%e\n', E*1E6);

fprintf(fid, '%d %d %d\n', tri');

%%% Write the correct receptor location
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Receptor Location"\n');
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', size(receptor_xy,1));
fprintf(fid, 'PASSIVEVARLIST = [3-13]\n');
fprintf(fid, '%e\n', xy(sensorIndex,1));
fprintf(fid, '%e\n', xy(sensorIndex,2));

%%% Write the correct source location
fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "Source Location"\n');
fprintf(fid, 'ZONETYPE = ORDERED\n');
fprintf(fid, 'I = %d\n', size(source_xy,1));
fprintf(fid, 'PASSIVEVARLIST = [3-13]\n');
fprintf(fid, '%e\n', xy(sourceIndex,1));
fprintf(fid, '%e\n', xy(sourceIndex,2));


%%% Source iteration with no noise and no regularization
% pre-determine number of interations
% determine the number of iterations in the directory

iRun = 1;
nIters = zeros(nRuns,1);
for i = 1 : nRuns
    while 1
        srcPath = fullfile(run_dir{i}, 'Source', ['Source_' int2str(nIters(i)) '.mat']);
        if exist(srcPath )~= 2
            break;
        end
        nIters(i) = nIters(i) + 1;
    end
end

f_all = nan(max(nIters),nRuns);
f_all_obs = nan(max(nIters),nRuns);
f_all_min = nan(max(nIters),nRuns);
g_norm_all = nan(max(nIters), nRuns);
g_norm_pr = nan(max(nIters), nRuns);
m_all = nan(max(nIters), nRuns);
r2_all = nan(max(nIters),nRuns);
r2_max = nan(max(nIters),nRuns);

B = spdiags(s2m,0,nNodes,nNodes);

for iRun = 1 : nRuns
    iter = 0 : (nIters(iRun)-1);
    for i = 1 : nIters(iRun)
        load( fullfile(run_dir{iRun}, 'Source', ['Source_Correct.mat']), 'signal_o')
        signal_o_temp = signal_o;
        o_mean = mean(mean(signal_o_temp));  % IS THIS THE CORRECT AVERAGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        load( fullfile(run_dir{1}, 'Source', ['Source_Correct.mat']), 'signal_o')
        load( fullfile(run_dir{iRun}, 'Source', ['Source_' num2str(iter(i)) '.mat']), 's','f','g','r2')
        f_all(i,iRun) = f;
        f_all_obs(i,iRun) = f-theta(iRun)/2*(s'*(B*s));
        f_all_min(i,iRun) = sum(sum((signal_o-signal_o_temp).^2))/2;
        m_all(i,iRun) = dot(s2m,s);
        g_norm_all(i,iRun) = norm(g,Inf);
        g_norm_pr(i,iRun) = norm(projgr(s, g, lb, ub, nbd),Inf);
        r2_all(i,iRun) = r2;
        %r2_max(i,iRun) = 1 - sum(sum((signal_o-signal_o_temp).^2))/sum(sum((signal_o-o_mean).^2));
        r2_max(i,iRun) = 1 - sum(sum((signal_o-signal_o_temp).^2))/sum(sum((signal_o_temp-o_mean).^2));
        fprintf(fid, '\nZONE\n');
        fprintf(fid, 'T = "%s (I=%d)"\n',run_name{iRun}, iter(i));
        fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
        fprintf(fid, 'NODES = %d\n', nNodes);
        fprintf(fid, 'ELEMENTS = %d\n', nTris);
        fprintf(fid, 'STRANDID = %d\n', iRun);
        fprintf(fid, 'SOLUTIONTIME = %d\n', iter(i));
        
        fprintf(fid, 'PASSIVEVARLIST = [4-13]\n');
        fprintf(fid, 'VARSHARELIST = ([1,2]=1)\n');
        fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
        
        fprintf(fid, '%e\n', s*1E6);
    end
end


for iRun = 1 : nRuns
    fprintf(fid, '\nZONE\n');
    fprintf(fid, 'T = "CONVERGENCE (%s)"\n',run_name{iRun} );
    fprintf(fid, 'ZONETYPE = ORDERED\n');
    fprintf(fid, 'I = %d\n', nIters(iRun));
    fprintf(fid, 'PASSIVEVARLIST = [1-3,9-13]\n');
    fprintf(fid, '%d\n', 0:(nIters(iRun)-1));
    fprintf(fid, '%e\n', f_all(1:nIters(iRun),iRun)/f_all(1,1));
    fprintf(fid, '%e\n', g_norm_all(1:nIters(iRun),iRun)/g_norm_all(1,1));
    fprintf(fid, '%e\n', m_all(1:nIters(iRun),iRun)/m_star);
    fprintf(fid, '%e\n', r2_all(1:nIters(iRun),iRun));
end

for iRun = 1 : nRuns
    fprintf(fid, '\nZONE\n');
    fprintf(fid, 'T = "CONVERGENCE: OBS (%s)"\n',run_name{iRun} );
    fprintf(fid, 'ZONETYPE = ORDERED\n');
    fprintf(fid, 'I = %d\n', nIters(iRun));
    fprintf(fid, 'PASSIVEVARLIST = [1-3,9-13]\n');
    fprintf(fid, '%d\n', 0:(nIters(iRun)-1));
    fprintf(fid, '%e\n', f_all_obs(1:nIters(iRun),iRun)/f_all_obs(1,1));
    fprintf(fid, '%e\n', g_norm_pr(1:nIters(iRun),iRun)/g_norm_pr(1,1));
    fprintf(fid, '%e\n', m_all(1:nIters(iRun),iRun)/m_star);
    fprintf(fid, '%e\n', r2_all(1:nIters(iRun),iRun));
end

for iRun = 1 : nRuns
    fprintf(fid, '\nZONE\n');
    fprintf(fid, 'T = "CONVERGENCE: LIMIT (%s)"\n',run_name{iRun} );
    fprintf(fid, 'ZONETYPE = ORDERED\n');
    fprintf(fid, 'I = %d\n', nIters(iRun));
    fprintf(fid, 'PASSIVEVARLIST = [1-3,6,7,9-13]\n');
    fprintf(fid, '%d\n', 0:(nIters(iRun)-1));
    fprintf(fid, '%e\n', f_all_min(1:nIters(iRun),iRun)/f_all_obs(1,1));
    fprintf(fid, '%e\n', r2_max(1:nIters(iRun),iRun));
end

for iRun = 1 : nRuns
    fprintf(fid, '\nZONE\n');
    fprintf(fid, 'T = "CONVERGENCE: ADJUSTED (%s)"\n',run_name{iRun} );
    fprintf(fid, 'ZONETYPE = ORDERED\n');
    fprintf(fid, 'I = %d\n', nIters(iRun));
    fprintf(fid, 'PASSIVEVARLIST = [1-3,6,7,9-13]\n');
    fprintf(fid, '%d\n', 0:(nIters(iRun)-1));
    fprintf(fid, '%e\n', f_all_obs(1:nIters(iRun),iRun)/f_all_obs(1,1)-f_all_min(1:nIters(iRun),iRun)/f_all_obs(1,1));
    fprintf(fid, '%e\n', r2_all(1:nIters(iRun),iRun)./r2_max(1:nIters(iRun),iRun));
end



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




cdf_threshold = [.999 .99 .95 .9];

for j = 1 : length(cdf_threshold)
    for iRun = 1 : nRuns
        iter = 0 : (nIters(iRun)-1);
        for i = 1 : nIters(iRun)
            load( fullfile(run_dir{iRun}, 'Source', ['Source_' num2str(iter(i)) '.mat']), 's')
            %f_all(i,iRun) = f;
            s(percentile_field>=cdf_threshold(j)) = 0;
            m_all(i,iRun) = dot(s2m,s);
            %g_norm_all(i,iRun) = max(abs(g));
            %r2_all(i,iRun) = r2;
            fprintf(fid, '\nZONE\n');
            fprintf(fid, 'T = "%s, CDF=%f (I=%d)"\n',run_name{iRun}, cdf_threshold(j), iter(i));
            fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
            fprintf(fid, 'NODES = %d\n', nNodes);
            fprintf(fid, 'ELEMENTS = %d\n', nTris);
            fprintf(fid, 'STRANDID = %d\n', iRun+j*nRuns);
            fprintf(fid, 'SOLUTIONTIME = %d\n', iter(i));
            
            fprintf(fid, 'PASSIVEVARLIST = [4-13]\n');
            fprintf(fid, 'VARSHARELIST = ([1,2]=1)\n');
            fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
            
            fprintf(fid, '%e\n', s*1E6);
        end
    end
    for iRun = 1 : nRuns
        fprintf(fid, '\nZONE\n');
        fprintf(fid, 'T = "CONVERGENCE (%s) CDF=%f"\n',run_name{iRun},cdf_threshold(j) );
        fprintf(fid, 'ZONETYPE = ORDERED\n');
        fprintf(fid, 'I = %d\n', nIters(iRun));
        fprintf(fid, 'PASSIVEVARLIST = [1-3,5,6,8,9-13]\n');
        fprintf(fid, '%d\n', 0:(nIters(iRun)-1));
        fprintf(fid, '%e\n', m_all(1:nIters(iRun),iRun)/m_star);
    end
end


fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "PDF"\n');
fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
fprintf(fid, 'NODES = %d\n', nNodes);
fprintf(fid, 'ELEMENTS = %d\n', nTris);
fprintf(fid, 'PASSIVEVARLIST = [3-9,11-13]\n');
fprintf(fid, 'VARSHARELIST = ([1,2]=1)\n');
fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
fprintf(fid, '%e\n', P);

fprintf(fid, '\nZONE\n');
fprintf(fid, 'T = "CDF"\n');
fprintf(fid, 'ZONETYPE = FETRIANGLE\n');
fprintf(fid, 'NODES = %d\n', nNodes);
fprintf(fid, 'ELEMENTS = %d\n', nTris);
fprintf(fid, 'PASSIVEVARLIST = [3-10,12,13]\n');
fprintf(fid, 'VARSHARELIST = ([1,2]=1)\n');
fprintf(fid, 'CONNECTIVITYSHAREZONE = 1\n');
%fprintf(fid, '%e\n', percentile_field);
fprintf(fid, '%e\n', 1-percentile_field);


for iRun = 1 : nRuns
    iter = 0 : (nIters(iRun)-1);
    for j = 1 : nSensors
        load( fullfile(run_dir{iRun}, 'Source', 'Source_Correct.mat'), 'signal_o')

            fprintf(fid, '\nZONE\n');
            fprintf(fid, 'T = "R# %d, %s, OBS"\n',j,run_name{iRun});
            fprintf(fid, 'ZONETYPE = ORDERED\n');
            fprintf(fid, 'I = %d\n', nt);
            
            fprintf(fid, 'PASSIVEVARLIST = [3-11]\n');
            fprintf(fid, '%e\n', xy(sensorIndex(j),1)*ones(nt,1));
            fprintf(fid, '%e\n', xy(sensorIndex(j),2)*ones(nt,1));
            
            fprintf(fid, '%e\n', t/60);
            fprintf(fid, '%e\n', signal_o(:,j)*1E6);
        for i = 1 : nIters(iRun)
            load( fullfile(run_dir{iRun}, 'Source', ['Source_' num2str(iter(i)) '.mat']), 'signal_c')

            fprintf(fid, '\nZONE\n');
            fprintf(fid, 'T = "R# %d, %s, (I=%d)"\n',j,run_name{iRun}, iter(i));
            fprintf(fid, 'ZONETYPE = ORDERED\n');
            fprintf(fid, 'I = %d\n', nt);
            
            fprintf(fid, 'PASSIVEVARLIST = [3-11]\n');
            fprintf(fid, '%e\n',  xy(sensorIndex(j),1)*ones(nt,1));
            fprintf(fid, '%e\n',  xy(sensorIndex(j),2)*ones(nt,1));
            
            fprintf(fid, '%e\n', t/60);
            fprintf(fid, '%e\n', signal_c(:,j)*1E6);
        end
    end
end

fclose(fid);
srcPath = fullfile(outPath, [outName '.dat']);
dstPath = fullfile(outPath, [outName '.plt']);
system(['preplot ' srcPath ' ' dstPath]);







