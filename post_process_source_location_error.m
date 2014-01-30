clc
close
clear

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

% run_dir = {fullfile(data_dir, 'run_058');
%     fullfile(data_dir, 'run_059');
%     fullfile(data_dir, 'run_060');
%     fullfile(data_dir, 'run_061');
%     fullfile(data_dir, 'run_062');
%     fullfile(data_dir, 'run_063');
%     fullfile(data_dir, 'run_064');
%     fullfile(data_dir, 'run_065');
%     fullfile(data_dir, 'run_066');
%     fullfile(data_dir, 'run_067'); };

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
load( fullfile(domainDir, 'Domain.mat') );
xy =xy - ones(nNodes,1)*min(xy);
s2m = zeros(nNodes,1); 
for i = 1 : nTris
    s2m(tri(i,:)) = s2m(tri(i,:)) + det([ones(3,1) xy(tri(i,:),:)]);
end
s2m = s2m/6;  

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

load(fullfile(run_dir{1}, 'Source', 'Source_Correct.mat'));
idxE = find(E>0);
src_label = 'HJCFBIAEDG';
for iRun = 1 : 10
    fprintf('iRun = %d\n', iRun);
   fprintf('--------------------------------------------------\n') 
   fprintf('| S | x* [m] | y* [m] |  x [m] |  y [m] |  r [m] |\n')
   fprintf('|---|--------|--------|--------|--------|--------|\n')

for iSource = 1 : length(idxE)
    
    %fprintf('| %s | %6f | %6f | %6f | %6f | %6f |\n',char(64+iSource),)
    load(fullfile(run_dir{iRun}, 'Source', 'Source_200.mat'));
    s( sum((xy-ones(nNodes,1)*xy(idxE(iSource),:)).^2,2) > 40^2 ) = 0;
    [vs,idxs] = max(s);
% sqrt(sum((xy(idxE,:) - xy(idxs,:)).^2,2))
xy_s_i = xy(idxs,:);
r = sqrt(sum((xy(idxE(iSource),:) - xy_s_i).^2,2));
fprintf('| %s | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f |\n',src_label(iSource),xy(idxE(iSource),1),xy(idxE(iSource),2),xy_s_i(1),xy_s_i(2),r)
end
fprintf('--------------------------------------------------\n\n') 

end


load(fullfile(run_dir{1}, 'Source', 'Source_Correct.mat'));
idxE = find(E>0);
src_label = 'HJCFBIAEDG';
   fprintf('NO OBSERVATIONAL NOISE\n')
   fprintf('--------------------------------------------------\n') 
   fprintf('| S |   r_0  |   r_1  |  r_10  |  r_100 | r_1000 |\n')
   fprintf('|:-:|:------:|:------:|:------:|:------:|:------:|\n')
    for iSource = 1 : length(idxE)
    
   fprintf('| %s',src_label(iSource));

for iRun = 1 :5
    
    %fprintf('| %s | %6f | %6f | %6f | %6f | %6f |\n',char(64+iSource),)
    load(fullfile(run_dir{iRun}, 'Source', 'Source_200.mat'));
    s( sum((xy-ones(nNodes,1)*xy(idxE(iSource),:)).^2,2) > 25^2 ) = 0;
    [vs,idxs] = max(s);
% sqrt(sum((xy(idxE,:) - xy(idxs,:)).^2,2))
xy_s_i = xy(idxs,:);
r = sqrt(sum((xy(idxE(iSource),:) - xy_s_i).^2,2));
fprintf(' | %6.2f',r)
end
fprintf(' |\n')


    end
fprintf('--------------------------------------------------\n\n') 


fprintf('10%% OBSERVATIONAL NOISE\n')
   fprintf('--------------------------------------------------\n') 
   fprintf('| S |   r_0  |   r_1  |  r_10  |  r_100 | r_1000 |\n')
   fprintf('|:-:|:------:|:------:|:------:|:------:|:------:|\n')
    for iSource = 1 : length(idxE)
    
   fprintf('| %s',src_label(iSource));

for iRun = 6 :10
    
    %fprintf('| %s | %6f | %6f | %6f | %6f | %6f |\n',char(64+iSource),)
    load(fullfile(run_dir{iRun}, 'Source', 'Source_200.mat'));
    s( sum((xy-ones(nNodes,1)*xy(idxE(iSource),:)).^2,2) > 25^2 ) = 0;
    [vs,idxs] = max(s);
% sqrt(sum((xy(idxE,:) - xy(idxs,:)).^2,2))
xy_s_i = xy(idxs,:);
r = sqrt(sum((xy(idxE(iSource),:) - xy_s_i).^2,2));
fprintf(' | %6.2f',r)
end
fprintf(' |\n')


    end
fprintf('--------------------------------------------------\n\n') 



%% Find the area of half max volumetric emission rate for each source
area = nan(10,10);
num_elems = nan(10,10);
for iRun = 1 :10
%iRun = 10;

%src_radius = [35, 60, 100, 60, 60, 35, 40, 70, 60, 70];

load(fullfile(run_dir{1}, 'Source', 'Source_Correct.mat'));

idxE = find(E>0);
src_label = 'HJCFBIAEDG';

for iSource = 1 :1
    load(fullfile(run_dir{iRun}, 'Source', 'Source_200.mat'));
    %load(fullfile(run_dir{iRun}, 'Source', 'Source_Correct.mat'));
    %s = E;
s( sum((xy-ones(nNodes,1)*xy(idxE(iSource),:)).^2,2) > 25^2 ) = 0;
[vs,idxs] = max(s);
load(fullfile(run_dir{iRun}, 'Source', 'Source_200.mat'));
%load(fullfile(run_dir{iRun}, 'Source', 'Source_Correct.mat'));
%s = E;
s( s> vs ) = 0;
s( s< vs/2 ) = 0;
%s( sum((xy-ones(nNodes,1)*xy(idxs,:)).^2,2) > src_radius(iSource)^2) = 0;



    
%trisurf(tri,xy(:,1),xy(:,2),s,'edgecolor','interp','facecolor','interp')
%view(2)
%axis image
%colorbar

area(iSource,iRun) = 0;
num_elems(iSource,iRun) = 0;
for i = 1 : nTris
    if any(s(tri(i,:)) > 0)
        area(iSource,iRun) = area(iSource,iRun) + det([ones(3,1), xy(tri(i,:),:)])/2;
        num_elems(iSource,iRun) = num_elems(iSource,iRun) + 1;
    end
end
%src_label(iSource)
end
   
end

% Find all nodes within a reasonable distance from the candidate source
%   that are above the half max.

fprintf('------------------------------------------------------------\n')
fprintf('| S |   B = 0  |   B = 1  |  B = 10  |  B = 100 | B = 1000 |\n')

for i = 1 : 1
    fprintf('|:-:|:--------:|:--------:|:--------:|:--------:|:--------:|\n')
    fprintf('| %s', src_label(i))
    for j = 1 : 5
        fprintf(' | %8.0f', area(i,j))
    end
    fprintf(' |\n');
    fprintf('|  '); 
    for j = 1 : 5
        fprintf(' | %8d', num_elems(i,j))
    end
    fprintf(' |\n');
    
end
fprintf('------------------------------------------------------------\n')
fprintf('\n');
fprintf('------------------------------------------------------------\n')
fprintf('| S |   B = 0  |   B = 1  |  B = 10  |  B = 100 | B = 1000 |\n')

for i = 1 : 1
    fprintf('|:-:|:--------:|:--------:|:--------:|:--------:|:--------:|\n')
    fprintf('| %s', src_label(i))
    for j = 6 : 10
        fprintf(' | %8.0f', area(i,j))
    end
    fprintf(' |\n');
    fprintf('|  '); 
    for j = 6 : 10
        fprintf(' | %8d', num_elems(i,j))
    end
    fprintf(' |\n');
    
end
fprintf('------------------------------------------------------------\n')
fprintf('\n')




