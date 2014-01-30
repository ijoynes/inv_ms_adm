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


f_all = nan(max(nIters),nRuns);
g_norm_all = nan(max(nIters), nRuns);
m_all = nan(max(nIters), nRuns);
r2_all = nan(max(nIters),nRuns);

df2 = nan( max(nIters), nRuns );
df4 = nan( max(nIters), nRuns );
df6 = nan( max(nIters), nRuns );
dg6 = nan( max(nIters), nRuns );
dm6 = nan( max(nIters), nRuns );
dr6 = nan( max(nIters), nRuns );

for iRun = 1 : nRuns
    iter = 0 : (nIters(iRun)-1);
    for i = 1 : nIters(iRun)

        load( fullfile(run_dir{iRun}, 'Source', ['Source_' num2str(iter(i)) '.mat']), 's','f','g','r2')
        f_all(i,iRun) = f;
        
        m_all(i,iRun) = dot(s2m,s);
        g_norm_all(i,iRun) = norm(g,Inf);
        
        r2_all(i,iRun) = r2;
    end
end

%df2
for iRun = 1 : nRuns
  for i = 1 : nIters(iRun)
    if i <= 1
      df2(i,iRun) = -1.5*f_all(i,iRun) + 2*f_all(i+1,iRun) -0.5*f_all(i+2,iRun);
    elseif i >= nIters(iRun)
      df2(i,iRun) = (-1.5*f_all(i,iRun) + 2*f_all(i-1,iRun) -0.5*f_all(i-2,iRun))/-1.0;
    else
      df2(i,iRun) = -0.5*f_all(i-1,iRun) + 0.5*f_all(i+1,iRun);
    end
  end
end
df2=df2/f(1,1);


%df4
for iRun = 1 : nRuns
  for i = 1 : nIters(iRun)
    if i <= 2
      df4(i,iRun) = -25.0/12.0*f_all(i,iRun) + 4.0*f_all(i+1,iRun) -3.0*f_all(i+2,iRun) ...
                   +4.0/3.0*f_all(i+3,iRun) -1.0/4.0*f_all(i+4,iRun);
    elseif i >= nIters(iRun)-1
      df4(i,iRun) = (-25.0/12.0*f_all(i,iRun) + 4.0*f_all(i-1,iRun) -3.0*f_all(i-2,iRun) ...
                   +4.0/3.0*f_all(i-3,iRun) -1.0/4.0*f_all(i-4,iRun))/-1.0;
    else
      df4(i,iRun) = 1.0/12.0*f_all(i-2,iRun)-2.0/3.0*f_all(i-1,iRun)+2.0/3.0*f_all(i+1,iRun) ...
                   -1.0/12.0*f_all(i+2,iRun);
    end
  end
end
df4=df4/f(1,1);

%df6
for iRun = 1 : nRuns
  for i = 1 : nIters(iRun)
    if i <= 3
      df6(i,iRun) = -49.0/20.0 * f_all(i,iRun)   +      6.0 * f_all(i+1,iRun) ...
                    -15.0/2.0  * f_all(i+2,iRun) + 20.0/3.0 * f_all(i+3,iRun) ...
                    -15.0/4.0  * f_all(i+4,iRun) +  6.0/5.0 * f_all(i+5,iRun)  ...
                    - 1.0/6.0  * f_all(i+6,iRun);
    elseif i >= nIters(iRun)-2
      df6(i,iRun) = (-49.0/20.0 * f_all(i,iRun)   +      6.0 * f_all(i-1,iRun) ...
                    -15.0/2.0  * f_all(i-2,iRun) + 20.0/3.0 * f_all(i-3,iRun) ...
                    -15.0/4.0  * f_all(i-4,iRun) +  6.0/5.0 * f_all(i-5,iRun)  ...
                    - 1.0/6.0  * f_all(i-6,iRun))/-1.0;
    else
      df6(i,iRun) = -1.0/60.0 * f_all(i-3, iRun) + 3.0/20.0 * f_all(i-2, iRun) ...
                    -3.0/4.0 * f_all(i-1, iRun) + 1.0/60.0 * f_all(i+3, iRun) ...
                    -3.0/20.0 * f_all(i+2, iRun) + 3.0/4.0  * f_all(i+1, iRun);
    end
  end
end
df6=df6/f(1,1);

%dg6
for iRun = 1 : nRuns
  for i = 1 : nIters(iRun)
    if i <= 3
      dg6(i,iRun) = -49.0/20.0 * g_norm_all(i,iRun)   +      6.0 * g_norm_all(i+1,iRun) ...
                    -15.0/2.0  * g_norm_all(i+2,iRun) + 20.0/3.0 * g_norm_all(i+3,iRun) ...
                    -15.0/4.0  * g_norm_all(i+4,iRun) +  6.0/5.0 * g_norm_all(i+5,iRun)  ...
                    - 1.0/6.0  * g_norm_all(i+6,iRun);
    elseif i >= nIters(iRun)-2
      dg6(i,iRun) = (-49.0/20.0 * g_norm_all(i,iRun)   +      6.0 * g_norm_all(i-1,iRun) ...
                    -15.0/2.0  * g_norm_all(i-2,iRun) + 20.0/3.0 * g_norm_all(i-3,iRun) ...
                    -15.0/4.0  * g_norm_all(i-4,iRun) +  6.0/5.0 * g_norm_all(i-5,iRun)  ...
                    - 1.0/6.0  * g_norm_all(i-6,iRun))/-1.0;
    else
      dg6(i,iRun) = -1.0/60.0 * g_norm_all(i-3, iRun) + 3.0/20.0 * g_norm_all(i-2, iRun) ...
                    -3.0/4.0 * g_norm_all(i-1, iRun) + 1.0/60.0 * g_norm_all(i+3, iRun) ...
                    -3.0/20.0 * g_norm_all(i+2, iRun) + 3.0/4.0  * g_norm_all(i+1, iRun);
    end
  end
end

%dr6
for iRun = 1 : nRuns
  for i = 1 : nIters(iRun)
    if i <= 3
      dr6(i,iRun) = -49.0/20.0 * r2_all(i,iRun)   +      6.0 * r2_all(i+1,iRun) ...
                    -15.0/2.0  * r2_all(i+2,iRun) + 20.0/3.0 * r2_all(i+3,iRun) ...
                    -15.0/4.0  * r2_all(i+4,iRun) +  6.0/5.0 * r2_all(i+5,iRun)  ...
                    - 1.0/6.0  * r2_all(i+6,iRun);
    elseif i >= nIters(iRun)-2
      dr6(i,iRun) = (-49.0/20.0 * r2_all(i,iRun)   +      6.0 * r2_all(i-1,iRun) ...
                    -15.0/2.0  * r2_all(i-2,iRun) + 20.0/3.0 * r2_all(i-3,iRun) ...
                    -15.0/4.0  * r2_all(i-4,iRun) +  6.0/5.0 * r2_all(i-5,iRun)  ...
                    - 1.0/6.0  * r2_all(i-6,iRun))/-1.0;
    else
      dr6(i,iRun) = -1.0/60.0 * r2_all(i-3, iRun) + 3.0/20.0 * r2_all(i-2, iRun) ...
                    -3.0/4.0 * r2_all(i-1, iRun) + 1.0/60.0 * r2_all(i+3, iRun) ...
                    -3.0/20.0 * r2_all(i+2, iRun) + 3.0/4.0  * r2_all(i+1, iRun);
    end
  end
end


%dm6
for iRun = 1 : nRuns
  for i = 1 : nIters(iRun)
    if i <= 3
      dm6(i,iRun) = -49.0/20.0 * m_all(i,iRun)   +      6.0 * m_all(i+1,iRun) ...
                    -15.0/2.0  * m_all(i+2,iRun) + 20.0/3.0 * m_all(i+3,iRun) ...
                    -15.0/4.0  * m_all(i+4,iRun) +  6.0/5.0 * m_all(i+5,iRun)  ...
                    - 1.0/6.0  * m_all(i+6,iRun);
    elseif i >= nIters(iRun)-2
      dm6(i,iRun) = (-49.0/20.0 * m_all(i,iRun)   +      6.0 * m_all(i-1,iRun) ...
                    -15.0/2.0  * m_all(i-2,iRun) + 20.0/3.0 * m_all(i-3,iRun) ...
                    -15.0/4.0  * m_all(i-4,iRun) +  6.0/5.0 * m_all(i-5,iRun)  ...
                    - 1.0/6.0  * m_all(i-6,iRun))/-1.0;
    else
      dm6(i,iRun) = -1.0/60.0 * m_all(i-3, iRun) + 3.0/20.0 * m_all(i-2, iRun) ...
                    -3.0/4.0 * m_all(i-1, iRun) + 1.0/60.0 * m_all(i+3, iRun) ...
                    -3.0/20.0 * m_all(i+2, iRun) + 3.0/4.0  * m_all(i+1, iRun);
    end
  end
end

semilogy(iter,(dm6(:,6)))
plot(iter,dm6(:,6))
