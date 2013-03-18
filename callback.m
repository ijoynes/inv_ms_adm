function [stop,user_data] = callback(x,iter,state,user_data,opts)

global inputDir simDir c_k c_star m m_star f f_obs df df_obs s x

switch state
  case 'init'
    fprintf('------------------------------------------------------------------------------------------------------------\n');
    fprintf('|       Date/Time      | Iter |      f       |  inf_norm(g)  |      r^2      |      r       |     m/m*     |\n');
    fprintf('------------------------------------------------------------------------------------------------------------\n');
    nx = length(x);
    niters = opts.maxits + 1;
    
    user_data.m_star = m_star;
    user_data.its = 0;
    user_data.f = zeros(1, niters);
    user_data.f_obs = zeros(1, niters);
    user_data.r  = zeros(1, niters);
    user_data.r2 = zeros(1, niters);
    user_data.m = zeros(1, niters);

    user_data.x = zeros(nx, niters);
    user_data.s = zeros(nx, niters);
    user_data.g = zeros(nx, niters);
    user_data.g_obs = zeros(nx, niters);
    
    
  case 'iter'

    if iter.it == 1
        user_data.f(1) = f_0;
        user_data.f_obs(1) = f_obs_0;
        user_data.r(1) = r;
        user_data.r(1) = r2;
        user_data.m(1) = m;
        user_data.m_norm(1) = m/m_star;
        user_data.x
    end

    s = volumetric_emission_rate(x);
    c = c_k;
    r  = compute_correlation_coefficient(c_k, c_star);
    r2 = compute_coefficient_of_determination(c_k, c_star);
    m_norm = m/m_star;

    user_data.f(iter.it+1) = iter.f;
    user_data.f_obs(iter.it+1) = f_obs;
    user_data.r(iter.it+1) = r;
    user_data.r2(iter.it+1) = r2;


    user_data.x(:,iter.it+1) = x;
    user_data.its = user_data.its + 1;

    
    %m = dot(spaceIntWeight, E_approx);


    fprintf('| %s | %4d | %8.6e |  %8.6e |  %9.6e | %8.6e | %8.6e |\n', datestr(now), iter.it, iter.f, max(abs(iter.g)),r2,r,m_norm);
    s = x;
    f = iter.f;
    g = iter.g;

    file_num = generate_file_num(nPass, nt);
    results_path = fullfile(iter_dir, [iter_label, file_num, '.mat']);
    save(results_path,'f','f_obs', 'df', 'df_obs', 's', 'x' 'c','r','r2','m','m_norm');

  case 'done'
    fprintf('------------------------------------------------------------------------------------------------------------\n');
    fprintf('\n');
    user_data.f = user_data.f(1:user_data.its);
    user_data.x = user_data.x(:,1:user_data.its);
end

stop = 0;
