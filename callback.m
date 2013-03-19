function [stop,user_data] = test_callback(x,iter,state,user_data,opts)

global inputDir simDir c_k c_star m m_star

switch state
  case 'init'
    fprintf('------------------------------------------------------------------------------------------------------------\n');
    fprintf('|       Date/Time      | Iter |      f       |  inf_norm(g)  |      r^2      |      r       |     m/m*     |\n');
    fprintf('------------------------------------------------------------------------------------------------------------\n');

    user_data.f = zeros(1,opts.maxits);
    user_data.x = zeros(length(x),opts.maxits);
    user_data.its = 0;
    
  case 'iter'
    user_data.f(iter.it) = iter.f;
    user_data.x(:,iter.it) = x;
    user_data.its = user_data.its + 1;

    file_num = generate_file_num(nPass, nt);
    results_path = fullfile(iter_dir, [iter_label, file_num, '.mat']);

    s = volumetric_emission_rate(x);
    c = c_k;
    r  = compute_correlation_coefficient(c_k, c_star);
    r2 = compute_coefficient_of_determination(c_k, c_star);

    m = dot(spaceIntWeight, E_approx);
    m_norm = m/m_star;

    fprintf('| %s | %4d | %8.6e |  %8.6e |  %9.6e | %8.6e | %8.6e |\n', datestr(now), iter.it, iter.f, max(abs(iter.g)),r2,r,m_norm);
    s = x;
    f = iter.f;
    g = iter.g;
    save(results_path,'f','f_obs', 'df', 'df_obs', 's', 'x', 'c','r','r2','m','m_norm');

  case 'done'
    fprintf('------------------------------------------------------------------------------------------------------------\n');
    fprintf('\n');
    user_data.f = user_data.f(1:user_data.its);
    user_data.x = user_data.x(:,1:user_data.its);
end

stop = 0;
