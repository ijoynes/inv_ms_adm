function [stop,user_data] = callback(x,iter,state,user_data,opts)

global_vars

switch state
  case 'init'
    % fprintf('------------------------------------------------------------------------------------------------------------\n');
    % fprintf('|       Date/Time      | Iter |      f       |    |proj g|   |      r^2      |      r       |     m/m*     |\n');
    % fprintf('------------------------------------------------------------------------------------------------------------\n');
    fprintf('Start of emission source parameter inversion...\n\n');
    fprintf('-------------------------------------------------------------------------------------------------------\n');
    fprintf('|       Date/Time      | Iter |      f      |   |proj g|   |      r^2     |      r      |    m/m*     |\n');
    fprintf('-------------------------------------------------------------------------------------------------------\n');

    n = length(x);
    user_data.f = zeros(1, opts.maxits+1);
    user_data.g = zeros(n, opts.maxits+1);
    user_data.g_proj = zeros(n, opts.maxits+1);
    user_data.x = zeros(n, opts.maxits+1);
    user_data.its = 1;
    
  case 'iter'
    f = iter.f;
    g = iter.g;

    user_data.f(   iter.it + 1) = f;
    user_data.g(:, iter.it + 1) = g;
    user_data.x(:, iter.it + 1) = x;
    user_data.its = user_data.its + 1;

    s = volumetric_emission_rate(x);

    % Compute the projected gradient from boundary clipping.
    g_proj = projgr(x,g,lb,ub,nbd);

    % The observational component of the objective function and gradient
    % is computed by subtracting the regularization term.
    f_obs = f - 0.5*reg_par*((s-s_b)'*(B_inv*(s-s_b)));
    g_obs = g - reg_par*(B_inv*(s-s_b));

    r  = compute_correlation_coefficient(c, c_star);
    r2 = compute_coefficient_of_determination(c, c_star);

    m = dot(space_int_wgt, s);
    m_norm = m/m_star;

    fprintf('| %s | %4d | %8.5e |  %8.5e |  %9.5e | %8.5e | %8.5e |\n', datestr(now), iter.it, f, norm(g_proj,Inf),r2,r,m_norm);
    
    file_num = generate_file_num(iter.it, max_iters + 1);
    results_path = fullfile(iter_dir, [iter_label, file_num, '.mat']);
    save(results_path, 'f', 'f_obs', 'g', 'g_proj', 'g_obs', 's', 'x', 'c', 'r', 'r2', 'm', 'm_norm');

  case 'done'
    % fprintf('------------------------------------------------------------------------------------------------------------\n');
    fprintf('-------------------------------------------------------------------------------------------------------\n');
    fprintf('\n');
    
    % Save the initial values to the user_data array
    user_data.f(1) = f_0;
    user_data.g(:,1) = g_0;
    user_data.g_proj(:,1) = g_proj_0;
    user_data.x(:,1) = x_0;

    % Trim the user_data arrays if the optimization routine exits before
    % the maximum number of iteration is reached.
    user_data.f = user_data.f(1:user_data.its);
    user_data.g = user_data.g(:,1:user_data.its);
    user_data.g_proj = user_data.g_proj(:,1:user_data.its);
    user_data.x = user_data.x(:,1:user_data.its);
    user_data.iters = (0:user_data.its-1)';

end

stop = 0;
