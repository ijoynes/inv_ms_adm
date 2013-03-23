function [stop,user_data] = callback(x,iter,state,user_data,opts)

global_vars

switch state
  case 'init'
    % fprintf('------------------------------------------------------------------------------------------------------------\n');
  % fprintf('|       Date/Time      | Iter |      f       |    |proj g|   |      r^2      |      r       |     m/m*     |\n');
    % fprintf('------------------------------------------------------------------------------------------------------------\n');
    fprintf('Start of emission source parameter inversion...\n\n');
    fprintf('-------------------------------------------------------------------------------------------------------\n');
    fprintf('|       Date/Time      | Iter |    fk/f0    | |pgk|/|pg0|  |      r^2     |      r      |    m/m*     |\n');
    fprintf('-------------------------------------------------------------------------------------------------------\n');
    n = length(x);
    nIters = opts.maxits+1;
    
    user_data.its = 1;
    user_data.m_star = m_star;
    user_data.s_star = s_star;
    user_data.s_b    = s_b;
    user_data.c_star = c_star;

    user_data.f      = nan(1, nIters);
    user_data.f_obs  = nan(1, nIters);
    user_data.r      = nan(1, nIters);
    user_data.r2     = nan(1, nIters);
    user_data.m      = nan(1, nIters);

    user_data.x      = nan(n, nIters);
    user_data.s      = nan(n, nIters);
    user_data.g      = nan(n, nIters);
    user_data.g_proj = nan(n, nIters);
    user_data.g_obs  = nan(n, nIters);
    user_data.c      = nan(nt, nReceptors, nIters);

  case 'iter'
    s = volumetric_emission_rate(x);
    f = iter.f;
    g = iter.g;

    % The observational component of the objective function and gradient
    % is computed by subtracting the regularization term.
    f_obs = f - 0.5*reg_par*((s-s_b)'*(B_inv*(s-s_b)));
    g_obs = g - reg_par*(B_inv*(s-s_b));

    % Compute the projected gradient from boundary clipping.
    g_proj = projgr(x,g,lb,ub,nbd);

    r  = compute_correlation_coefficient(c, c_star);
    r2 = compute_coefficient_of_determination(c, c_star);
    m = dot(space_int_wgt, s);
    m_norm = m/m_star;

    user_data.its = user_data.its + 1;
    user_data.f(   iter.it + 1) = f;
    user_data.f_obs(:, iter.it+1) = f_obs;
    user_data.g(:, iter.it + 1) = g;
    user_data.g_obs(:, iter.it+1) = g_obs;
    user_data.g_proj(:, iter.it+1) = g_proj;
    user_data.x(:, iter.it + 1) = x;
    user_data.m(:,iter.it+1) = m;
    user_data.r(iter.it+1) = r;
    user_data.r2(iter.it+1) = r2;
    user_data.x(:,iter.it+1) = x;
    user_data.s(:,iter.it+1) = s;
    user_data.c(:,:,iter.it+1) = c;

    fprintf('| %s | %4d | %8.5e |  %8.5e |  %9.5e | %8.5e | %8.5e |\n', datestr(now), iter.it, f/f_0, norm(g_proj,Inf)/norm(g_proj_0,Inf) ,r2,r,m_norm);
    
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
    user_data.f_obs(1) = f_obs_0;
    user_data.g_obs(:,1) = g_obs_0;
    user_data.r(1) = r_0;
    user_data.r2(1) = r2_0;
    user_data.m(1) = m_0;
    user_data.s(:,1) = s_0;
    user_data.c(:,:,1) = c_0;
    
    % Trim the user_data arrays if the optimization routine exits before
    % the maximum number of iteration is reached.
    n = user_data.its;
    user_data.iters = (0:n-1)';

    user_data.f      = user_data.f(1:n);
    user_data.f_obs  = user_data.f_obs(1:n);

    user_data.x      = user_data.x(:,1:n);
    user_data.s      = user_data.s(:,1:n);

    user_data.g      = user_data.g(:,1:n);
    user_data.g_obs  = user_data.g_obs (:,1:n);
    user_data.g_proj = user_data.g_proj(:,1:n);
    
    user_data.r      = user_data.r(1:n);
    user_data.r2     = user_data.r2(1:n);
    user_data.m      = user_data.m(1:n);

end

stop = 0;
