function [cpol,dpol] = bdotzero(bdrift_zero,VaF,VaB,a,b,z,Rb,Ra, d_upper, d_lower, par)
    cellfun(@(x) assignin('caller', x, par.(x)), fieldnames(par));


    % Pre-allocation
    cpol   = zeros(I,J,Nz);
    dpol   = zeros(I,J,Nz);
    
    foc_resid = @(x,Va,dpositive) fone(b(bi),a(aj),z(zk),x,Va,dpositive,Rb(bi,aj,zk),par);
    
    N = length(VaF);
    
    % positive deposit
    I_dpos = foc_resid(0,VaF,1) < 0; % indicator
    lb_dpos = zeros(N,1); % lower bound
    I_updrift = 1;   % direction of FD for Va: use VaF for Va?
    

    % % example:
    % aroots = 1:4;
    % broots = -(5:8);
    % croots = zeros(4,1)
    % roots_coef = [aroots', broots', croots']
    % arrayfun(@(i) roots(roots_coef(i,:)), inds, 'UniformOutput', false)
    % cellfun(@max,arrayfun(@(i) roots(roots_coef(i,:)), inds, 'UniformOutput', false))

    dmax_quad_coefs = [chi1/(2*a), (chi0+1*ones(N,1)), -(Rb.*b + (1-xi)*w*z)];
    ub_dpos = cellfun(@max,arrayfun(@(i) roots(dmax_quad_coefs(i,:)), 1:N, 'UniformOutput', false))*(1-e-12);
    
    
    % case of deposit zero: inaction region (up-drift)
    I_dzero = (foc_resid(0,VaF(bi,aj,zk),0) <= 0).*(~I_dpos); % indicator for d = 0
    % d = 0;

    % negative deposit indicator
    I_dneg = (~I_dpos).*(~I_dzero);
    
    % negative desposit with dupper < dlower (case I)
    I_dneg_case1 = (~I_dpos).*(dupper < dlower);
    I_updrift = 1;

    lb_case1 = dlower;
    ub_case1 = 0;

    % negative deposit with updrift
    % Negative d but up-drift region
    I_dneg_case2up = (foc_resid(dupper,VaF,0) < 0);
    I_updrift = 1;
    
    lb_case1 = dupper;
    ub_case1 = 0;

    

end