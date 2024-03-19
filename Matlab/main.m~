%==========================================================================
% This code solves the household block of a continuous time model with
% idiosyncratic risk, and two-assets with a convex adjustment cost on 
% transfers between the two, as in Kaplan, Moll & Violante (2018)
%
% Written by: Soroush Sabet & Patrick Schneider
% Date: February 2024
%==========================================================================

clearvars; close; clc;

tic;

% -------------------------------------------------------------------------
% Read parameters

par = parameters;
cellfun(@(x) assignin('base', x, par.(x)), fieldnames(par));

% Check the parameter values are valid
assert(rho - ra > 0, 'ra too large, more than compensates discount rate')
assert(rho - rb_pos > 0, 'rb_pos too large, more than compensates discount rate')
assert(chi0 < 1,'chi0 large, not interesting!');
assert(ra*chi1 < 1-chi0, "ra too large, illiquid wealth accumulates faster than the fastest withdrawal rate that's optimal under adjustment costs -> desired infinite accumulation");
assert(rb_neg * bmin + (1-xi)*w*z(1) > 0, "Liquid flow must be positive at all points in the state space! Consider raising w, z(low), or moving the borrowing constraint closer to zero");

% -------------------------------------------------------------------------
% Prepare endogenous state grids

grids = makegrids(par);
cellfun(@(x) assignin('base', x, grids.(x)), fieldnames(grids));

% -------------------------------------------------------------------------
% Pre-computations

% Grid of asset returns
Rb = rb_pos.*(bbb>0) + rb_neg.*(bbb<0);
Ra = repmat(ones(I,J)*ra,1,1,Nz);
 
% Two particular (arrays of) points on the d-domain
d_zerodrift = -(Ra .* aaa + xi * w * zzz); % The d policy that ensures zero a-drift
d_lower     = (chi0-1)/chi1.*aaa; % The d policy at which marginal adjustment cost equals withdrawal

% Poisson transition matrix between exogenous states
Bswitch     = kron(la_mat, speye(I*J));

% -------------------------------------------------------------------------
% Solve stationary policies & value
% Initial guess
v0   = U((1-xi)*w*zzz + Ra.*aaa + Rb.*bbb,par)/rho; 
vNew = v0; % The most recent guess
v0_alt = v0;

% Main loop
dist = ones(maxit,1); slowCount = 0;
for n=1:maxit
    
    try
        % Relabel the current value guess
        vOld = vNew;
        
        % Update policies & value
        [c,d,sb,sa,vNew,A] = updateHousehold(vOld,Rb,Ra,d_zerodrift,d_lower,Bswitch,grids,par);
        
        % Check convergence & iterate
        vChange = vOld - vNew;
        
        dist(n) = max(max(max(abs(vChange))));
        disp(['Value Function, Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
        if dist(n)<crit
            disp('Value Function Converged, Iteration = ')
            disp(n)
            break
        end 
    catch
        fprintf('\nTripped - looking for better guess\n')
        vNew = v0_alt; %Reset the starting guess
        par.Delta = .5; %Slow things down
        slowCount = 1;
    end

    if slowCount >0
        slowCount = slowCount + 1;
        if slowCount > 10
            fprintf('\nGoing fast again\n')
            v0_alt = vNew;
            slowCount = 0;
            par.Delta = 1e8;
        end
    end
end
toc

% -------------------------------------------------------------------------
% Solve stationary distribution

g  = stationaryDistribution(A, a, b, par);

% to check the stationary distribution/measure
g_next = (speye(Nz*I*J) - A'*(1000))\g(:);
max(abs(g_next- g(:)))

% -------------------------------------------------------------------------
% Make plots of the stationary policies & distributions
stationaryFigures;
