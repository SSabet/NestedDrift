clearvars; close; clc;

tic;

% -------------------------------------------------------------------------
% Read parameters

par = parameters;
cellfun(@(x) assignin('base', x, par.(x)), fieldnames(par));

% Check the parameter values are valid
assert(chi0 < 1,'chi0 large, not interesting!');
assert(ra*chi1 < 1-chi0, 'ra too large, illiquid wealth accumulates to infinity');

% -------------------------------------------------------------------------
% Prepare endogenous state grids

grids = makegrids(par);
cellfun(@(x) assignin('base', x, grids.(x)), fieldnames(grids));

% -------------------------------------------------------------------------
% Pre-computations

% Grid of asset returns
Rb = rb_pos.*(bbb>0) + rb_neg.*(bbb<0);
Ra = repmat(ones(I,J)*ra,1,1,Nz);

% SOROUSH - why turn off the diminishing returns at the top end like this??

% raa = ones(1,J)*ra;
% % if ra>>rb, impose tax on ra*a at high a, otherwise some households
% % accumulate infinite illiquid wealth (not needed if ra is close to or less than rb)
% tau = 10; raa = ra.*(1 - (1.33.*amax./a).^(1-tau)); % plot(a,raa.*a);
% just for clarity; otherwise premultiplying aaa by raa works as well, using implicit expansion
 
% Two particular (arrays of) points on the d-domain
d_zerodrift = -(Ra .* aaa + xi * w * zzz); % The d policy that ensures zero a drift
d_lower     = (chi0-1)/chi1.*aaa; % The d policy at which marginal adjustment cost equals withdrawal

% Poisson transition matrix between exogenous  states
Bswitch     = kron(la_mat, speye(I*J));

% -------------------------------------------------------------------------
% Solve stationary policies & value

% Initial guess
v0   = U((1-xi)*w*zzz + Ra.*aaa + Rb.*bbb,par)/rho; vNew = v0;

% Main loop
dist = ones(maxit,1);
for n=1:maxit
    
    % Relabel the current value guess
    vOld = vNew;   
    
    % Update policies & value
    [c,d,sb,sa,vNew,A] = updateHousehold(vOld,Rb,Ra,d_zerodrift,d_lower,Bswitch,grids,par);
    
    % Check convergence & iterate
    Vchange = vOld - vNew;
    
    dist(n) = max(max(max(abs(Vchange))));
    disp(['Value Function, Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end 
end
toc

% -------------------------------------------------------------------------
% Solve stationary distribution

g = stationary_dist(A, a, b, par);

% -------------------------------------------------------------------------
% Make plots of the stationary policies & distributions
stationary_figures;
