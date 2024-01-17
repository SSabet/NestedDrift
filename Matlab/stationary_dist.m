function g = stationary_dist(A, da, db, par)
M          = par.I*par.J*par.Nz;
AT         = A';

% Fix one value so matrix isn't singular:
vec        = zeros(M,1);
iFix       = 10;
vec(iFix)  = 1;
AT(iFix,:) = [zeros(1,iFix-1),1,zeros(1,M-iFix)];

% Solve system:
g_stacked  = AT\vec;

% Normalise total distibution to 1
g_sum      = g_stacked'*ones(M,1)*da*db;
g_stacked  = g_stacked./g_sum;

% Reshape results
g          = reshape(g_stacked,par.I,par.J,par.Nz);
end