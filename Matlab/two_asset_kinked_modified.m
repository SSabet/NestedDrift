clearvars; close; clc;

%% Parameters
tic;
% read parameters
parameters;
% alternatively: read them from json (but json doesn't recognize ";", so
% reshape)
% par = readstruct("parameters.json");
% la_mat = reshape(la_mat, 2, 2);

cellfun(@(x) assignin('base', x, par.(x)), fieldnames(par));

% check the parameter values
assert(chi0 < 1,'chi0 large, not interesting!');
assert(ra*chi1 < 1-chi0, 'ra too large, illiquid wealth accumulates to infinity');

%% Grids
% grid types: you can choose from  
% 'Linear' = equispaced,  
% 'NL_symmetric' = nonlinear with concentration of points on both ends,
% 'NL_lefts' = nonlinear grid with concentration of points at the left end (lower bound),
% 'NL_rights' = nonlienar grid with concentration of points at the right end
bgrid_type = 'NL_symmetric';
bgrid_type = 'Linear';
agrid_type = 'NL_leftskewed';
grids = makegrids(bmin, bmax, I, bgrid_type, amin, amax, J, agrid_type, z);
cellfun(@(x) assignin('base', x, grids.(x)), fieldnames(grids));

da = (amax-amin)/(J-1);
db = (bmax-bmin)/(I-1);

% matrix of liquid returns
Rb = rb_pos.*(bbb>0) + rb_neg.*(bbb<0);

raa = ones(1,J)*ra;
% % if ra>>rb, impose tax on ra*a at high a, otherwise some households
% % accumulate infinite illiquid wealth (not needed if ra is close to or less than rb)
% tau = 10; raa = ra.*(1 - (1.33.*amax./a).^(1-tau)); % plot(a,raa.*a);

% matrix of illiquid returns
% just for clarity; otherwise premultiplying aaa by raa works as well, using implicit expansion
Ra = repmat(ones(I,1)*raa, 1,1,Nz); 

%% Utility function, marginal utility and its inverse
if gamma==1
    U = @(c) log(c);
else
    U = @(c) c.^(1-gamma)/(1-gamma);
end

MU = @(c) c.^(-gamma);
MU_inv = @(x) x.^(-1/gamma);

%% Preallocation
VbF = zeros(I,J,Nz);
VbB = zeros(I,J,Nz);
VaF = zeros(I,J,Nz);
VaB = zeros(I,J,Nz);
Va = zeros(I,J,Nz);
Vb = zeros(I,J,Nz);
c_B = zeros(I,J,Nz);
c_F = zeros(I,J,Nz);

dist = ones(maxit,1);
%% Precomputations, and Initial guess
% two particular (array of) points on the d-domain
d_zerodrift = -(Ra.*aaa + xi*w*zzz);
d_lower = (chi0-1)/chi1.*aaa;

% precomputing the Poisson transition matrix
Bswitch = kron(la_mat, speye(I*J));

% initial guess
v0 = U((1-xi)*w*zzz + Ra.*aaa + Rb.*bbb)/rho;

%% Main Loop
v = v0;

for n=1:maxit
    
    V = v;   

    % Prepare value derivatives for use in FOC calculations
    VbF(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))./(bbb(2:I,:,:)-bbb(1:I-1,:,:)); % forward difference wrt b
    %VbB(2:I,:,:)   = (V(2:I,:,:)-V(1:I-1,:,:))./(bbb(2:I,:,:)-bbb(1:I-1,:,:)); % backward difference wrt b
    VbB(2:I,:,:)   = VbF(1:I-1,:,:); % backward difference wrt b

    VaF(:,1:J-1,:) = (V(:,2:J,:)-V(:,1:J-1,:))./(aaa(:,2:J,:)-aaa(:,1:J-1,:)); % forward difference wrt a
    %VaB(:,2:J,:)   = (V(:,2:J,:)-V(:,1:J-1,:))./(aaa(:,2:J,:)-aaa(:,1:J-1,:)); % backward difference wrt a
    VaB(:,2:J,:)   = VaF(:,1:J-1,:); % backward difference wrt a
    
    % Build consumption and deposit policies conditional on assumed 
    %   liquid and illiquid drift directions 
    c_B(2:I,:,:)   = MU_inv(VbB(2:I,:,:));
    %c_F(1:I-1,:,:) = MU_inv(VbF(1:I-1,:,:));
    c_F(1:I-1,:,:) = c_B(2:I,:,:);

    % Vb = (V(2:I,:,:)-V(1:I-1,:,:))./(bbb(2:I,:,:)-bbb(1:I-1,:,:));
    % Va = (V(:,2:J,:)-V(:,1:J-1,:))./(aaa(:,2:J,:)-aaa(:,1:J-1,:));
    % c = MU_inv(Vb);

    d_BB = two_asset_kinked_FOC(VaB, VbB, a, par.chi0, par.chi1);
    d_BF = two_asset_kinked_FOC(VaF, VbB, a, par.chi0, par.chi1);
    d_FB = two_asset_kinked_FOC(VaB, VbF, a, par.chi0, par.chi1);
    d_FF = two_asset_kinked_FOC(VaF, VbF, a, par.chi0, par.chi1);

    % Impose the zero drift deposit in cases where the illiquid asset will overrun its state bounds otherwise
    d_FF(:,J,:) = d_zerodrift(:,J,:);
    d_BF(:,J,:) = d_zerodrift(:,J,:);
    d_FB(:,1,:) = d_zerodrift(:,1,:);
    d_BB(:,1,:) = d_zerodrift(:,1,:); 

    % Idenfity whether the conditional deposit policy creates drift consistent with the conditions used
    I_BB = (d_BB < d_zerodrift);
    I_BF = (d_BF > d_zerodrift);
    I_FF = (d_FF > d_zerodrift);
    I_FB = (d_FB < d_zerodrift);

    % Build deposit policies for forward and backward liquid drift cases, using only consistent deposit policies
    d_B = d_BF.*I_BF + d_BB.*I_BB + d_zerodrift.*(~I_BB .* ~I_BF);
    d_B(1,:,:) = 0; % This case will never be used, just over-writing a nan when updating d
    
    d_F = d_FF.*I_FF + d_FB.*I_FB + d_zerodrift.*(~I_FB .* ~I_FF);
    d_F(I,:,:) = 0; % This policy will never be used, just over-writing a nan when updating d
    
    % Build backward liquid drift policy, and an indicator for when it's consitent with itself
    sb_B = (1-xi)*w*zzz + Rb.*bbb - d_B -...
        two_asset_kinked_cost(d_B,aaa, chi0, chi1) - c_B;
    
    % at lower b-boundary don't use Vb_B; if Vb_F dosn't work, leave it to
    % the next part that deals with b-drift zero
    I_B = sb_B < 0; 
    I_B(1,:,:) = 0;

    % Build forward liquid drift policy, and an indicator for when it's consitent with itself
    sb_F = (1-xi)*w*zzz + Rb.*bbb - d_F -...
        two_asset_kinked_cost(d_F,aaa, chi0, chi1) - c_F;
    
    % ...and vice versa for bmax:
    I_F = (sb_F > 0) .* (I_B==0); % Give precedence to the backward drift if there's a clash
    I_F(I,:,:) = 0;
    
    % Find consumption and deposit policies for the case of zero liquid drift
    I_0 = 1 - I_B - I_F;
    d_0 = bdotzero(I_0,VaF,VaB,a,b,z,Rb,Ra,d_zerodrift,d_lower,MU, par);
    c_0 = (1-xi)*w*zzz + Rb.*bbb - d_0 - two_asset_kinked_cost(d_0,aaa, chi0, chi1);
    
    % Build unconditional policies
    c  = c_F.*I_F + c_B.*I_B + c_0.*I_0;
    d  = d_F.*I_F + d_B.*I_B + d_0.*I_0;
    sb = (1-xi)*w*zzz + Rb.*bbb - c - d - two_asset_kinked_cost(d,aaa, chi0, chi1);
    sa = Ra.*aaa + xi .* w .* zzz + d;
        
    % if (sum(sum(sum(I_0<0)))>0)
    %     fprintf('There are elements with both I_F & I_B positive. V not concave in b? \n')
    %     [iii, jjj, kkk] =  ind2sub(size(I_0), find(I_0 < 0));
    %     not_ccv_indexes = [find(I_0 < 0), iii, jjj, kkk];    % 
    % end
    % 
    % if (sum(sum(sum(VbB<0)))>0)
    %     fprintf('There are elements with VbB < 0. V not concave in b? \n') 
    %     [iii, jjj, kkk] =  ind2sub(size(VbB), find(VbB < 0));
    %      neg_Vb_indexes = [find(VbB < 0), iii, jjj, kkk];
    % end
    
    u  = U(c);

    % Build transition matrix matrix
    % A  = driftMatrixLiquid(sb,db,db,par) + driftMatrixIlliquid(sa,da,da,par) + Bswitch;
    A  = driftMatrixLiquid2(sb,bbb,par) + driftMatrixIlliquid2(sa,aaa,par) + Bswitch;
    A2 = driftMatrix(sb, bbb, sa,aaa, par) + Bswitch;

    % if max(abs(sum(A,2)))>10^(-12)
    % 
    %     disp('Improper Transition Matrix')
    %     [ii, jj, kk] =  ind2sub(size(V), find(abs(sum(A,2))>10^(-12)));
    %     Improper_entries = [find(abs(sum(A,2)) > 10^(-12)), ii, jj, kk];
    % 
    %     break
    % end
    
    B = (1/Delta + rho)*speye(I*J*Nz) - A;
    
    u_stacked = reshape(u,I*J*Nz,1);
    V_stacked = reshape(V,I*J*Nz,1);
    
    vec = u_stacked + V_stacked/Delta;
    
    V_stacked = B\vec; %SOLVE SYSTEM OF EQUATIONS
        
    V = reshape(V_stacked,I,J,Nz);   
    
    
    Vchange = V - v;
    v = V;
    
    dist(n) = max(max(max(abs(Vchange))));
    disp(['Value Function, Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
    
end
toc

%% %%%%%%%%%%%%%%%%%%%%%%%%
% STATIONARY DISTRIBUTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g=stationary_dist(A, da, db, par);
g = stationary_dist2(A, a, b, par);
%% %%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make plots of the stationary policies & distributions
% stationary_figures;
