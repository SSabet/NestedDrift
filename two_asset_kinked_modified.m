clear all; close all; clc;

%% Parameters
%par = readstruct("parameters.json");
parameters;

cellfun(@(x) assignin('base', x, par.(x)), fieldnames(par));

% check the parameter values
assert(chi0 < 1,'chi0 large, not interesting!');
assert(ra*chi1 < 1-chi0, 'ra too large, illiquid wealth accumulates to infinity');

la_mat = reshape(la_mat, 2, 2);
%% Grids
b = linspace(bmin,bmax,I)';
db = (bmax-bmin)/(I-1);

a = linspace(amin,amax,J);
da = (amax-amin)/(J-1);

[aa, bb] = meshgrid(a,b);
[aaa, bbb, zzz] = meshgrid(a,b,z);

Bswitch = [
    speye(I*J)*la_mat(1,1), speye(I*J)*la_mat(1,2);
    speye(I*J)*la_mat(2,1), speye(I*J)*la_mat(2,2)];

%matrix of liquid returns
Rb = rb_pos.*(bbb>0) + rb_neg.*(bbb<0);
%raa = ra.*ones(1,J);


%if ra>>rb, impose tax on ra*a at high a, otherwise some households
%accumulate infinite illiquid wealth (not needed if ra is close to or less than rb)
tau = 10; raa = ra.*(1 - (1.33.*amax./a).^(1-tau)); 
% plot(a,raa.*a);
%matrix of illiquid returns
Ra(:,:,1) = ones(I,1)*raa;
Ra(:,:,2) = ones(I,1)*raa;

d_upper = -(Ra.*aaa + xi*w*zzz);
d_lower = (chi0-1)/chi1.*aaa;
%% Utility Function
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
c_B = zeros(I,J,Nz);
c_F = zeros(I,J,Nz);
c = zeros(I,J,Nz);
updiag = zeros(I*J,Nz);
lowdiag = zeros(I*J,Nz);
centdiag = zeros(I*J,Nz);
AAi = cell(Nz,1);
BBi = cell(Nz,1);

%% INITIAL GUESS
tic;
% v0 = (((1-xi)*w*zzz + ra.*aaa + rb_neg.*bbb).^(1-gamma))/(1-gamma)/rho;
v0 = (((1-xi)*w*zzz + Ra.*aaa + Rb.*bbb).^(1-gamma))/(1-gamma)/rho;
%v0 = (((1-xi)*w*zzz - (d_upper+two_asset_kinked_cost(d_upper,aaa, chi0, chi1))+ Rb.*bbb).^(1-gamma))/(1-gamma)/rho;
%v0 = (((1-xi)*w*z c_B(2:Izz - (max(d_upper,d_lower)+two_asset_kinked_cost(max(d_upper,d_lower),aaa, chi0, chi1))+ Rb.*bbb).^(1-gamma))/(1-gamma)/rho;

v = v0;

for n=1:maxit
    
    V = v;   

    % Prepare value derivatives for use in FOC calculations
    VbF(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db; % forward difference wrt b
    VbB(2:I,:,:)   = (V(2:I,:,:)-V(1:I-1,:,:))/db; % backward difference wrt b
    VaF(:,1:J-1,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da; % forward difference wrt a
    VaB(:,2:J,:)   = (V(:,2:J,:)-V(:,1:J-1,:))/da; % backward difference wrt a

    % Build consumption and deposit policies conditional on assumed 
    %   liquid and illiquid drift directions 
    c_B(2:I,:,:)   = MU_inv(VbB(2:I,:,:));
    c_F(1:I-1,:,:) = MU_inv(VbF(1:I-1,:,:));
    
    d_BB = two_asset_kinked_FOC(VaB, VbB, a, par.chi0, par.chi1);
    d_BF = two_asset_kinked_FOC(VaF, VbB, a, par.chi0, par.chi1);
    d_FB = two_asset_kinked_FOC(VaB, VbF, a, par.chi0, par.chi1);
    d_FF = two_asset_kinked_FOC(VaF, VbF, a, par.chi0, par.chi1);

    % Impose the zero drift deposit in cases where the illiquid asset will overrun its state bounds otherwise
    d_FF(:,J,:) = d_upper(:,J,:);
    d_BF(:,J,:) = d_upper(:,J,:);
    d_FB(:,1,:) = d_upper(:,1,:);
    d_BB(:,1,:) = d_upper(:,1,:); 

    % Idenfity whether the conditional deposit policy creates drift consistent with the conditions used
    I_BB = (d_BB < d_upper);
    I_BF = (d_BF > d_upper);
    I_FF = (d_FF > d_upper);
    I_FB = (d_FB < d_upper);

    % Build deposit policies for forward and backward liquid drift cases, using only consistent deposit policies
    d_B = d_BF.*I_BF + d_BB.*I_BB + d_upper.*(~I_BB .* ~I_BF);
    d_B(1,:,:) = 0; % This case will never be used, just over-writing a nan for stability
    
    d_F = d_FF.*I_FF + d_FB.*I_FB + d_upper.*(~I_FB .* ~I_FF);
    d_F(I,:,:) = 0; % This policy will never be used, just over-writing a nan for stability
    
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
    d_0 = bdotzero(I_0,VaF,VaB,a,b,z,Rb,Ra,d_upper,d_lower,par);
    c_0 = (1-xi)*w*zzz + Rb.*bbb - d_0 - two_asset_kinked_cost(d_0,aaa, chi0, chi1);
    
    % Build unconditional policies
    c  = c_F.*I_F + c_B.*I_B + c_0.*I_0;
    d  = d_F.*I_F + d_B.*I_B + d_0.*I_0;
    sb = (1-xi)*w*zzz + Rb.*bbb - c - d - two_asset_kinked_cost(d,aaa, chi0, chi1);
    sa = Ra.*aaa + xi .* w .* zzz + d;
        
    if (sum(sum(sum(I_0<0)))>0)
        fprintf('There are elements with both I_F & I_B positive. V not concave in b? \n')
        [iii, jjj, kkk] =  ind2sub(size(I_0), find(I_0 < 0));
        not_ccv_indexes = [find(I_0 < 0), iii, jjj, kkk];
        
        check_idx = 5102:5105;
        [d_B(check_idx); d_F(check_idx)];
        [d_BB(check_idx); d_BF(check_idx); d_FB(check_idx); d_FF(check_idx)];
    end

    if (sum(sum(sum(VbB<0)))>0)
        fprintf('There are elements with VbB < 0. V not concave in b? \n') 
        [iii, jjj, kkk] =  ind2sub(size(VbB), find(VbB < 0));
         neg_Vb_indexes = [find(VbB < 0), iii, jjj, kkk];
    end
    
    u  = c.^(1-gamma)/(1-gamma);

    % Build transition matrix matrix
    A  = driftMatrixLiquid(sb,db,db,par) + driftMatrixIlliquid(sa,da,da,par) + Bswitch;
    
    if max(abs(sum(A,2)))>10^(-12)
        
        disp('Improper Transition Matrix')
        [ii, jj, kk] =  ind2sub(size(V), find(abs(sum(A,2))>10^(-12)));
        Improper_entries = [find(abs(sum(A,2)) > 10^(-12)), ii, jj, kk];
        find(abs(sum(A,2))>10^(-12));
        
        BB = driftMatrixLiquid(sb,db,db,par);
        AA = driftMatrixIlliquid(sa,da,da,par);
        find(abs(sum(BB,2))>10^(-12));
        find(abs(sum(AA,2))>10^(-12));

        break
    end
    
    B = (1/Delta + rho)*speye(I*J*Nz) - A;
    
    u_stacked = reshape(u,I*J*Nz,1);
    V_stacked = reshape(V,I*J*Nz,1);
    
    vec = u_stacked + V_stacked/Delta;
    
    V_stacked = B\vec; %SOLVE SYSTEM OF EQUATIONS
        
    V = reshape(V_stacked,I,J,Nz);   
    
    
    Vchange = V - v;
    v = V;
    
    b_dist(:,n) = max(abs(Vchange(:,:,2)),[],2);
    a_dist(:,n) = max(max(abs(Vchange),[],3),[],1);
    ab_dist(:,:,n) = max(abs(Vchange),[],3);
   
    dist(n) = max(max(max(abs(Vchange))));
    disp(['Value Function, Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
    
end
toc

% d = I_B.*d_B + Id_F.*d_F;
m = d + xi*w*zzz + Ra.*aaa;
s = (1-xi)*w*zzz + Rb.*bbb - d - two_asset_kinked_cost(d,aaa, chi0, chi1) - c;

sc = (1-xi)*w*zzz + Rb.*bbb - c;
sd = - d - two_asset_kinked_cost(d,aaa, chi0, chi1);

subplot(2,2,1)
plot(a,d(:,:,1),a,zeros(J,1),'k--')
xlabel('Illiquid Wealth, a')

subplot(2,2,2)
plot(a,m(:,:,1),a,zeros(J,1),'k--')
xlabel('Illiquid Wealth, a')

subplot(2,2,3)
plot(b,c(:,:,1))
xlabel('Liquid Wealth, b')

subplot(2,2,4)
plot(b,s(:,:,1),b,zeros(I,1),'k--')
xlabel('Liquid Wealth, b')

% Saving policies, liquid asset, fixing a
figure;
subplot(2,2,1)
j = 1;
plot(b, squeeze(sb(:,j,:)), 'LineWidth', 3)
title("Saving policy for liquid asset, fixing a ="+a(j))
yline(0, '--r', 'zero');
xlabel("b")
legend('Low-income','High-income')

subplot(2,2,2)
j = floor(J/4);
plot(b, squeeze(sb(:,j,:)), 'LineWidth', 3)
title("Saving policy for liquid asset, fixing a ="+a(j))
yline(0, '--r', 'zero');
xlabel("b")
legend('Low-income','High-income')

subplot(2,2,3)
j = floor(3*J/4);
plot(b, squeeze(sb(:,j,:)), 'LineWidth', 3)
title("Saving policy for liquid asset, fixing a ="+a(j))
yline(0, '--r', 'zero');
xlabel("b")
legend('Low-income','High-income')

subplot(2,2,4)
j = J;
plot(b, squeeze(sb(:,j,:)), 'LineWidth', 3)
title("Saving policy for liquid asset, fixing a ="+a(j))
yline(0, '--r', 'zero');
xlabel("b")
legend('Low-income','High-income')

% Saving policies, illiquid asset, fixing b
figure;
subplot(2,2,1)
i = 1;
plot(a, squeeze(sa(i,:,:)), 'LineWidth', 3)
title("Saving policy for illiquid asset, fixing b ="+b(i))
yline(0, '--r', 'zero');
xlabel("a")
legend('Low-income','High-income')

subplot(2,2,2)
plot(a, squeeze(sa(i,:,:)), 'LineWidth', 3)
i = floor(I/4);
title("Saving policy for illiquid asset, fixing b ="+b(i))
yline(0, '--r', 'zero');
xlabel("a")
legend('Low-income','High-income')

subplot(2,2,3)
plot(a, squeeze(sa(i,:,:)), 'LineWidth', 3)
i = floor(3*I/4);
title("Saving policy for illiquid asset, fixing b ="+b(i))
yline(0, '--r', 'zero');
xlabel("a")
legend('Low-income','High-income')

subplot(2,2,4)
plot(a, squeeze(sa(i,:,:)), 'LineWidth', 3)
i = I;
title("Saving policy for illiquid asset, fixing b ="+b(i))
yline(0, '--r', 'zero');
xlabel("a")
legend('Low-income','High-income')


% Deposit policies across liquid asset, fixing a
figure;
subplot(2,2,1)
j = 1;
plot(b, squeeze(d(:,j,:)), 'LineWidth', 3)
title("Deposit policy for liquid asset, fixing a ="+a(j))
yline(0, '--r', 'zero');
xlabel("b")
legend('Low-income','High-income')

subplot(2,2,2)
j = floor(J/4);
plot(b, squeeze(d(:,j,:)), 'LineWidth', 3)
title("Deposit policy for liquid asset, fixing a ="+a(j))
yline(0, '--r', 'zero');
xlabel("b")
legend('Low-income','High-income')


subplot(2,2,3)
j = floor(3*J/4);
plot(b, squeeze(d(:,j,:)), 'LineWidth', 3)
title("Deposit policy for liquid asset, fixing a ="+a(j))
yline(0, '--r', 'zero');
xlabel("b")
legend('Low-income','High-income')


subplot(2,2,4)
j = J;
plot(b, squeeze(d(:,j,:)), 'LineWidth', 3)
title("Deposit policy for liquid asset, fixing a ="+a(j))
yline(0, '--r', 'zero');
xlabel("b")
legend('Low-income','High-income')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATIONARY DISTRIBUTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

M  = I*J*Nz;
AT = A';

% Fix one value so matrix isn't singular:
vec        = zeros(M,1);
iFix       = 1657;
vec(iFix)  = 0.01;
AT(iFix,:) = [zeros(1,iFix-1),1,zeros(1,M-iFix)];

% Solve system:
g_stacked  = AT\vec;
g_sum      = g_stacked'*ones(M,1)*da*db;
g_stacked  = g_stacked./g_sum;

% Reshape results
g   = reshape(g_stacked,I,J,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make plots of the stationary policies & distributions
stationary_figures;
