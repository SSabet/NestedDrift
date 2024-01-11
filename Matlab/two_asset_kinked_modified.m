clear all; close all; clc;

%% Parameters
par = readstruct("parameters.json");
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
raa = ra.*ones(1,J);

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
v0 = (((1-xi)*w*zzz + ra.*aaa + rb_neg.*bbb).^(1-gamma))/(1-gamma)/rho;
v0 = (((1-xi)*w*zzz - (d_upper+two_asset_kinked_cost(d_upper,aaa, chi0, chi1))...
    + Rb.*bbb).^(1-gamma))/(1-gamma)/rho;
v = v0;

for n=1:maxit
    fprintf('Iteration %d. \n', n);
    V = v;   
    %DERIVATIVES W.R.T. b
    VbF(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db; % forward difference
    %VbF(I,:,:) = ((1-xi)*w*zzz(I,:,:) + Rb(I,:,:).*bmax).^(-gamma); %state constraint boundary condition
    
    VbB(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db; % backward difference
    %VbB(1,:,:) = ((1-xi)*w*zzz(1,:,:) + Rb(1,:,:).*bmin).^(-gamma); %state constraint boundary condition

    %DERIVATIVES W.R.T. a
    VaF(:,1:J-1,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da;  % forward difference
    VaB(:,2:J,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da; % backward difference
    
    c_B(2:I,:,:) = MU_inv(VbB(2:I,:,:));
    c_F(1:I-1,:,:) = MU_inv(VbF(1:I-1,:,:));
    [~, d_BB] = update_policy(VaB, VbB, MU_inv, aaa, d_upper, par);
    [~,d_BF] = update_policy(VaF, VbB, MU_inv, a, d_upper, par);
    [~, d_FB] = update_policy(VaB, VbF, MU_inv, aaa, d_upper, par);
    [~,d_FF] = update_policy(VaF, VbF, MU_inv, a, d_upper, par);
    
    d_FF(:,J,:) = d_upper(:,J,:);
    d_BF(:,J,:) = d_upper(:,J,:);
    d_FB(:,1,:) = d_upper(:,1,:);
    d_BB(:,1,:) = d_upper(:,1,:); 

    d_F = d_FF.*+d_FB.*(d_FB < d_upper);
    d_B = d_BF.*(d_BF > d_upper)+d_BB.*(d_BB < d_upper);
    
    sb_B = (1-xi)*w*zzz + Rb.*bbb - d_B -...
        two_asset_kinked_cost(d_B,aaa, chi0, chi1) - c_B;
    % at lower b-boundary don't use Vb_B; if Vb_F dosn't work, leave it to
    % the next part that deals with b-drift zero
    sb_B(1,:,:) = 0; 
    I_B = sb_B < 0;

    sb_F = (1-xi)*w*zzz + Rb.*bbb - d_F -...
        two_asset_kinked_cost(d_F,aaa, chi0, chi1) - c_F;
    % vice versa for bmax:
    sb_F(I,:,:) = 0;
    I_F = sb_F > 0;
    
    I_0 = 1 - I_B - I_F;

    if (sum(sum(sum(I_0<0)))>0)
        fprintf('There are elements with both I_F & I_B positive. V not concave in b?')
    end

    % to Patrick: from this point on we should apply your procedure to any
    % points where I_0 is one, the c_0 vector will be the result of your
    % code for sb = 0

    %c_0 = (1-xi)*w*zzz + Rb.*bbb;
    [c_0, d_0] = bdotzero(I_0,VaF,VaB,a,b,z,Rb,par);
    c = c_F.*I_F + c_B.*I_B + c_0.*I_0;
    u = c.^(1-gamma)/(1-gamma);
    
    %CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
    X = -Ic_B.*sc_B/db -Id_B.*sd_B/db;
    Y = (Ic_B.*sc_B - Ic_F.*sc_F)/db + (Id_B.*sd_B - Id_F.*sd_F)/db;
    Z = Ic_F.*sc_F/db + Id_F.*sd_F/db;
    
    for i = 1:Nz
        centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
    end

    lowdiag(1:I-1,:) = X(2:I,1,:);
    updiag(2:I,:) = Z(1:I-1,1,:);
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
        updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end
    
    for nz=1:Nz
        BBi{nz}=spdiags(centdiag(:,nz),0,I*J,I*J)+spdiags([updiag(:,nz);0],1,I*J,I*J)+spdiags([lowdiag(:,nz);0],-1,I*J,I*J);
    end
    
    BB = [BBi{1}, sparse(I*J,I*J); sparse(I*J,I*J), BBi{2}];
    
    %CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
    % to Patrick: we shouldn't need the boundary lines below
    dB = Id_B.*dBB + Id_F.*dFB;
    dF = Id_B.*dBF + Id_F.*dFF;
    MB = min(dB,0);
    MF = max(dF,0) + xi*w*zzz + Ra.*aaa;
    MB(:,J,:) = xi*w*zzz(:,J,:) + dB(:,J,:) + Ra(:,J,:).*amax; %this is hopefully negative
    MF(:,J,:) = 0;
    chi = -MB/da;
    yy =  (MB - MF)/da;
    zeta = MF/da;
    
    %MATRIX AAi
    for nz=1:Nz
        %This will be the upperdiagonal of the matrix AAi
        AAupdiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
        for j=1:J
            AAupdiag=[AAupdiag;zeta(:,j,nz)];
        end
        
        %This will be the center diagonal of the matrix AAi
        AAcentdiag= yy(:,1,nz);
        for j=2:J-1
            AAcentdiag=[AAcentdiag;yy(:,j,nz)];
        end
        AAcentdiag=[AAcentdiag;yy(:,J,nz)];
        
        %This will be the lower diagonal of the matrix AAi
        AAlowdiag=chi(:,2,nz);
        for j=3:J
            AAlowdiag=[AAlowdiag;chi(:,j,nz)];
        end
        
        %Add up the upper, center, and lower diagonal into a sparse matrix
        AAi{nz} = spdiags(AAcentdiag,0,I*J,I*J)+spdiags(AAlowdiag,-I,I*J,I*J)+spdiags(AAupdiag,I,I*J,I*J);
        
    end
    
    AA = [AAi{1}, sparse(I*J,I*J); sparse(I*J,I*J), AAi{2}];
    
    A = AA + BB + Bswitch;
    
    if max(abs(sum(A,2)))>10^(-12)
        disp('Improper Transition Matrix')
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

d = Id_B.*d_B + Id_F.*d_F;
m = d + xi*w*zzz + Ra.*aaa;
s = (1-xi)*w*zzz + Rb.*bbb - d - two_asset_kinked_cost(d,aaa, chi0, chi1) - c;

sc = (1-xi)*w*zzz + Rb.*bbb - c;
sd = - d - two_asset_kinked_cost(d,aaa, chi0, chi1);

subplot(1,2,1)
plot(a,d(:,:,1),a,zeros(J,1),'k--')
xlabel('Illiquid Wealth, a')

subplot(1,2,2)
plot(a,m(:,:,1),a,zeros(J,1),'k--')
xlabel('Illiquid Wealth, a')

subplot(1,2,1)
plot(b,c(:,:,1))
xlabel('Liquid Wealth, b')

subplot(1,2,2)
plot(b,s(:,:,1),b,zeros(I,1),'k--')
xlabel('Liquid Wealth, b')




%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATIONARY DISTRIBUTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECONSTRUCT TRANSITION MATRIX WITH SIMPLER UPWIND SCHEME
X = -min(s,0)/db;
Y = min(s,0)/db - max(s,0)/db;
Z = max(s,0)/db;

for i = 1:Nz
    centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
end

lowdiag(1:I-1,:) = X(2:I,1,:);
updiag(2:I,:) = Z(1:I-1,1,:);
for j = 2:J
    lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
    updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
end

for nz=1:Nz
    BBi{nz}=spdiags(centdiag(:,nz),0,I*J,I*J)+spdiags([updiag(:,nz);0],1,I*J,I*J)+spdiags([lowdiag(:,nz);0],-1,I*J,I*J);
end

BB = [BBi{1}, sparse(I*J,I*J); sparse(I*J,I*J), BBi{2}];

%CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
chi = -min(m,0)/da;
yy =  min(m,0)/da - max(m,0)/da;
zeta = max(m,0)/da;


%MATRIX AAi
for nz=1:Nz
    %This will be the upperdiagonal of the matrix AAi
    AAupdiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
    for j=1:J
        AAupdiag=[AAupdiag;zeta(:,j,nz)];
    end
    
    %This will be the center diagonal of the matrix AAi
    AAcentdiag= yy(:,1,nz);
    for j=2:J-1
        AAcentdiag=[AAcentdiag;yy(:,j,nz)];
    end
    AAcentdiag=[AAcentdiag;yy(:,J,nz)];
    
    %This will be the lower diagonal of the matrix AAi
    AAlowdiag=chi(:,2,nz);
    for j=3:J
        AAlowdiag=[AAlowdiag;chi(:,j,nz)];
    end
    
    %Add up the upper, center, and lower diagonal into a sparse matrix
    AAi{nz} = spdiags(AAcentdiag,0,I*J,I*J)+spdiags(AAlowdiag,-I,I*J,I*J)+spdiags(AAupdiag,I,I*J,I*J);
    
end

AA = [AAi{1}, sparse(I*J,I*J); sparse(I*J,I*J), AAi{2}];
A = AA + BB + Bswitch;

M = I*J*Nz;
AT = A';
% Fix one value so matrix isn't singular:
vec = zeros(M,1);
iFix = 1657;
vec(iFix)=.01;
AT(iFix,:) = [zeros(1,iFix-1),1,zeros(1,M-iFix)];

% Solve system:
g_stacked = AT\vec;
g_sum = g_stacked'*ones(M,1)*da*db;
g_stacked = g_stacked./g_sum;


g(:,:,1) = reshape(g_stacked(1:I*J),I,J);
g(:,:,2) = reshape(g_stacked(I*J+1:I*J*2),I,J);

%%%%%%%%%%%
% FIGURES %
%%%%%%%%%%%

figure(1)
set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(b,a,c(:,:,1)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Consumption, Low Type')

subplot(1,2,2)
surf(b,a,c(:,:,2)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Consumption, High Type')
%print -depsc consumption.eps

figure(2)
set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(b,a,d(:,:,1)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Deposits, Low Type')

subplot(1,2,2)
surf(b,a,d(:,:,2)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Deposits, High Type')
%print -depsc deposits.eps


figure(3)
set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(b,a,s(:,:,1)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Liquid Savings, Low Type')

subplot(1,2,2)
surf(b,a,s(:,:,2)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Liquid Savings, High Type')
%print -depsc sav.eps

figure(4)
set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(b,a,m(:,:,1)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
title('Illiquid Savings, Low Type')

subplot(1,2,2)
surf(b,a,m(:,:,2)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Illiquid Savings, High Type')
%print -depsc ill_sav.eps


figure(5)
icut = 25; jcut=20;
bcut=b(1:icut); acut=a(1:jcut);
gcut = g(1:icut,1:jcut,:);

set(gcf,'PaperPosition',[0 0 30 10])
subplot(1,2,1)
surf(bcut,acut,gcut(:,:,1)')
set(gca,'FontSize',16)
view([-10 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin b(icut)])
ylim([amin a(jcut)])
title('Stationary, Low Type')

subplot(1,2,2)
surf(bcut,acut,gcut(:,:,2)')
set(gca,'FontSize',16)
view([-10 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin b(icut)])
ylim([amin a(jcut)])
title('Stationary Distribution, High Type')
%print -depsc stat_dist.eps


dgrid = linspace(-0.05,0.05,100)
figure(6)
plot(dgrid,two_asset_kinked_cost(dgrid,ones(size(dgrid)), chi0, chi1),'Linewidth',2)
set(gca,'FontSize',16)
ylabel('Cost, % of Illiquid Stock')
xlabel('Deposits/Withdrawal, % of Illiquid Stock')
%print -depsc adj_cost.eps