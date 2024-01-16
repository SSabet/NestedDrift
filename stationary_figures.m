%%%%%%%%%%%
% FIGURES %
%%%%%%%%%%%

 % Consumption policies
 
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

% Deposit policies
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

% Liquid saving policies
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

% Illiquid saving policies
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

% Distributions
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

% % Adjustment cost illustration
% dgrid = linspace(-0.05,0.05,100);
% figure(6)
% plot(dgrid,two_asset_kinked_cost(dgrid,ones(size(dgrid)), chi0, chi1),'Linewidth',2)
% set(gca,'FontSize',16)
% ylabel('Cost, % of Illiquid Stock')
% xlabel('Deposits/Withdrawal, % of Illiquid Stock')
