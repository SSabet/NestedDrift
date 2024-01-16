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
surf(b,a,sb(:,:,1)')
set(gca,'FontSize',16)
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Liquid Savings, Low Type')

subplot(1,2,2)
surf(b,a,sb(:,:,2)')
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
surf(b,a,sa(:,:,1)')
set(gca,'FontSize',16)
view([-70 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
title('Illiquid Savings, Low Type')

subplot(1,2,2)
surf(b,a,sa(:,:,2)')
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

% Saving policies, liquid asset, fixing a
figure(7);
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
figure(8);
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
figure(9);
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
