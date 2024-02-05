% ==============================================================================
% Plot figures
% ==============================================================================

% ------------------------------------------------------------------------------
% Consumption policies
 
figure(1)
set(gcf,'PaperPosition',[0 0 16 9])
subplot(1,2,1)
surf(b,a,c(:,:,1)','EdgeAlpha',0.1)
set(gca,'FontSize',9)
xlabel('Liquid, b')
ylabel('Illiquid, a')
xlim([bmin bmax])
ylim([amin amax])
title('Consumption, Low Type')

subplot(1,2,2)
surf(b,a,c(:,:,2)','EdgeAlpha',0.1)
set(gca,'FontSize',9)
xlabel('Liquid, b')
ylabel('Illiquid, a')
xlim([bmin bmax])
ylim([amin amax])
title('Consumption, High Type')

print('plot_consumption','-dpng')

% ------------------------------------------------------------------------------
% Deposit policies

figure(2)
set(gcf,'PaperPosition',[0 0 16 9])
subplot(1,2,1)
surf(b,a,d(:,:,1)','EdgeAlpha',0.2)
set(gca,'FontSize',9)
view([-70 30])
xlabel('Liquid, b')
ylabel('Illiquid, a')
xlim([bmin bmax])
ylim([amin amax])
title('Deposits, Low Type')

subplot(1,2,2)
surf(b,a,d(:,:,2)','EdgeAlpha',0.2)
set(gca,'FontSize',9)
view([-70 30])
xlabel('Liquid, b')
ylabel('Illiquid, a')
xlim([bmin bmax])
ylim([amin amax])
title('Deposits, High Type')

print('plot_deposit','-dpng')

% ------------------------------------------------------------------------------
% Liquid saving policies

figure(3)
set(gcf,'PaperPosition',[0 0 16 9])
subplot(1,2,1)
surf(b,a,sb(:,:,1)','EdgeAlpha',0.2)
view([-40 20])
set(gca,'FontSize',9)
xlabel('Liquid, b')
ylabel('Illiquid, a')
xlim([bmin bmax])
ylim([amin amax])
title('Liquid Savings, Low Type')

subplot(1,2,2)
surf(b,a,sb(:,:,2)','EdgeAlpha',0.2)
view([-40 20])
set(gca,'FontSize',9)
xlabel('Liquid, b')
ylabel('Illiquid, a')
xlim([bmin bmax])
ylim([amin amax])
title('Liquid Savings, High Type')

print('plot_liquidSaving','-dpng')

% ------------------------------------------------------------------------------
% Illiquid saving policies

figure(4)
set(gcf,'PaperPosition',[0 0 16 9])
subplot(1,2,1)
surf(b,a,sa(:,:,1)','EdgeAlpha',0.2)
set(gca,'FontSize',9)
view([-40 20])
xlabel('Liquid, b')
ylabel('Illiquid, a')
title('Illiquid Savings, Low Type')

subplot(1,2,2)
surf(b,a,sa(:,:,2)','EdgeAlpha',0.2)
set(gca,'FontSize',9)
view([-40 20])
xlabel('Liquid, b')
ylabel('Illiquid, a')
xlim([bmin bmax])
ylim([amin amax])
title('Illiquid Savings, High Type')

print('plot_illiqduisSaving','-dpng')

% ------------------------------------------------------------------------------
% Distributions

figure(5)
set(gcf,'PaperPosition',[0 0 16 9])

subplot(1,2,1)
surf(b,a,g(:,:,1)','EdgeAlpha',0.1)
set(gca,'FontSize',9)
view([-70 30])
xlabel('Liquid, b')
ylabel('Illiquid, a')
title('Stationary Distribution, Low Type')

subplot(1,2,2)
surf(b,a,g(:,:,2)',EdgeAlpha=0.1)
set(gca,'FontSize',9)
view([-70 30])
xlabel('Liquid, b')
ylabel('Illiquid, a')
xlim([bmin bmax])
ylim([amin amax])
title('Stationary Distribution, High Type')

print('plot_distribution','-dpng')

% ------------------------------------------------------------------------------
% Phase diagrams overlaid with density contour

figure(10)
set(gcf,'PaperPosition',[0 0 16 9])
gsig = g .* (g>1e-10);
sc = 3; % Determines density of arrows on plot
subplot(1,2,1)
quiver(bbb(1:sc:end,1:sc:end,1), aaa(1:sc:end,1:sc:end,1), sb(1:sc:end,1:sc:end,1), sa(1:sc:end,1:sc:end,1),0)
hold on
contour(b,a,gsig(:,:,1)')
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Phase diagram, Low Type')
hold off

subplot(1,2,2)
quiver(bbb(1:sc:end,1:sc:end,2), aaa(1:sc:end,1:sc:end,2), sb(1:sc:end,1:sc:end,2), sa(1:sc:end,1:sc:end,2),0)
hold on
contour(b,a,gsig(:,:,2)')
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin bmax])
ylim([amin amax])
title('Phase diagram, High Type')
hold off

print('plot_phaseDiagrams','-dpng')
