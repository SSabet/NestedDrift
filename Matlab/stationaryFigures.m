% ==============================================================================
% Plot figures
% ==============================================================================

% ------------------------------------------------------------------------------
% Consumption policies
 
figure()
set(gcf,'Units','centimeters','Position',[20 10 16 9]);

subplot(1,2,1)
surf(grids.b,grids.a,sol.c(:,:,1)','EdgeAlpha',0.1)
set(gca,'FontSize',9)
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('Low Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

subplot(1,2,2)
surf(grids.b,grids.a,sol.c(:,:,2)','EdgeAlpha',0.1)
set(gca,'FontSize',9)
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('High Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

exportgraphics(gcf,'plot_consumption.pdf','BackgroundColor','none')

% ------------------------------------------------------------------------------
% Deposit policies

figure()
set(gcf,'Units','centimeters','Position',[20 10 16 9]);

subplot(1,2,1)
surf(grids.b,grids.a,sol.d(:,:,1)','EdgeAlpha',0.2)
set(gca,'FontSize',9)
view([-70 30])
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('Deposits, Low Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

subplot(1,2,2)
surf(grids.b,grids.a,sol.d(:,:,2)','EdgeAlpha',0.2)
set(gca,'FontSize',9)
view([-70 30])
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('Deposits, High Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

exportgraphics(gcf,'plot_deposit.pdf','BackgroundColor','none')

% ------------------------------------------------------------------------------
% Liquid saving policies

figure()
set(gcf,'Units','centimeters','Position',[20 10 16 9]);

subplot(1,2,1)
surf(grids.b,grids.a,sol.sb(:,:,1)','EdgeAlpha',0.2)
view([-40 20])
set(gca,'FontSize',9)
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('Liquid Savings, Low Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex', 'TickLabelInterpreter', 'latex');

subplot(1,2,2)
surf(grids.b,grids.a,sol.sb(:,:,2)','EdgeAlpha',0.2)
view([-40 20])
set(gca,'FontSize',9)
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('Liquid Savings, High Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

exportgraphics(gcf,'plot_liquidSaving.pdf','BackgroundColor','none')

% ------------------------------------------------------------------------------
% Illiquid saving policies

figure()
set(gcf,'Units','centimeters','Position',[20 10 16 9]);

subplot(1,2,1)
surf(grids.b,grids.a,sol.sa(:,:,1)','EdgeAlpha',0.2)
set(gca,'FontSize',9)
view([-40 20])
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
title('Illiquid Savings, Low Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

subplot(1,2,2)
surf(grids.b,grids.a,sol.sa(:,:,2)','EdgeAlpha',0.2)
set(gca,'FontSize',9)
view([-40 20])
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('Illiquid Savings, High Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

exportgraphics(gcf,'plot_illiqduisSaving.pdf','BackgroundColor','none')

% ------------------------------------------------------------------------------
% Distributions

figure()
set(gcf,'Units','centimeters','Position',[20 10 16 9]);

subplot(1,2,1)
surf(grids.b,grids.a,sol.g(:,:,1)','EdgeAlpha',0.1)
set(gca,'FontSize',9)
view([-70 30])
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
zlim([0,inf])
title('Stationary Distribution, Low Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

subplot(1,2,2)
surf(grids.b,grids.a,sol.g(:,:,2)',EdgeAlpha=0.1)
set(gca,'FontSize',9)
view([-70 30])
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
zlim([0,inf])
title('Stationary Distribution, High Type','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex');

exportgraphics(gcf,'plot_distribution.pdf','BackgroundColor','none')

% ------------------------------------------------------------------------------
% Phase diagrams overlaid with density contour

figure()
set(gcf,'Units','centimeters','Position',[20 10 16 9]);

gsig = sol.g .* (sol.g>1e-10);
sc   = 3; % Determines density of arrows on plot

subplot(1,2,1)
quiver(grids.bbb(1:sc:end,1:sc:end,1), grids.aaa(1:sc:end,1:sc:end,1), sol.sb(1:sc:end,1:sc:end,1), sol.sa(1:sc:end,1:sc:end,1),0)
hold on
contour(grids.b,grids.a,gsig(:,:,1)')
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('Phase diagram, Low Type','Interpreter','latex')
hold off
set(gca, 'TickLabelInterpreter', 'latex');

subplot(1,2,2)
quiver(grids.bbb(1:sc:end,1:sc:end,2), grids.aaa(1:sc:end,1:sc:end,2), sol.sb(1:sc:end,1:sc:end,2), sol.sa(1:sc:end,1:sc:end,2),0)
hold on
contour(grids.b,grids.a,gsig(:,:,2)')
xlabel('b','Interpreter','latex')
ylabel('a','Interpreter','latex')
xlim([par.bmin par.bmax])
ylim([par.amin par.amax])
title('Phase diagram, High Type','Interpreter','latex')
hold off
set(gca, 'TickLabelInterpreter', 'latex');

exportgraphics(gcf,'plot_phaseDiagrams.pdf','BackgroundColor','none')
