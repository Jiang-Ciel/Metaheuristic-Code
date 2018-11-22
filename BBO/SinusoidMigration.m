function SinusoidMigration
N = 50;
k = 1 : N;
mu = ( 1 - cos(pi*k / N )) / 2;
lambda = 1 - mu;
close all
figure; 
plot(k, mu, 'b-', 'LineWidth',2.0); hold on
plot(k, lambda, 'r--', 'LineWidth',2.0)
legend('emigration \mu', 'immigration \lambda')
box off
set(gca,'FontSize',14); set(gcf,'Color','White'); set(gca,'Box','on');
set(gca, 'Xtick', [])
xlabel('fitness')
ylabel('rate')