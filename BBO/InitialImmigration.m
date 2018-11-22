function InitialImmigration
t = 1 : 0.01 : 3*pi/2;
close all; SetPlotOptions
plot(t, sin(t));
xlabel('Fitness');
ylabel('Immigration Rate');
set(gca, 'xtick', []);
set(gca, 'ytick', [])
axis([1 3*pi/2 -1 1]);
set(gca, 'box', 'off');