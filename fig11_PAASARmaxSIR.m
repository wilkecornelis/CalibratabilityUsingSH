close all
clear

% Import data
load('sim_results/results_PAASAR.mat')

% Plot
figure
plot(SIR_max, 'k.-');
title('Maximum SIR at each date index')
ylabel('Maximum SIR (dB)')
xlabel('Date index')
grid on


