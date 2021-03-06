clear
close all

%% Plot gain estimate statistics vs SIR

% Import data
load('sim_results/results_IntOnly_M1.mat');
g_phaseM1 = g_phase;
RMSGerrM1 = g_rmse;
sortorderM1 = sortorder;
SIR_M1 = SIR;

load('sim_results/results_IntOnly_M2.mat');
RMSGerrM2 = g_rmse;
sortorderM2 = sortorder;
SIR_M2 = SIR;

load('sim_results/results_IntOnly_M4.mat');
g_phaseM4 = g_phase;
RMSGerrM4 = g_rmse;
sortorderM4 = sortorder;
SIR_M4 = SIR;

load('sim_results/results_IntOnly_rand.mat');
g_phase_rand = g_phase;
RMSGerr_rand = g_rmse;
sortorder_rand = sortorder;
SIR_rand = SIR;

% % Plot
figure
semilogy(10 * log10(SIR_M1), RMSGerrM1(sortorderM1), 'ko-');
hold on
semilogy(10 * log10(SIR_M2), RMSGerrM2(sortorderM2), 'kx-');
semilogy(10 * log10(SIR_M4), RMSGerrM4(sortorderM4), 'k-.');
semilogy(10 * log10(SIR_rand), RMSGerr_rand(sortorder_rand), 'rx');
hold off

title('RMS Gain error vs SIR')
xlabel('SIR (dB)');
ylabel('RMS gain error')
xlim([0,50])
legend('M = 1', 'M = 2', 'M = 4', 'Random');
grid on

% Plot
% figure
% semilogy(10 * log10(SIR_M1), RMSGerrM1(sortorderM1), 'ko-');
% hold on
% semilogy(10 * log10(SIR_M2), RMSGerrM2(sortorderM2), 'kx-');
% semilogy(10 * log10(SIR_M4), RMSGerrM4(sortorderM4), 'k-.');
% hold off
