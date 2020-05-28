close all 
clear

%% Import data

% Results for M=1
load('sim_results/results_IntOnly_M1.mat');
GerrM1_p = g_phase;     % Gain phase errors
GerrM1_m = g_mag;       % Gain magnitude errors

% Results for M=2
load('sim_results/results_IntOnly_M2.mat');
GerrM2_p = g_phase;     % Gain phase errors
GerrM2_m = g_mag;       % Gain magnitude errors

% Results for M=4
load('sim_results/results_IntOnly_M4.mat');
GerrM4_p = g_phase;     % Gain phase errors
GerrM4_m = g_mag;       % Gain magnitude errors


%% Plot statistics of gain estimates vs interferer boresight angle and tile size
int_angles = 10:1:85;       % interferer angles in deg


figure
plot(int_angles, mean(GerrM1_m), 'k-',...
     int_angles, mean(GerrM1_p), 'k--');
set(gca, 'FontSize', 16);
title('Gain error with 1-by-1 tiles');
xlabel('Interferer boresight angle');
ylabel('Gain error');
legend('gain magnitude', 'gain phase (rad)');
ylim([0 0.25])
grid on
 

figure
plot(int_angles, mean(GerrM2_m), 'k-',...
     int_angles, mean(GerrM2_p), 'k--');
set(gca, 'FontSize', 16);
title('Gain error with 2-by-2 tiles');
xlabel('Interferer boresight angle');
ylabel('Gain error');
legend('gain magnitude', 'gain phase (rad)');
ylim([0 0.25])
grid on
 

figure
plot(int_angles, mean(GerrM4_m), 'k-',...
     int_angles, mean(GerrM4_p), 'k--');
set(gca, 'FontSize', 16);
title('Gain error with 4-by-4 tiles');
xlabel('Interferer boresight angle');
ylabel('Gain error');
legend('gain magnitude', 'gain phase (rad)');
ylim([0 0.25])
grid on
