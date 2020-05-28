clear
close all

%% Import data

% Results for M=1
load('sim_results/results_IntOnly_M1.mat');
rxyM1_i = rxy_i;    % Tile-reference beam crosscorrelations excl. calibration source
rxxM1_i = rxx_i;    % Receive path autocorrelations excl. calibration source

% Results for M=2
load('sim_results/results_IntOnly_M2.mat');
rxyM2_i = rxy_i;
rxxM2_i = rxx_i;

% Results for M=4
load('sim_results/results_IntOnly_M4.mat');
rxyM4_i = rxy_i;
rxxM4_i = rxx_i;

%% Plot correlation statistics
int_angles = 10:1:85;

figure
plot(int_angles, mean(abs(rxyM1_i)), 'k-');
set(gca, 'FontSize', 16);
title('Average correlation magnitude with 1-by-1 tiles');
ylabel('mean |r_{xy}^{int}| ')
hold on
yyaxis right
plot(int_angles, mean(rxxM1_i), 'k--');
axis([10 90 0 max(mean(rxxM1_i))]);
set(gca, 'YColor', [0 0 0]);
ylabel('mean |r_{xx}^{int}| ')
xlabel('Interferer boresight angle')
legend('r_{xy}', 'r_{xx}');
yyaxis left
grid on
hold off

figure
plot(int_angles, mean(abs(rxyM2_i)), 'k-');
set(gca, 'FontSize', 16);
title('Average correlation magnitude with 2-by-2 tiles');
ylabel('mean |r_{xy}^{int}| ')
hold on
yyaxis right
plot(int_angles, mean(rxxM2_i), 'k--');
axis([10 90 0 max(mean(rxxM2_i))]);
set(gca, 'YColor', [0 0 0]);
ylabel('mean |r_{xx}^{int}| ')
xlabel('Interferer boresight angle')
legend('r_{xy}', 'r_{xx}');
yyaxis left
grid on
hold off

figure
plot(int_angles, mean(abs(rxyM4_i)), 'k-');
set(gca, 'FontSize', 16);
title('Average correlation magnitude with 4-by-4 tiles');
ylabel('mean |r_{xy}^{int}|')
hold on
yyaxis right
plot(int_angles, mean(rxxM4_i), 'k--');
axis([10 90 0 max(mean(rxxM4_i))]);
set(gca, 'YColor', [0 0 0]);
ylabel('mean |r_{xx}^{int}|')
xlabel('Interferer boresight angle')
legend('r_{xy}', 'r_{xx}');
yyaxis left
grid on
hold off
