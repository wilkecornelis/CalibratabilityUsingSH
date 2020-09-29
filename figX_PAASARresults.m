close all
clear

% Import data
load('sim_results/results_PAASAR.mat')

% The matrices containing the results are initialised for all satellites on all dates.
% However, since not all satellites are used for calibration, we first have to
% filter out the calibration results. The default entry is NaN.

% Filter results by finding indices of results that correspond to calibrations (!= NaN).
% Stack results in vectors.

SIR_cal_ind = find(isnan(SIR) == 0);

SIR_vec = reshape(SIR,[],1);
SIR_vec = SIR_vec(SIR_cal_ind);

g_rmse_vec = reshape(g_rmse,[],1);
g_rmse_vec = g_rmse_vec(SIR_cal_ind);

g_phase_vec = reshape(g_rmse,[],1);
g_phase_vec = g_phase_vec(SIR_cal_ind);

g_mag_vec = reshape(g_mag,[],1);
g_mag_vec = g_mag_vec(SIR_cal_ind);

%% plot

% Sort results according to SIR in ascending order
[SIR_sorted, sortorder] = sort(SIR_vec);
% 
% figure
semilogy((SIR_sorted),rad2deg(g_phase_vec(sortorder)), 'k.-');
title('RMS gain error as a function of SIR')
ylabel('RMS gain error')
xlabel('SIR (dB)')
grid on


