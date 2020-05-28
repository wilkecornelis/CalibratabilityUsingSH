clear
close all

%% Plot gain estimate statistics vs SNR

% % Import data
load('sim_results\results_NoiseOnlyT1');
g_rmse_T1 = g_rmse;
LsignalT1 = Lsignal;

load('sim_results\results_NoiseOnlyT2');
g_rmse_T2 = g_rmse;
LsignalT2 = Lsignal;

load('sim_results\results_NoiseOnlyT3');
g_rmse_T3 = g_rmse;
LsignalT3 = Lsignal;

% Plot
figure
semilogy(10*log10(snr_vals),g_rmse_T1(:,1),'k-.')
hold on
semilogy(10*log10(snr_vals),g_rmse_T2(:,1),'k--')
semilogy(10*log10(snr_vals),g_rmse_T3(:,1),'k-')
hold off

legend(['Nsamples = ' num2str(LsignalT1)], ['Nsamples = ' num2str(LsignalT2)], ['Nsamples = ' num2str(LsignalT3)],'Location','southwest');
title('RMS Gain error vs SNR')
xlabel('SNR (dB)')
ylabel('RMS gain error')
grid on
grid off
grid on