close all
clear all

% Import data
load('sim_results/results_IntAndNoise.mat')


%% Plot gain error vs SNR and interferer boresight angle
angles = 10:1:85;

figure
imagesc(angles,10*log10(snr_vals),squeeze(mean(g_mag)).')
title('Gain magnitude error')
xlabel('Interferer boresight angle (deg)')
ylabel('SNR (dB)')
set(gca, 'YDir','normal')
colorbar
colormap(jet);

figure
imagesc(angles,10*log10(snr_vals),squeeze(mean(g_phase)).')
title('Gain phase error')
xlabel('Interferer boresight angle (deg)')
ylabel('SNR (dB)')
set(gca, 'YDir','normal')
colorbar
colormap(jet);