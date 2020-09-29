% This scripts plots the results from the LOFAR HBA application example
%
% Nelis Wilke 11 September 2020

clear all
close all

% Import data
load('sim_results\HBA.mat')

%% Plot LOFAR HBA postage stamp image calibrated
l = -0.15:0.0005:0.2;
m = 0.05:0.0005:0.25;

figure
P = imagesc(l,m,(skymapcal./(max(max(skymapcal)))));
title('LOFAR HBA calibrated skymap');
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
colormap(jet);
colorbar;


%% Plot LOFAR HBA postage stamp image uncalibrated
figure
P = imagesc(l,m,(skymapuncal./(max(max(skymapuncal)))));
title('LOFAR HBA uncalibrated skymap');
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
colormap(jet);
colorbar;

%% Plot gain estimate results

figure
plot(angle(g_iter))
title('Gain phases of all antennas vs iteration')
xlabel('Iteration')
ylabel('Gain phase (rad)')

figure
plot(abs(g_iter))
title('Gain magnitudes of all antennas vs iteration')
xlabel('Iteration')
ylabel('Gain magnitude')




