clear all
close all

%% Plot LOFAR LBA skymap

% Import data
data_sel = 3;

load(['sim_results/LBA',num2str(data_sel),'.mat'])

%% Plot LOFAR LBA calibrated skymap
l = -1:0.005:1;
m = -1:0.005:1;

figure
P = pcolor(l,m,(skymapcal./(max(max(skymapcal)))));
title('LOFAR skymap calibrated');
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
set(gca,'YDir', 'normal', 'XDir', 'reverse');
set(P, 'EdgeColor', 'none');
colormap(jet);
colorbar;
daspect([1 1 1]);

%% Plot LOFAR LBA uncalibrated skymap
figure
P = pcolor(l,m,(skymapuncal./(max(max(skymapuncal)))));
title('LOFAR skymap uncalibrated');
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
set(gca,'YDir', 'normal', 'XDir', 'reverse');
set(P, 'EdgeColor', 'none');
colormap(jet);
colorbar;
daspect([1 1 1]);

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


