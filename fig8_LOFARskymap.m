clear
close all

% Import data
load('data/LOFARskymap.mat')


%% Plot LOFAR skymap

l = -1:0.01:1;
m = -1:0.01:1;
P = pcolor(l,m,(skymap./(max(max(skymap)))));
title('LOFAR skymap');
xlabel('l');
ylabel('m');
set(P, 'EdgeColor', 'none');
colormap(jet);
colorbar;
daspect([1 1 1]);


