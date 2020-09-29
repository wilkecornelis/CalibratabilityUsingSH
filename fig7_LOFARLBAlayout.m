close all
clear

% Import data
load('data/LOFARLBAantpos.mat')


%% Plot LOFAR LBA antenna positions

plot(antpos(:,1),antpos(:,2),'*')
daspect([1 1 1])
title('LOFAR LBA antenna positions')
ylabel('y (m)')
xlabel('x (m)')
xlim([-50 50])
ylim([-50 50])
set(gca,'YTick',-50:10:50)
set(gca,'XTick',-50:10:50)
grid on
