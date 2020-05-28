clear
close all

% Import data
load('data/AF.mat')

%% Plot station pattern

% lm resolution
res = 500;
l_range = linspace(-1,1,res);
m_range = linspace(-1,1,res);

P = pcolor(l_range,m_range,20*log10(abs(AF)/max(max(abs(AF)))));
hold on
pointsize = 350;

% RGB for green
RGB = [0 255 0]/256;
scatter(0,0,pointsize,RGB,'*');

% RGB for black
RGB = [0 0 0]/256;
scatter(0.45,0,pointsize,RGB,'*');

set(P, 'EdgeColor', 'none');
caxis([-100 0]);
title('Station pattern with source positions indicated')
xlabel('l');
ylabel('m');
daspect([1 1 1]);
colorbar;
colormap(jet);
ax = gca;
ax.YDir = 'normal';


