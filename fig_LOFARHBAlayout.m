close all
clear

% Import data
load('data/LOFARHBAantpos.mat')

%% Plot LOFAR LBA antenna positions

% Define square representing a tile
x1=-2.575;
x2=2.575;
y1=-2.575;
y2=2.575;
x_0 = [x1, x2, x2, x1, x1];
y_0 = [y1, y1, y2, y2, y1];
square = [x_0.',y_0.'];

% Rotate square 
rot_theta = 75.635;
rot_mat = [cosd(rot_theta) -1*sind(rot_theta); sind(rot_theta) cosd(rot_theta)];

square_rot = rot_mat*square.';
square_rot = square_rot.';

% Plot tiles
for tile_id = 1:length(r_tile)
  tile = (square_rot+r_tile(tile_id,:));
  plot(tile(:,1),tile(:,2),'bl-','LineWidth', 2)
  hold on
end

daspect([1 1 1])
title('LOFAR HBA tile positions')
ylabel('y (m)')
xlabel('x (m)')
xlim([-35 35])
ylim([-35 35])
set(gca,'YTick',-50:10:50)
set(gca,'XTick',-50:10:50)
grid on
