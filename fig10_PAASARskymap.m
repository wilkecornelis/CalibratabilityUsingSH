clear
close all

% load results
load('sim_results/results_PAASARskymap.mat')

% Plot
l_range = -1:0.01:1;
m_range = -1:0.01:1;

figure
pointsize = 55;
scatter(-1*(lmn_nvsat(:,2)),-1*(lmn_nvsat(:,1)),'o');
hold on
scatter(-1*(lmn_calsat(:,2)),-1*(lmn_calsat(:,1)), pointsize, P_c,'filled');
text(-1*(lmn_calsat(:,2)-0.075),-1*(lmn_calsat(:,1)),num2str(round(SIR(1,1:length(lmn_calsat)).',1)),'Color','black');
xlabel('l')
ylabel('m')
xlim([-1 1])
ylim([-1 1])
daspect([1 1 1])
colormap(jet)
colorbar
grid on











