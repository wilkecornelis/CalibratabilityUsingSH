% This code tests self-holography on a LOFAR LBA dataset
% Calibration is applied to the supplied dataset and the SIR is calculated
% and displayed. The calibration results can be viewed by running
% figX_LOFARLBAresults.mat

% Nelis Wilke 11 September 2020
% Copyright 2020. See the LICENSE file at the top-level directory of this distribution.

close all
clear

%% Initialise constants
c = 2.9979245e8;                    % Speed of light in m/s
Nfreq = 512;                        % Number of subbands
freq = 43554688;                    % RCU 1 
lambda = c/freq;                    % Wavelegth in m
k = (2*pi)/lambda;

%% Import LBA config
% Three dataset selections available.
% select "1" for a dataset where CasA has altitude angle = 84 degrees(this is the dataset used for the results in the paper)
% Select "2" CasA has altitude angle = 46 degrees
% Select "3" CasA has altitude angle = 26 degrees (Extreme case)

load('data\LBA\LBA_config.mat')
data_sel = 1;

Sigma = Sigma(:,:,data_sel);
t_obs = t_obs(data_sel);

%% SH calibration
disp(['Dataset date (UTC): ', datestr(t_obs)]);
disp(['Dataset freq (Hz): ', int2str(freq)]);

% Use 3CR catalogoue
load('data\srclist3CR.mat');

% Select CasA 
srcsel = 324;
rasrc = srclist3CR(srcsel).alpha;
decsrc = srclist3CR(srcsel).delta;
epoch = true(length(srcsel), 1);

% Apparent lm coordinates of CasA 
[l_Cas,m_Cas] = radectolm(rasrc,decsrc,JulianDay(t_obs),stat_long,stat_lat);
% Calculate apparent altitude of CasA 
alt_Cas = acosd(sqrt(l_Cas^2 + m_Cas^2));
disp(['Calibration source altitude angle (degrees): ', num2str(alt_Cas)])
% Geometric delays of casA source
a_Cas = exp(-1i * k * r_ant(:,1:2) * [l_Cas,m_Cas].');

% Select CygnusA
srcsel = 283;
rasrc = srclist3CR(srcsel).alpha;
decsrc = srclist3CR(srcsel).delta;
epoch = true(length(srcsel), 1);
% Apparent lm coordinates of CygA source
[l_Cyg,m_Cyg] = radectolm(rasrc,decsrc,JulianDay(t_obs),stat_long,stat_lat);
% Calculate apparent altitude of CygA source
alt_Cyg = acosd(sqrt(l_Cyg^2 + m_Cyg^2));
% Geometric delays of CygA source
a_Cyg = exp(-1i * k * r_ant(:,1:2) * [l_Cyg,m_Cyg].');

% Set calibration source
cal_src = 'Cas';

if cal_src == 'Cas'
    l_c = l_Cas;
    m_c = m_Cas;
    a_c = a_Cas;
elseif cal_src == 'Cyg'
    l_c = l_Cyg;
    m_c = m_Cyg;
    a_c = a_Cyg;
end    

% Initialise calibration variables
Niter = 20;                                                   % Number of SH solving iterations
P = length(Sigma);                                            % Number of receive paths
sigma_cal = 1;                                                % Calibration source power
g_prev = ones(P,1);                                           % Initialise gains
g_iter = ones(Niter,P);

% Beamform Sigma to cal source
Sigma_bf = diag(a_c') *Sigma* diag(a_c);         % Apply beamforming weights
r2 = diag(Sigma_bf);                             % Measured autocorrelations

%% Start SH
for iter = 1:Niter
    
    % update reference beam weight vector
    w_ref = 1./ conj(g_prev);
    
    % Tile-reference beam crosscorrelations (r_xy)
    r1 = Sigma_bf*(w_ref);
    % Compile "y-vector"
    Y = [r1;r2];
    
    % Construct submatrices of the A-matrix using the current gain estimate
    a1 = ((g_prev)'*w_ref)*eye(P)*sigma_cal;
    a2 = diag(w_ref);
    a3 = diag(conj(g_prev)*sigma_cal);
    a4 = eye(P);
    
    % compile A-mat
    A = [a1 a2; a3 a4];
    
    % Solve unknown gains and system noise ("x-vector")
    x_est = A\Y;
    
    % update gain estimates
    g_est = x_est(1:P);
    % ensure proper phase referencing
    g_est = g_est / (g_est(1) / abs(g_est(1)));
    
    % apply averaging every second iteration
    if (mod(iter, 2) == 0)
        g_est = (g_est + g_prev) / 2;
        % ensure proper phase referencing
        g_est = g_est / (g_est(1) / abs(g_est(1)));
    end   
    % Store most recent gain estimates
    g_prev = g_est;
    
    % Store gain estimates for each iteration
    g_iter(iter,:) = g_est;
end

%% Calculate SIR

% Normalize estimated gains to unity magnitude
g_est = g_est./abs(g_est);
% Phase calibrate measured covariance matrix
Sigma_cal = diag(1 ./ (g_est)) * Sigma * diag(1 ./ g_est)';

% Estimate apparent calibration source power from measured acm
Rhat = Sigma_cal;
A = exp(-(2 * pi * 1i * freq / c) * (r_ant(:,1:2) * [l_c,m_c].'));             % Array response towards cal source
sigma_c = sigma_scale(data_sel)*real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * (Rhat(:)));    % Estimated flux
Sigma_c = sigma_c*(a_c*a_c');                                                           % Covariance matrix of cal source                              

% Demonstrate removal of cal source
% l = -1:0.005:1;
% m = -1:0.005:1;
% 
% skymapcal = acm2skyimage((Sigma_cal-Sigma_c), r_ant(:, 1), r_ant(:, 2), freq(freq_id), l, m);
% imagesc(l,m,(skymapcal.'));
% hold on
% plot(l_c,m_c,'x')
% colormap('jet')
% colorbar

% Beamform Sigma_c to cal source
Sigma_c_bf = diag(a_c')*Sigma_c*diag(a_c')';

% Measured array covariance matrix excluding noise (sky sources only)
Sigma_ci = Sigma_cal-diag(diag(Sigma_cal));         

% EVD analysis and Sigma_i reconstruct
[V,D] = eig(Sigma_ci);                            % Eigenvalue decomposition of Sigma exlcuding noise
eigv_sort = sort(real(diag(D)));            
eigv = diag(D) + abs(eigv_sort(1));               % Add abs of lowes eig value to all eig values. Force all eig values positive
Sigma_ci = V*diag(eigv)*V';                       % Reconstruct cov matrix for sky sources only
Sigma_i = Sigma_ci-Sigma_c;                       % Cov matrix of all interferers
Sigma_i_bf = diag(a_c')*Sigma_i*diag(a_c')';      % Interference covariance matrix beamformed to calibration source

% SIR calculation using (32) in journal paper
Rc = mean(abs((ones(P,1)'*Sigma_c_bf).' - diag(Sigma_c_bf)));
Ri = mean(abs((ones(P,1)'*Sigma_i_bf).' - diag(Sigma_i_bf)));

SIR = 10*log10(Rc/Ri);

disp(['SIR (dB): ',num2str(SIR)])

%% Save results for plotting
l = -1:0.005:1;
m = -1:0.005:1;
skymapcal = acm2skyimage(diag(1 ./ (g_est)) * Sigma * diag(1 ./ g_est)', r_ant(:, 1), r_ant(:, 2), freq, l, m);
skymapuncal = acm2skyimage( Sigma , r_ant(:, 1), r_ant(:, 2), freq, l, m);



save(['sim_results\LBA',num2str(data_sel),'.mat'],'skymapcal','skymapuncal','g_iter')






