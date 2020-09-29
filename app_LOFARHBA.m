% This code tests self-holography on a LOFAR HBA dataset
% Calibration is applied to the supplied dataset and the SIR is calculated
% and displayed. The calibration results can be viewed by running
% figX_LOFARHBAresults.mat
%
% Nelis Wilke 11 September 2020
% Copyright 2020. See the COPYRIGHT file at the top-level directory of this distribution.

close all
clear

%% Initialise constants
c = 2.9979245e8;                     % Speed of light in m/s     
Nfreq = 512;                         % Number of subbands
freq = 1e8:1e8/Nfreq:2e8-1e-3;       % RCU 5     
freq_id = 450;                       % Subband ID
freq(freq_id);                       % Subband select
lambda = c/freq(freq_id);            % Wavelegth in m
k = (2*pi)/lambda;                  


%% Import HBA data
load('data\HBA\HBA_config.mat')
%% SH calibration

% Import source position
disp(['Dataset date (UTC): ', datestr(t_obs(freq_id))]);
disp(['Dataset freq (Hz): ', int2str(freq(freq_id))]);

% Use 3CR catalogoue
load('data\srclist3CR.mat');         

% Select cal source
srcsel = 283;   % CygA                        
rasrc = srclist3CR(srcsel).alpha;       
decsrc = srclist3CR(srcsel).delta; 
epoch = true(length(srcsel), 1);

% Apparent lm coordinates of cal source
[l_c,m_c] = radectolm(rasrc,decsrc,JulianDay(t_obs(freq_id)),stat_long,stat_lat);

% Calculate apparent altitude of cal source
alt_c = acosd(sqrt(l_c^2 + m_c^2));
disp(['Calibration source altitude angle (degrees): ', num2str(alt_c)])

% Geometric delays of cal source
a_c = exp(-1i * k * r_tile(:,1:2) * [l_c,m_c].');

% Initialise calibration variables
Niter = 20;                                                   % Number of SH solving iterations
P = length(Sigma);                                            % Number of receive paths
sigma_cal = 1;                                                % Calibration source power
g_prev = ones(P,1);                                           % Initialise gains
g_iter = ones(Niter,P);

Sigma_bf = diag(a_c') * Sigma * diag(a_c);         % Apply beamforming weights

% Measured autocorrelations
r2 = diag(Sigma_bf);    

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

%% Estimate SIR

% Normalize estimated gains to unity magnitude
g_est = g_est./abs(g_est);
% Phase calibrate measured covariance matrix
Sigma_cal = diag(1 ./ (g_est)) * Sigma * diag(1 ./ g_est)';

% Estimate apparent calibration source power from measured acm
Rhat = Sigma_cal;
A = exp(-(2 * pi * 1i * freq(freq_id) / c) * (r_tile(:,1:2) * [l_c,m_c].'));
sigma_c = 0.8*real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * (Rhat(:)));
% Calculate model covariance matrix for calibration source
Sigma_c = sigma_c*(a_c*a_c');

% Demonstrate removal of cal source
% l = -0.15:0.0005:0.2;
% m = 0.05:0.0005:0.25;
% skymapcal = acm2skyimage((Sigma_cal-Sigma_c), r_tile(:,1), r_tile(:,2), freq(freq_id), l, m);
% figure
% imagesc(l,m,(skymapcal.'));
% hold on
% plot(l_c,m_c,'x')
% colormap('jet')
% colorbar

% Model array covariance matrix for cal source
Sigma_c_bf = diag(a_c')*Sigma_c*diag(a_c')';

% Measured array covariance matrix excluding noise (sky sources only)
Sigma_ci = Sigma_cal-diag(diag(Sigma_cal));         

% EVD analysis and Sigma_i reconstruct
[V,D] = eig(Sigma_ci);                   % Eigenvalue decomposition of Sigma exlcuding noise
eigv_sort = sort(real(diag(D)));            
eigv = diag(D) + abs(eigv_sort(1));      % Add abs of lowes eig value to all eig values. Force all eig values positive
Sigma_ci = V*diag(eigv)*V';              % Reconstruct cov matrix for sky sources only
Sigma_i = Sigma_ci-Sigma_c;              % Cov matrix of all interferers
Sigma_i_bf = diag(a_c')*Sigma_i*diag(a_c')';

% SIR calculation using (33) in journal paper
Rc = mean(abs((ones(P,1)'*Sigma_c_bf).' - diag(Sigma_c_bf)));
Ri = mean(abs((ones(P,1)'*Sigma_i_bf).' - diag(Sigma_i_bf)));

SIR = 10*log10(Rc/Ri);
disp(['SIR (dB): ', num2str(SIR)])

%% Save results for plotting
l = -0.15:0.005:0.2;
m = 0.05:0.005:0.25;
skymapcal = acm2skyimage(diag(1 ./ (g_est)) * Sigma * diag(1 ./ g_est)', r_tile(:, 1), r_tile(:, 2), freq(freq_id), l, m);
skymapuncal = acm2skyimage( Sigma , r_tile(:, 1), r_tile(:, 2), freq(freq_id), l, m);

save('sim_results\HBA.mat','skymapcal','skymapuncal','g_iter')





