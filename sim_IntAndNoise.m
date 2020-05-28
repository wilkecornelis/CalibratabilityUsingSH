%% Test impact of SIR on self-holography with varying tile sizes
%
% This script tests self-holography on a 24-by-24 URA for 
% different separation between the calibration source and an
% interfering source. The SNR of the signal measurements are also varied to
% determine the performance as a function of interference and noise level,
%
% Nelis Wilke, 28 March 2020


close all
clear 

%% setup variables
c = 2.99792e8;              % speed of light (m/s)
f = 1e9;                    % Observation frequency (Hz)
lambda = c/f;               % Wavelength (m)
k = 2 * pi / lambda;        % Wave number in 1/m
d = lambda/2;               % Inter-element spacing in m

B = 10e6;                              % Integration bandwidth in Hz
tau = 0.0001;                          % Integration time in s
Lsignal = B * tau;                     % Number of time samples

Niter = 2;                              % Number of solving iterations
Nelem = 576;                            % Total number of elements in station
int_angles = 10:1:85;                   % interferer angles in deg
snr_vals = (1:1:150)/10000;             % SNR values to investigate
Nangles = length(int_angles);           % Number of interferer angles to investigate
Nsnr = length(snr_vals);                % Number of SNR values to investigate
P = 576;                                % Number of receive paths

% Import EEP
[~, ~, phi0] = readColData('data/EEP/phi0.txt',2,1,1);     % Dipole E-field at 1GHz (phi = 0 deg)

% True (complex) gain values
g_true = ones(P, 1) + 0.1 * (randn(P, 1) + 1i * randn(P, 1));
% Ensure proper phase referencing
g_true = g_true / (g_true(1) / abs(g_true(1)));
% Amplitude normalization to facilitate comparison
g_true = g_true / mean(abs(g_true));
G = diag(g_true);

% Initialise
g_mag = zeros(P,Nangles,Nsnr);               % Gain magnitude error after last iteration
g_phase = zeros(P,Nangles,Nsnr);             % Gain phase error after last iteration

%% Setup reference coordinates of each tile (position of first element of each tile with respect to station)
ind = 1;
D_elem = zeros(P, 2);
for m = 0:sqrt(P)-1
    for n = 0:sqrt(P)-1
        D_elem(ind, :) = d .* [m n];
        ind = ind + 1;
    end
end

%% Loop over interferer angles
for int_angle_idx = 1:Nangles
    tic
    disp(['Interferer angle is ' num2str(int_angles(int_angle_idx)) ' degrees']);
    
    % Setup source distribution
    l_s = zeros(2,2);
    l_s(1,:) = [0 0];
    l_s(2,:) = [sind(int_angles(int_angle_idx)) .* cosd(0), sind(int_angles(int_angle_idx)).*sind(0)];
    
    % Value of EEP at source locations
    EEP(1) = 1;
    EEP(2) = phi0(int_angles(int_angle_idx)+91);
    
    for snr_idx =1:Nsnr
        disp(['SNR index is ' num2str(snr_idx) ' of ' num2str(length(snr_vals))]);
        
        % Set SNR
        SNR = snr_vals(snr_idx);
        
        %% Start SH calibration
        % Initialise calibration variables
        g_est_prev = ones(P,1);          % Initial gain estimates
        sigmas = [1;4];                  % source powers
        
        for iter = 1:Niter
            % Generate a signal measurement for this iteration
            [Sigma] = Correlator_IntAndNoise(sigmas,SNR,Lsignal,l_s,EEP,D_elem,P,d,k);
            % Apply true gains
            Sigma = (G*Sigma*G');
            
            % Estimate gains using SH
            [g_est,r1] = Self_Holography_iter(iter,g_est_prev,Sigma);
            g_est_prev = g_est;
        end
        g_mag(:,int_angle_idx,snr_idx) = abs(abs(g_est) - abs(g_true));
        g_phase(:,int_angle_idx,snr_idx) = abs(angle(g_est) - angle(g_true));
        
    end
    toc
end
save(['sim_results/results_IntAndNoise.mat'], ...
    'g_mag','g_phase','snr_vals');


