% This script tests self-holography on a 24-by-24 URA
% for varying integration times and instantaneous SNR at element level.

% Nelis Wilke, 28 March 2020
% Copyright 2020. See the LICENSE file at the top-level directory of this distribution.

close all
clear

%% setup variables
c = 2.99792e8;              % speed of light (m/s)
f = 1e9;                    % Observation frequency (Hz)
lambda = c/f;               % Wavelength (m)
k = 2 * pi / lambda;        % Wave number in 1/m
d = lambda/2;               % Inter-element spacing in m
B = 10e6;                   % Integration bandwidth in Hz

Niter = 3;                    % Number of SH solving iterations
Nruns = 100;                  % Number of calibration runs

snr_vals = (1:1:150)/10000;                 % SNR values to investigate
intTime_vals = [0.0001,0.0002,0.0004];      % Integration times to investigate

Nsnr = length(snr_vals);                    % Number of SNR values to investigate
Ntimes = length(intTime_vals);              % Number of integration times to investigate

P = 576;                                    % Total number of elements in the array

%% Setup antenna element coordinates
ind = 1;
D_elem = zeros(P, 2);
for m = 0:sqrt(P)-1
    for n = 0:sqrt(P)-1
        D_elem(ind, :) = d .* [m n];
        ind = ind + 1;
    end
end

% Direction cosine of calibration source
lm_s(1,:) = [0 0];
% Value of EEP at source locations
EEP = 1;
% true (complex) gain values
g_true = ones(P, 1) + 0.1 * (randn(P, 1) + 1i * randn(P, 1));
% ensure proper phase referencing
g_true = g_true / (g_true(1) / abs(g_true(1)));
% amplitude normalization to facilitate comparison
g_true = g_true / mean(abs(g_true));
G = diag(g_true);

%% Loop over different integration durations
g_est_run = zeros(P,Nruns);
g_rmse = zeros(Nsnr,Ntimes);
for intTime_idx = 1:Ntimes
    tau = intTime_vals(intTime_idx);                 % Integration time in s
    Lsignal = B * tau;                               % Number of time samples
    
    %% Loop over different SNR values
    for snr_idx =1:Nsnr
        SNR = snr_vals(snr_idx);
        disp(['SNR index is ' num2str(snr_idx) ' of ' num2str(Nsnr)]);
        
        for run = 1:Nruns      
            
            %% Start SH calibration           
            % Initialise calibration variables
            g_est_prev = ones(P,1);               % Initial gain estmates
            sigma = 1;                           % Power of calibration source           
            for iter = 1:Niter
                % Generate a signal measurement for this iteration
                [Sigma] = Correlator_NoiseOnly(sigma,SNR,Lsignal,lm_s,EEP,D_elem,P,d,k);   
                % Apply true gains
                Sigma = (G*Sigma*G');
                % Estimate gains using SH
                [g_est,r1] = Self_Holography_iter(iter,g_est_prev,Sigma);
                g_est_prev = g_est;
            end            
            % Store gain estimates for this run.
            g_est_run(:,run) = g_est;
        end
        g_est = mean(g_est_run,2);
        g_rmse(snr_idx,1) = sqrt(mean(abs(g_true - g_est).^2));
    end
    save(['sim_results/results_NoiseOnlyT',num2str(intTime_idx),'.mat'], ...
    'g_rmse','Lsignal','snr_vals');
end






