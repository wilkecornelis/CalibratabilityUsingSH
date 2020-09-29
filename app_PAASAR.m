% This code simulates SH calibration on the PAASAR concept using GPS satellites from the European
% and American constellations as calibration sources.

% Nelis Wilke 18 May 2020
% Copyright 2020. See the LICENSE file at the top-level directory of this distribution.

clear
close all

%% Setup constants
Ndates = 72;                % Number of date recordings
Nsats = 56;                 % Number of satellites

f = 1544.1e6;                % Observation frequency (Hz)
c = 2.99792e8;               % speed of light (m/s)
lambda = c/f;                % Wavelength (m)
d = lambda/2;                % Inter-element spacing in m
k = 2 * pi / lambda;         % Wave number in 1/m
P = 50;                      % Total number of elements in array

B = 10e6;                      % Integration bandwidth in Hz
tau = 0.001;                   % Integration time in s
Lsignal = B * tau;             % Number of time samples

Nruns = 1;                     % Number of calibration runs
Niter = 10;                    % Number of SH solving iterations

% true (complex) gain values
g_true = ones(P, 1)+ 0.1 * (randn(P, 1) + 1i * randn(P, 1));
% ensure proper phase referencing
g_true = g_true / (g_true(1) / abs(g_true(1)));
% amplitude normalization to facilitate comparison
g_true = g_true / mean(abs(g_true));
G = diag(g_true);

% Import satellite positions
load('data/PAASAR/GPSsat_altaz.mat');
load('data/PAASAR/GPSsat_lm.mat');

% Import coordinates of PAASAR antenna elements
antpos = load('data/PAASAR/PAASARantpos.mat');
antpos = antpos.pos;

% Import normal vectors of PAASAR triangles. A triangle corresponds to a
% segment of the icosahedron sphere.
tri_n = load('data/PAASAR/PAASARnormals.mat');
tri_n = tri_n.tri_n;

%% Initialise

SIR = NaN(Nsats,Ndates);         % SIRs for all sources on each date
SIR_max =  zeros(Ndates,1);        % Maximum SIR on each date

g_rmse = NaN(Nsats,Ndates);      % RMS gain estimation error
g_phase = NaN(Nsats,Ndates);     % Gain phase estimation error
g_mag = NaN(Nsats,Ndates);       % Gain magnitude estimation error

%% Start simulation

% Loop over all dates
for date_id = 1:Ndates
    disp(' ')
    disp(['Date ID: ' num2str(date_id) ' of ' num2str(Ndates)]);
    disp(' ')
    
    % clear variables 
    clear lmn_sel
    clear a_sat
    clear EEP
    clear EEP_tri
    clear alt_sel
    clear az_sel
    clear sat_cart
    clear calsat_ids
    
    % Identify indices of satellites that are above the horizon. These sources
    % are further referred to as "selected sources". "Satellites" and
    % "sources" are used interchangeably throughout the code.
    sat_above_ids = find(isnan(sat_altaz(:,date_id)) == 0);
    
    % store altaz angles of selected source (_sel refers to selected)
    altaz_sel = [rad2deg(sat_altaz(sat_above_ids,date_id,1)) rad2deg(sat_altaz(sat_above_ids,date_id,2))];
    % store direction cosine lmn coordinates of selected sources in a single matrix
    lmn_sel = [sat_lm(sat_above_ids,date_id,1) sat_lm(sat_above_ids,date_id,2) sqrt(1 - sat_lm(sat_above_ids,date_id,1).^2 - sat_lm(sat_above_ids,date_id,2).^2)];
    % calculate geometric delay vectors of selected sources
    a_sat = exp(2*pi*-1j*f*(antpos*lmn_sel.')/c);
    
    % Calculate EEP attenuation for each triangle.
    % EEP = cos of the angle between normal vector of each triangle and unit
    % vector pointing to each source.
    EEP_tri = cosd(acosd(tri_n.'*lmn_sel.'));
    % Set EEP to zero for sources that are below horizons of antennas
    EEP_tri(EEP_tri<0)=0;
    % Set antenna EEPs
    EEP(1:10,:) =  repmat(EEP_tri(1,:),10,1);
    EEP(11:20,:) = repmat(EEP_tri(2,:),10,1);
    EEP(21:30,:) = repmat(EEP_tri(3,:),10,1);
    EEP(31:40,:) = repmat(EEP_tri(4,:),10,1);
    EEP(41:50,:) = repmat(EEP_tri(5,:),10,1);
    
    % Apply EEP attenuation to magnitudes of satellite delay vectors
    a_sat = a_sat.*sqrt(EEP);
    
    % Get indices of sources that are visible to all antennas simultaneously
    calsat_ids = find(sum(EEP_tri>0)==5);
    
    % Loop over sources that are simultaneously visible to all antennas.
    % These sources will further be referred to as "visible sources".    
    for idx = 1:length(calsat_ids)
        disp(['working on source: ' num2str(calsat_ids(idx)) ' of ' num2str(length(sat_above_ids))]);

        % Define beamforming weights for current source
        a_bf = a_sat(:,calsat_ids(idx));
        
        % Covariance matrix of full sky
        Sigma_s = a_sat(:,:)*a_sat(:,:)';
        
        % Covariance matrix of calibration source only
        Sigma_c = (a_sat(:,calsat_ids(idx))*a_sat(:,calsat_ids(idx))');
        % Beamform to calibration source
        Sigma_c = diag(a_bf')*Sigma_c*diag(a_bf);
        % Calculate numerator of SIR
        num = mean(abs((ones(P,1)'*(Sigma_c)).'-diag(Sigma_c)));
        
        % Covariance matrix of interferers
        Sigma_i = diag(a_bf')*Sigma_s*diag(a_bf')' - Sigma_c;
        % Calculate denominator of SIR
        den = mean(abs((ones(P,1).'*(Sigma_i)).'-diag(Sigma_i)));
        
        % Calculate SIR (dB) for this source on this date
        SIR(calsat_ids(idx),date_id) = 10*log10(num/den);
%         disp(['SIR: ' num2str(SIR)])
        
        %% SH calibration
        % Calibrate only when SIR > 6 dB
        if SIR(calsat_ids(idx),date_id) > 6  
            % Initialise calibration variables
            g_rmse_run = zeros(P,Nruns);
            g_phase_run = zeros(P,Nruns);
            g_mag_run = zeros(P,Nruns);

            % Set SNR in dB for signal model
            SNR = 2;            
            % Set source powers for signal model
            sigmas = ones(length(sat_above_ids),1);   

            for run = 1:Nruns
                % Initialise gain estimates
                g_prev = ones(P,1);     
                for iter = 1:Niter
                    % Generate a signal measurement for this iteration
                    [Sigma] = Correlator_PAASAR(sigmas,10^(SNR/10),Lsignal,lmn_sel,EEP,P,antpos,k);
                    % Apply true gains and beamforming weights
                    Sigma = (G*Sigma*G');
                    % Beamform to calibration source
                    Sigma = diag(a_bf')*Sigma*diag(a_bf);

                    % Estimate gains using SH
                    [g_est,r1] = Self_Holography_iter(iter,g_prev,Sigma);
                    g_prev = g_est;
                end
                % Gain errors for this run
                g_rmse_run(:,run) = sqrt(mse(abs(g_true-g_est)));
                g_phase_run(:,run) = abs(angle(g_true)-angle(g_est));
                g_mag_run(:,run) = abs(abs(g_true)-abs(g_est));
            end
            % Gain errors averaged over all runs
            g_rmse(calsat_ids(idx),date_id) = mean(mean(g_rmse_run,2));    % RMS gain error
            g_phase(calsat_ids(idx),date_id) = mean(mean(g_phase_run,2));  % Phase error
            g_mag(calsat_ids(idx),date_id) = mean(mean(g_mag_run,2));      % Magnitude error
        end
    end
    SIR_max(date_id) = max(SIR(:,date_id));
end


% save(['sim_results/results_PAASAR.mat'], ...
%     'g_rmse','g_mag', 'g_phase','SIR','SIR_max');










