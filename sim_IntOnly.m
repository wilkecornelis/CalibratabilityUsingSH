
% This script tests self-holography on a 24-by-24 URA with varying tile
% sizes for different separation between the calibration source and
% interfering sources. As the aim of this script is to assess the impact of
% interference on calibration using self-holography, the simulation is done
% for zero instrumental noise. To completely isolate the effects of
% interference, exepcted values are used for the signal measurements. 
%
% Nelis Wilke, 18 November 2019
% reviewed by Stefan J. Wijnholds, 18 November 2019
% Copyright 2020. See the LICENSE file at the top-level directory of this distribution.

close all;
clear;

%% setup variables
c = 2.99792e8;              % speed of light (m/s)
f = 1e9;                    % Observation frequency (Hz)
lambda = c/f;               % Wavelength (m)
k = 2 * pi / lambda;        % Wave number in 1/m
d = lambda/2;               % Inter-element spacing in m

Niter = 10;                 % Number of solving iterations
Nelem = 576;                % Total number of elements in station
int_angles = 10:1:85;       % interferer angles in deg
Nangles = length(int_angles);

%% Import EEP
[~, ~, phi0] = readColData('data/EEP/phi0.txt',2,1,1);     % Antenna E-field at 1GHz (phi = 0 deg)

%% Source signals
% List of tile sizes to investigate
% Note that M here denotes the number of elements on a side of a tile
% instead of the total number of elements in the tile
M_list = [1 2 4];

%% Loop over different tile sizes
for tile_id = 1:length(M_list)
    
    M = M_list(tile_id);    % Number of elements on a side of a tile
    P = Nelem / M^2;        % Number of tiles (or receive paths)
    
    % true (complex) gain values
    g_true = ones(P, 1);% + 0.1 * (randn(P, 1) + 1i * randn(P, 1));
    % ensure proper phase referencing
    g_true = g_true / (g_true(1) / abs(g_true(1)));
    % amplitude normalization to facilitate comparison
    g_true = g_true / mean(abs(g_true));
    G = diag(g_true);
    
    % Initialise
    rxx = zeros(P,Nangles);                 % Autocorrelations of receive paths (with calibration source included)
    rxy = zeros(P,Nangles);                 % Tile-reference beam crosscorrelations (with calibration source included)
    rxx_c = zeros(P,Nangles);               % Autocorrelations of receive paths (with calibration source excluded)
    rxy_c = zeros(P,Nangles);               % Tile-reference beam crosscorrelations (with calibration source excluded)
    
    rxx_i = zeros(P,Nangles);               % Autocorrelations of receive paths (with calibration source excluded)
    rxy_i = zeros(P,Nangles);               % Tile-reference beam crosscorrelations (with calibration source excluded)
    
    g_mag = zeros(P,Nangles);               % Gain magnitude error after last iteration
    g_phase = zeros(P,Nangles);             % Gain phase error after last iteration
    g_rmse = zeros(1,Nangles);              % RMS gain errors at each interferer angle
    SIR = zeros(1,Nangles);                 % SIR at each interferer angle
    
    
    %% Setup reference coordinates of each tile (position of first element of each tile with respect to station)
    antenna_idx = 1;
    D_tile = zeros(P, 2);
    for m = 0:sqrt(P)-1
        for n = 0:sqrt(P)-1
            D_tile(antenna_idx, :) = M * d .* [m n];
            antenna_idx = antenna_idx + 1;
        end
    end
    %% Loop over intereferer angles
    for int_angle_idx = 1:Nangles
        disp(['Interferer angle is ' num2str(int_angles(int_angle_idx)) ' degrees']);
        
        % Setup direction cosines of sources
        lm(1,:) = [0 0];
        lm(2,:) = [sind(int_angles(int_angle_idx)) .* cosd(0), sind(int_angles(int_angle_idx)).*sind(0)];
        
        % Value of EEP at source locations
        EEP(1) = 1;                                     % EEP in the direction of the calibration source
        EEP(2) = phi0(int_angles(int_angle_idx)+91);    % EEP in the direction of the interfering source
        
        % Set source powers
        sigmas = [1;2];
        
        % Initialise
        a_c = zeros(P,1);
        a_i = zeros(P,1);
      
        % Calculate station covariance matrix
        for tile_ind = 1:P          % Loop over tiles
            elem_idx =1;            % Element index with respect to station
            for m = 0:M-1
                for n = 0:M-1
                    D_elem(elem_idx,:) = d.*[m n] + D_tile(tile_ind,:);             % Position of element with reference to the station
                    elem_idx = elem_idx +1;
                end
            end
            a_c(tile_ind,1) = sum(exp((-1i*k)*((D_elem) * lm(1,:).')));    % Gemotric delay of calibration source at current tile
            a_i(tile_ind,1) = sum(exp((-1i*k)*((D_elem) * lm(2,:).')));    % Geometric delay of interferer at current tile 
        end
        % Calculate station covariance matrix
        Sigma = sigmas(1)*(a_c*a_c')+(sigmas(2)*EEP(2))*(a_i*a_i');                % Array convariance matric of calibration source
        Sigma_i = (sigmas(2)*EEP(2))*(a_i*a_i');                         % Array convariance matric of interfereing source
              
        %% Start SH calibration
        
        % Initial gain estimates
        g_est_prev = ones(P,1);         % Estimate of previous iteration
        g_est = ones(P,1);              % Estimate of current iteration
        rxy_iter = 0;                   % Crosscorrelations after last iteration
               
        for iter = 1:Niter
            % Estimate gains using SH
            [g_est,rxy_iter] = Self_Holography_iter(iter,g_est_prev,Sigma);
            g_est_prev = g_est;
        end
        
        % Apply true gains to intereferer covariance matrix
        Sigma_i = (G*Sigma_i*G');            
        
        rxx(:,int_angle_idx) = abs(diag(Sigma));
        rxy(:,int_angle_idx) = rxy_iter;
        
        rxy_c(:,int_angle_idx) = (Sigma - Sigma_i)*ones(P,1);
        rxx_c(:,int_angle_idx) = abs(diag(Sigma-Sigma_i));
        
        rxy_i(:,int_angle_idx) = Sigma_i*ones(P,1);
        rxx_i(:,int_angle_idx) = abs(diag(Sigma_i));
        
        g_mag(:,int_angle_idx) = (abs(abs(g_true)-abs(g_est)));
        g_phase(:,int_angle_idx) = (abs(angle(g_true)-angle(g_est)));
        g_rmse(1,int_angle_idx) = sqrt(mean(abs(g_true - g_est).^2));
        
        % Calculate SIR
        SIR(1,int_angle_idx) = (mean(abs(rxy_c(:,int_angle_idx)-rxx_c(:,int_angle_idx)))) ./ (mean(abs(rxy_i(:,int_angle_idx)-rxx_i(:,int_angle_idx))));       
    end

    % Sort SIR in ascending order
    [SIR, sortorder] = sort(SIR);
    
    % Save results
    save(['sim_results/results_IntOnly_M' num2str(M_list(tile_id)) '.mat'], ...
        'rxy_i','rxx_i','g_mag', 'g_phase','g_rmse','SIR','sortorder');
end






