function [g_iter,r1] = Self_Holography(iter,g_prev,Sigma)
%
%   [g_iter,r1] = Self_Holography(iter,g_prev,Sigma,a_cal)  
%
%   This function estimates the receive path gains of an array using a
%   single self-holography iteration
%
%
%   arguments
%   iter        : Current iteration number
%   g_prev      : Estimated gains of previous iteration
%   Sigma       : Array covariance matrix for current iteration, beamformed
%                 calibration source
%
%   returns
%   g_iter      : Estimated gains for current iteration  
%   r1          : Crosscorrelation of reference signal and individual
%                 antenna signals


%% Start SH calibration

% Initialise calibration variables
P = length(Sigma);
sigma_cal = 1;              % Calibration source power
r2 = diag(Sigma);

% update reference beam weight vector
w_ref = 1./ conj(g_prev);  % Beamform to calibration source
    
% update crosscorrelation vector
r1 = Sigma'*w_ref;       % Tile-reference beam crosscorrelations (r_xy)

% Compile "y-vector"
R = [r1;r2];

% Construct submatrices of the A-matrix using the current gain estimate
a1 = ((g_prev)'*w_ref)*eye(P)*sigma_cal;
a2 = diag(w_ref);
a3 = diag(conj(g_prev)*sigma_cal);
a4 = eye(P);

% compile A-mat
A = [a1 a2; a3 a4];

% Solve unknown gain and system noise vector ("x-vector")
x_est = A\R;

% update gain estimates
g_est = x_est(1:P);
% ensure proper phase referencing
g_est = g_est / (g_est(1) / abs(g_est(1)));
% perform amplitude normalization to facilitate comparison
g_est = g_est / mean(abs(g_est));

% apply averaging every second iteration
if (mod(iter, 2) == 0)
    g_est = (g_est + g_prev) / 2;
    % ensure proper phase referencing
    g_est = g_est / (g_est(1) / abs(g_est(1)));
    % perform amplitude normalization to facilitate comparison
    g_est = g_est / mean(abs(g_est));
end
g_iter = g_est;
 

end

