function [Sigma] = Correlator(sigma,SNR,Lsignal,l,EEP,D_elem,P,d,k)

%   [Sigma] = Correlator(Lsignal,l,EEP,D_tile,P,M,d)
%
%   This function calculates the array covariance matrix for a
%   an array pointed at broadside. This function is written specifically
%   for the "sim_NoiseOnly" simulation. The function assumes a tile size
%   larger than M=1

%   arguments
%   sigma           : Power of calibration source
%   SNR             : SNR of timesample measurement
%   Lsignal         : Integer indicating the number of timesamples
%   l               : Direction cosine of calibration source
%   EEP             : EEP towards calibration source
%   D_elem          : Ntiles x 2 containing x and y coordinates of antennas
%   P               : Number of receive paths (number of tiles)
%   d               : Inter-element spacing in m
%   k               : Wave number in 1/m
%
%   returns
%   Sigma           : Array covariance matrix


% Nelis Wilke 10 February 2020
% Copyright 2020. See the LICENSE file at the top-level directory of this distribution.

% calibration source signal
S = (sqrt(sigma).*(randn(1, Lsignal) + 1i * randn(1, Lsignal)))/sqrt(2);

% Generate noise signals at element level
n_t_elem = ((1 / sqrt(SNR)) .* (randn(P, Lsignal) + 1i * randn(P, Lsignal)))/sqrt(2);

% Array response matrix
A = exp((-1i*k)*((D_elem) * l.'));

% Array signal matrix
s_arr = ((A*S) + n_t_elem);


Sigma = (1/Lsignal) * (s_arr * s_arr');

end


