function [Sigma] = Correlator_tiles(sigmas,SNR,Lsignal,lm,EEP,D_elem,P,d,k)

%   [Sigma] = Correlator(Lsignal,l,EEP,D_tile,P,M,d)
%
%   This function calculates the array covariance matrix for a given source
%   distribution. The array can be subdivided into tiles by setting M > 1
%   Note that this function can only be applied to a regular dense array with variable tile size.

%   arguments
%   sigmas          : Vector containing source powers
%   SNR             : SNR of timesample measurement
%   Lsignal         : Integer indicating the number of timesamples
%   lm               : Nsrc x 2 direction cosine matrix of sources
%   EEP             : Nsrc x 1 vector containg gains towards sources
%   D_elem          : Nant x 2 containing x and y coordinates of antennas
%   P               : Number of receive paths (number of tiles)
%   M               : Number of elements on a side of a tile
%   d               : Inter-element spacing in m
%   k               : Wave number in 1/m
%
%   returns
%   Sigma           : Array covariance matrix

% CRW 10 February 2020

% Calibration source signal
s_c = (sqrt(sigmas(1)).*(randn(1, Lsignal) + 1i * randn(1, Lsignal)))/sqrt(2);
% Interferer source signal
s_i = (sqrt(sigmas(2)).*sqrt(EEP(2)).*(randn(1, Lsignal) + 1i * randn(1, Lsignal)))/(sqrt(2));

% Stack signals in single matrix
S(1,:) = s_c;
S(2,:) = s_i;
  
% Generate noise signals at tile level
n_t_elem = (sqrt(1 /(SNR)) .* (randn(P, Lsignal) + 1i * randn(P, Lsignal)))/sqrt(2);

% Array response matrix
A = exp((-1i*k)*((D_elem) * lm.'));

% Output signal for this tile
s_elem = (A*S) + n_t_elem;

% Calculate station covariance matrix
Sigma = (1/Lsignal) * (s_elem * s_elem');

end


