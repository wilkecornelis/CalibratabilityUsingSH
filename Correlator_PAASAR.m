function [Sigma] = Correlator(sigmas,SNR,Lsignal,lmn,EEP,P,D,k)

%   [Sigma] = Correlator(Lsignal,l,EEP,D_tile,P,M,d)
%
%   This function calculates the array covariance matrix for PAASAR..

%   arguments
%   sigmas          : Vector containing source powers
%   SNR             : SNR of timesample measurement
%   Lsignal         : Integer indicating the number of timesamples
%   lmn             : Nsrc x 3 direction cosine matrix of sources
%   EEP             : Nsrc x 1 vector containg gains towards sources
%   D               : Nant x 2 containing x and y coordinates of antennas
%   P               : Number of receive paths (number of tiles)
%   k               : Wavenumber in 1/m

%   returns
%   Sigma           : Array covariance matrix

% CRW 10 February 2020     

% Matrix of source signals
s_t = (sqrt(sigmas).*((randn(length(lmn), Lsignal) + 1i * randn(length(lmn), Lsignal))))/sqrt(2);  

% Matrix of noise signals
n_t = (1 / sqrt(SNR)) .* (randn(P, Lsignal) + 1i * randn(P, Lsignal));     

% Array response matrix
A = exp((-1i*k)*(D * lmn.')).*sqrt(EEP);

% matrix of element output signals
s_elem = A * s_t + n_t;

% Calculate station covariance matrix
Sigma = (1/Lsignal) * (s_elem * s_elem');
end


