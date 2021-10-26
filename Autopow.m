function S=autopow(f,U,sigma)
% syntax: function S=autopow(f,U,sigma)
% Autopower spectral density function of turbulence
% Input:
%   f: frequency (Hz)
%   U: the (10 minute) mean wind speed (m/s)
%   sigma: standard deviation (m/s)
% Output:
%   S: autopower spectral density (m^2/s)


% Kaimal spectrum according to the IEC standard       
% turbulence scale parameter (H >= 60m)
Lambda1=42;
% isotropic integral scale parameter
L1=8.1*Lambda1;
S=sigma.^2.*4.*L1./U./(1+6.*f.*L1./U).^(5/3);
