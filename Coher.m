function Coh=coher(f,r,Vhub,zhub,option)
% syntax: function Coh=coher(f,r,Vhub,zhub,option)
% Coherency function
% of longitudinal wind velocity fluctuations
% Implemented options:
%   1. von Karman coherency according to the IEC standard
%   2. Kaimal coherency according to the IEC standard
% Input:
%   f: frequency (Hz)
%   r: distance (in projection of rotor plane)
%   Vhub: the 10 minute average wind speed at hub height (m/s)
%   zhub: the hub height of the wind turbine (m)
% Output:
%   Coh: coherency (-)
% Author: Wim Bierbooms
% Last update: 17-4-'96
%

if (option==1)
% von Karman coherency according to the IEC standard       
% turbulence scale parameter
  Lambda1=min(0.678*zhub,20.3);
% isotropic integral scale parameter
  L1=3.484*Lambda1;
  x=2.*pi.*((f.*r./Vhub).^(2)+(0.12.*r./L1).^(2)).^(0.5);
  % exception for x=0
  Coh=ones(size(x));  
  Ind=find(x~=0);
  Coh(Ind)=2^(1/6)./gamma(5/6).*(x(Ind).^(5/6).*besselk(5/6,x(Ind))-...
                           0.5.*x(Ind).^(11/6).*besselk(1/6,x(Ind)));
elseif (option==2)
% Kaimal coherency according to the IEC standard       
% turbulence scale parameter
  Lambda1=min(0.7*zhub,42);
% isotropic integral scale parameter
  Lc=8.1*Lambda1;
  x=((f.*r./Vhub).^(2)+(0.12.*r./Lc).^(2)).^(0.5);
  Coh=exp(-12.*x);
else
  error 'option not implemented in COHER'
end  
