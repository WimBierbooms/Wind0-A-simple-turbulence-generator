function [t,UC]=wind0(yr,zr,U,sigma,N,deltat,fmax)
% syntax: function [t,UC]=wind0(yr,zr,U,sigma,N,deltat,fmax)
% simulation of a turbulent wind field
%
% INPUT:
%    yr, zr: specification of RECTANGULAR GRID of rotorpoints; row vectors
%    U: mean wind velocity (m/s)
%    sigma: standard deviation (m/s)
%    N: number of time points (including zero); N must be a power of 2 
%    deltat: time step (s)
%    fmax: maximum frequency spectrum (Hz); spectrum is cut-off above fmax
% OUTPUT:
%    t: time (s)
%    UC: turbulent wind velocities (m/s)

% number of points in rotor plane
Ny=length(yr);
Nz=length(zr);
Np=Ny*Nz;

% y and z coordinates of all rotor points in one column vector
Yr=reshape(yr'*ones(1,Nz),Np,1);
Zr=reshape(ones(Ny,1)*zr,Np,1);

r=zeros(Np,Np);
for i=1:Np
   for j=i+1:Np
      % distances between points 
      r(i,j)=sqrt((Yr(i)-Yr(j))^2+(Zr(i)-Zr(j))^2);
      r(j,i)=r(i,j);
   end
end

% time vector
t=[0:N-1]'*deltat;
% period
T=N*deltat;
% frequency step
deltaf=1/T;
% discretized frequencies
k=[1:N/2-1]';
f=k.*deltaf;
% autopower spectral density (one-sided)
Sa=Autopow(f,U,sigma);
% spectrum is cut-off above fmax by application of window
Index=find(f>fmax);
if ~isempty(Index)
  Nw=Index(1);
  w=zeros(N/2-1,1);
  W=window('hann',2*Nw+1);w(1:Nw+1)=W(Nw+1:2*Nw+1);
  Sa=w.*Sa;
end
% renormalize Sa to variance
Sa=sigma^2/(sum(Sa)/T)*Sa;

% Fouriercoefficients points in rotor plane
ak=zeros(Np,N/2-1);
bk=zeros(Np,N/2-1);
for k=1:N/2-1
   Coh=Coher(f(k),r,U,50,2);
   % Choleski decomposition
   L=sqrt(Sa(k)/T)*chol(Coh)';
   % vector of unit variance normal random numbers
   ran=randn(Np,1);
   ak(:,k)=L*ran;
   ran=randn(Np,1);
   bk(:,k)=L*ran;
end

% complex notation
i=sqrt(-1);
UC=zeros(N,Np);
for j=1:Np
   C=ak(j,:)'-i*bk(j,:)';
   C=1/2*[0;C;0;rot90(C')];
   % inverse FFT
   uc=N*ifft(C);
   if any(abs(imag(uc)) >= 1e-7*abs(uc) & abs(imag(uc)) >= 1e-12)
     max(abs(uc))
     max(imag(uc))
     error('imag too large uc')
   end
   UC(:,j)=real(uc);
end
% reshape UC: separate indices for y and z
UC=reshape(UC,N,Ny,Nz);
