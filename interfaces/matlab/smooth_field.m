function sv = smooth_field( g, v, lambda_stop, lambda_pass, order )

% Usage: sv = smooth_field( g, v, lambda_stop, lambda_pass, order )
%
% Applies Fourier smoothing to quantities on a periodic grid g. (g is an
% output of the load_brick or load_domain functions.) All wavelengths shorter
% than lambda_stop are removed from the data and all wavelengths longer
% than lambda_pass are preserved. The quantitiy v can lie on any node,
% cell, edge or face mesh associated with the grid g. Order specifies
% the indexing of v. If it is not given yxz ordering ([2 1 3], meshgrid
% compatible) is assumed.
%
% Written by:
%   Kevin J. Bowers, Ph.D.
%   Plasma Physics Group (X-1)
%   Los Alamos National Lab
%   March 2004 -  Adapted from V4PIC version 1 smooth3d.m

if nargin<5, order = [ 2 1 3 ]; end

% Extract grid parameters
nx = g(2); ny = g(3); nz = g(4);
dx = g(6); dy = g(7); dz = g(8);

% Put data into xyz order
sv = ipermute(v,order);

% Throw out periodic copies
sz = size(sv);
sv = sv(1:nx,1:ny,1:nz);

% Extract relevant parameters from the input data
k_pass = 2*pi/lambda_pass;
k_stop = 2*pi/lambda_stop;

% Make the spatial frequency grid
% Note: k > pi/L is equivalent to negative k
kx = 2*pi*(0:nx-1)/nx; kx = kx - 2*pi*(kx>pi); kx = kx/dx;
ky = 2*pi*(0:ny-1)/ny; ky = ky - 2*pi*(ky>pi); ky = ky/dy;
kz = 2*pi*(0:nz-1)/nz; kz = kz - 2*pi*(kz>pi); kz = kz/dz;
[ mkx, mky, mkz ] = ndgrid( kx, ky, kz );

% Setup the low pass filter response
mkr2 = mkx.^2 + mky.^2 + mkz.^2; % Mesh of physical k magnitude
kp2 = k_pass^2;
ks2 = k_stop^2;
hk = (mkr2<kp2)*1 + ...
     (mkr2>=kp2 & mkr2<=ks2).*(ks2-mkr2)/(ks2-kp2);

% Filter FFT the raw data
sv = real(ifftn(hk.*fftn(sv)));

% Put back periodic copies
if sz(1)==nx+1, sv(nx+1,:,:) = sv(1,:,:); end;
if sz(2)==ny+1, sv(:,ny+1,:) = sv(:,1,:); end;
if sz(3)==nz+1, sv(:,:,nz+1) = sv(:,:,1); end;

% Permute into order
sv = permute( sv, order );
