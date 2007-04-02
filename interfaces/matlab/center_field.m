function cv = center_field( g, v, method, order )

% Usage: cv = center_field( g, v, method, order )
%
% Node centers v on a periodic grid g. (g is an output of the load_brick or
% load_domain functions.) The quantitiy v can lie on any node, cell, edge or
% face mesh associated with the grid g. If method not 1, the quantity is node
% centered by averaging. If method is 1, the quantity is node centered by
% spectral methods. Order specifies the indexing of v. If it is not given yxz
% ordering ([2 1 3], meshgrid compatible) is assumed.
%
% Written by:
%   Kevin J. Bowers, Ph.D.
%   Plasma Physics Group (X-1)
%   Los Alamos National Lab
%   April 2004 -  Original version written

if nargin<3, method = 0; end
if nargin<4, order = [ 2 1 3 ]; end

% Extract grid parameters
nx = g(2); ny = g(3); nz = g(4);

% Put data into xyz order
cv = ipermute(v,order);

% Throw out periodic copies
sz = size(cv);
cv = cv(1:nx,1:ny,1:nz);

% Compute some indexing arrays to make life easier
ii = 1:nx; im1 = [ nx 1:nx-1 ];
jj = 1:ny; jm1 = [ ny 1:ny-1 ];
kk = 1:nz; km1 = [ nz 1:nz-1 ];

% Center the fields by averaging
if sz(1)==nx, cv = 0.5*( cv(im1,:,:) + cv ); end
if sz(2)==ny, cv = 0.5*( cv(:,jm1,:) + cv ); end
if sz(3)==nz, cv = 0.5*( cv(:,:,km1) + cv ); end

if method==1,
  % Compute averaging filter along each direction
  gx = ones(1,nx);
  gy = ones(1,ny);
  gz = ones(1,nz);
  if sz(1)==nx, gx = abs(cos(pi*(ii-1)/nx)); end
  if sz(2)==ny, gy = abs(cos(pi*(jj-1)/ny)); end
  if sz(3)==nz, gz = abs(cos(pi*(kk-1)/nz)); end

  % Compute the anti-averaging filter along each direction
  if sz(1)==nx && mod(nx,2)==0, gx(1+nx/2) = 1; end; % Avoid divide by zero
  if sz(2)==ny && mod(ny,2)==0, gy(1+ny/2) = 1; end;
  if sz(3)==nz && mod(nz,2)==0, gz(1+nz/2) = 1; end;
  gx = 1./gx;
  gy = 1./gy;
  gz = 1./gz;
  if sz(1)==nx && mod(nx,2)==0, gx(1+nx/2) = 0; end; % Handle lost info
  if sz(2)==ny && mod(ny,2)==0, gy(1+ny/2) = 0; end;
  if sz(3)==nz && mod(nz,2)==0, gz(1+nz/2) = 0; end;

  % Assemble the net filter response
  [ gx, gy, gz ] = ndgrid( gx, gy, gz );

  % Compute the spectrally correct centered value
  cv = real(ifftn(gx.*gy.*gz.*fftn(cv)));
end

% Put back periodic copies
cv(nx+1,:,:) = cv(1,:,:);
cv(:,ny+1,:) = cv(:,1,:);
cv(:,:,nz+1) = cv(:,:,1);

% Permute into output order
cv = permute( cv, order );
