function [ phi, ax, ay, az, rho ] = ...
  gauge_fields( g, ex, ey, ez, bx, by, bz, order );

% Usage: [ phi, ax, ay, az, rho ] = ...
%   gauge_fields( g, ex, ey, ez, bx, by, bz, order );
%
% Computes the scalar potential, the vector potential and the microscopic
% charge density associated with the fields ex, ey, ez, bx, by, bz
% on grid g (g is an output of the load_domain or load_brick functions).
% The potentials are in the Coulomb gauge and the integral of each potential
% over the region is zero. The input and output fields are on a Yee-mesh and
% in the order given. If order is not specified, yxz ordering ([2 1 3]) is
% assumed (load_domain, load_brick defaults and meshgrid compatible).
%
% The scalar potential calculation is based on the following:
%   grad PHI = -E
%   => div grad PHI = -div E
%   => laplacian PHI = -div E
%
% The vector potential calculation is based on the following:
%   div A = 0     (Coulomb gauge condition)
%   => A = curl G (The divergence of a curl is automatically zero)
%   => B = curl A = curl curl G = grad div G - laplacian G
% We have the freedom to pick div G = 0. Then:
%   laplacian G = -B
%   A = curl G
%
% The microscopic rho calculation is based on:
%   rho = eps0 div E
%   
% The Laplacian equation discretized with the usual 7-point stencil. This
% stencil is consistent with usual Yee mesh curl and divergence relations.
% The Laplacian is inverted via Fourier methods.
%
% Written by:
%   Kevin J. Bowers, Ph.D.
%   Plasma Physics Group (X-1)
%   Los Alamos National Lab
%   April 2004 - Adapted from V4PIC version 1 load_fields.m gauge calculation

if nargin<8, order = [ 2 1 3]; end

% Extract the relavant parameters from the grid
nx = g(2);
ny = g(3);
nz = g(4);

dx = g(6);
dy = g(7);
dz = g(8);

eps0 = g(10);

% Put the input fields into xyz order
if any(order~=[1 2 3]),
  ex = ipermute( ex, order );
  ey = ipermute( ey, order );
  ez = ipermute( ez, order );

  bx = ipermute( bx, order );
  by = ipermute( by, order );
  bz = ipermute( bz, order );
end

% Truncate periodic values of input fields
ex = ex(1:nx,1:ny,1:nz);
ey = ey(1:nx,1:ny,1:nz);
ez = ez(1:nx,1:ny,1:nz);

bx = bx(1:nx,1:ny,1:nz);
by = by(1:nx,1:ny,1:nz);
bz = bz(1:nx,1:ny,1:nz);

% Compute some indexing arrays to make life easier
ii = 1:nx; ip1 = [ 2:nx 1 ]; im1 = [ nx 1:nx-1 ];
jj = 1:ny; jp1 = [ 2:ny 1 ]; jm1 = [ ny 1:ny-1 ];
kk = 1:nz; kp1 = [ 2:nz 1 ]; km1 = [ nz 1:nz-1 ];

% Compute Laplacian kernel
gx=((2/dx)*sin((pi/nx)*(ii-1))).^2;
gy=((2/dy)*sin((pi/ny)*(jj-1))).^2;
gz=((2/dz)*sin((pi/nz)*(kk-1))).^2;
[ gx, gy, gz ] = ndgrid( gx, gy, gz );
kernel = gx + gy + gz; % Compute discretized equivalent of k^2
kernel(1,1,1) = 1;     % Avoid a divide by zero error at k=0
kernel = 1./kernel;    % Discretized equivalent of 1/k^2
kernel(1,1,1) = 0;     % Enforce integral phi d^3 x = 0

% Compute scalar potential
rho = ( ex - ex(im1,:,:) )*(1/dx) + ...
      ( ey - ey(:,jm1,:) )*(1/dy) + ...
      ( ez - ez(:,:,km1) )*(1/dz);
phi = real(ifftn(kernel.*fftn(rho)));

% Compute vector potential
gx = real(ifftn(kernel.*fftn(bx)));
gy = real(ifftn(kernel.*fftn(by)));
gz = real(ifftn(kernel.*fftn(bz)));
ax = ( gz-gz(:,jm1,:) )*(1/dy) - ( gy-gy(:,:,km1) )*(1/dz);
ay = ( gx-gx(:,:,km1) )*(1/dz) - ( gz-gz(im1,:,:) )*(1/dx);
az = ( gy-gy(im1,:,:) )*(1/dx) - ( gx-gx(:,jm1,:) )*(1/dy);

% Compute rho
rho = eps0*rho;

% Put periodic copies back into output
phi(nx+1,:,:) = phi(1,:,:);
phi(:,ny+1,:) = phi(:,1,:);
phi(:,:,nz+1) = phi(:,:,1);

ax(:,ny+1,:) = ax(:,1,:);
ay(:,:,nz+1) = ay(:,:,1);
az(nx+1,:,:) = az(1,:,:);

ax(:,:,nz+1) = ax(:,:,1);
ay(nx+1,:,:) = ay(1,:,:);
az(:,ny+1,:) = az(:,1,:);

rho(nx+1,:,:) = rho(1,:,:);
rho(:,ny+1,:) = rho(:,1,:);
rho(:,:,nz+1) = rho(:,:,1);

% Permute into order
if any(order~=[1 2 3]),
  phi = permute( phi, order );
  ax  = permute( ax,  order );
  ay  = permute( ay,  order );
  az  = permute( az,  order );
  rho = permute( rho, order );
end

