function [ g, ...
           rho,jx,jy,jz, ...
           ke,px,py,pz, ...
           txx,tyy,tzz, ...
           tyz,tzx,txy ...
         ] = load_domain_hydro(filename,order);

% Usage: [ g, ...             % Grid and species parameters
%          RHO, ...           % Charge density
%          Jx, Jy, Jz, ...    % Current density (charge flux)
%          K, ...             % Kinetic energy density
%          Px, Py, Pz, ...    % Momentum density (energy flux)
%          Txx, Tyy, Tzz, ... % Stress tensor diagonal
%          Tyz, Tzx, Txy ...  % Stress tensor off-diagonal
%        ] = load_domain_hydro(filename,order);
%
% g - Grid parameters:
%     g = [ nt nx ny nz ...
%           dt dx dy dz ...
%           cvac eps0 damp ...
%           x0 y0 z0 ...
%           spid spqm ...
%           rank ndom ]
%   where:
%     nt  nx, ny, nz - Time level and grid resolution
%     dt, dx, dy, dz - Time step and grid spacing
%     cvac           - Speed of light
%     eps0           - Permittivity of free space
%     damp           - Radiation damping parameter
%     x0, y0, z0     - Offset of this domain
%     spid spqm      - Species ID and the charge to mass ratio
%     rank           - ID of this domain
%     ndom           - Total number of domains (i.e. nproc)
%
% filename - Name of the hydro dump file to load.
%
% order - (optional) Indicates the desired field indexing.
%   [2 1 3] (default) - YXZ indexing (compatible with "meshgrid"
%                       and MATLAB's 3d plotting routines)
%   [1 2 3] - XYZ indexing (compatible with "ndgrid")
%   Other orderings are easy to figure out.
%
% Notes:
% - All values are time centered appropriately. However, as a
%   result, the Jx, Jy and Jz returned are not the actual
%   fields used to advance the simulation in time. Further, they
%   will not exactly satisfy local discretized continuity of
%   charge (global continuity of charge is satified). The Jx, Jy
%   and Jz used internally by the simulation do satisfy local 
%   continuity of charge exactly but are staggered in time.
% - All fields are in the appropriate units for the problem.
%   To understand the problem units, with f = f(x,p,t) as the
%   particle position-momentum phase space distribution,
%   gamma = sqrt( 1 + |p/(mc)|^2 ) and < a > = integral d^3 p a,
%     RHO = < q f >
%     Ji = < q pi f / (m gamma) >
%     KE = < m c^2 (gamma - 1) f >
%     Pi = < pi f >
%     Tij = < pi pj f / ( m gamma ) >
%   It is a straightforeward matter to construct the entire
%   relativistically correct stress energy tensor from these
%   quantities.
%
% Written by:
%   Kevin J. Bowers, Ph.D.
%   Plasma Physics Group (X-1)
%   Los Alamos National Lab
%   March 2004 -  Adapted from V4PIC version 1 load_species_fields.m

if nargin<2, order = [2 1 3]; end;

% Open the requested file

handle = fopen(filename);
if handle==-1, error('Could not open file'); end;

% Read binary compatibility information

cbit = fread(handle,1,'int8');
shsz = fread(handle,1,'int8');
isz  = fread(handle,1,'int8');
flsz = fread(handle,1,'int8');
dbsz = fread(handle,1,'int8');
mgcs = fread(handle,1,'uint16'); tsts = hex2dec('cafe');
mgci = fread(handle,1,'uint32'); tsti = hex2dec('deadbeef');
mgcf = fread(handle,1,'single');
mgcd = fread(handle,1,'double');

% Read the dump version and type

vers = fread(handle,1,'int32');
type = fread(handle,1,'int32');

% Read the metadata

nt   = fread(handle,1,'int32');
nx   = fread(handle,1,'int32');
ny   = fread(handle,1,'int32');
nz   = fread(handle,1,'int32');
dt   = fread(handle,1,'single');
dx   = fread(handle,1,'single');
dy   = fread(handle,1,'single');
dz   = fread(handle,1,'single');
x0   = fread(handle,1,'single');
y0   = fread(handle,1,'single');
z0   = fread(handle,1,'single');
cvac = fread(handle,1,'single');
eps0 = fread(handle,1,'single');
damp = fread(handle,1,'single');
rank = fread(handle,1,'int32');
npro = fread(handle,1,'int32');
spid = fread(handle,1,'int32');
spqm = fread(handle,1,'single');

% Read the field array header

elsz = fread(handle,1,'int32');
ndim = fread(handle,1,'int32');
dim0 = fread(handle,1,'int32');
dim1 = fread(handle,1,'int32');
dim2 = fread(handle,1,'int32');

% Check for file compatibility / corruption

if cbit~=8,      fclose(handle); error('Invalid cbit'); end
if shsz~=2,      fclose(handle); error('Invalid shsz'); end
if isz ~=4,      fclose(handle); error('Invalid isz');  end
if flsz~=4,      fclose(handle); error('Invalid flsz'); end
if dbsz~=8,      fclose(handle); error('Invalid dbsz'); end
if mgcs~=tsts,   fclose(handle); error('Invalid mgcs'); end
if mgci~=tsti,   fclose(handle); error('Invalid mgci'); end
if mgcf~=1,      fclose(handle); error('Invalid mgcf'); end
if mgcd~=1,      fclose(handle); error('Invalid mgcd'); end
if vers~=0,      fclose(handle); error('Invalid vers'); end
if type~=2,      fclose(handle); error('Invalid type'); end
if nt<0,         fclose(handle); error('Invalid nt');   end
if nx<1,         fclose(handle); error('Invalid nx');   end
if ny<1,         fclose(handle); error('Invalid ny');   end
if nz<1,         fclose(handle); error('Invalid nz');   end
if dt<0,         fclose(handle); error('Invalid dt');   end
if dx<0,         fclose(handle); error('Invalid dx');   end
if dy<0,         fclose(handle); error('Invalid dy');   end
if dz<0,         fclose(handle); error('Invalid dz');   end
if cvac<=0,      fclose(handle); error('Invalid cvac'); end
if rank<0,       fclose(handle); error('Invalid rank'); end
if rank>=npro,   fclose(handle); error('Invalid rank'); end
if npro<1,       fclose(handle); error('Invalid npro'); end
if elsz~=64,     fclose(handle); error('Invalid elsz'); end
if ndim~=3,      fclose(handle); error('Invalid ndim'); end
if dim0~=(nx+2), fclose(handle); error('Invalid dim0'); end
if dim1~=(ny+2), fclose(handle); error('Invalid dim1'); end
if dim2~=(nz+2), fclose(handle); error('Invalid dim2'); end

% Setup the grid parameters array

g = [ nt nx ny nz dt dx dy dz cvac eps0 damp x0 y0 z0 spid spqm rank npro ];

% Read the raw field data

shape = [nx+2 ny+2 nz+2];
nc = prod(shape);
data = fread(handle,[16,nc],'single');

rho = reshape( data(1, :), shape );
jx  = reshape( data(2, :), shape );
jy  = reshape( data(3, :), shape );
jz  = reshape( data(4, :), shape );
ke  = reshape( data(5, :), shape );
px  = reshape( data(6, :), shape );
py  = reshape( data(7, :), shape );
pz  = reshape( data(8, :), shape );
txx = reshape( data(9, :), shape );
tyy = reshape( data(10,:), shape );
tzz = reshape( data(11,:), shape );
tyz = reshape( data(12,:), shape );
tzx = reshape( data(13,:), shape );
txy = reshape( data(14,:), shape );

clear data
fclose(handle);

% Throw out ghost data and permute into order
rho = permute( rho(2:nx+2,2:ny+2,2:nz+2), order );
jx  = permute( jx( 2:nx+2,2:ny+2,2:nz+2), order );
jy  = permute( jy( 2:nx+2,2:ny+2,2:nz+2), order );
jz  = permute( jz( 2:nx+2,2:ny+2,2:nz+2), order );
ke  = permute( ke( 2:nx+2,2:ny+2,2:nz+2), order );
px  = permute( px( 2:nx+2,2:ny+2,2:nz+2), order );
py  = permute( py( 2:nx+2,2:ny+2,2:nz+2), order );
pz  = permute( pz( 2:nx+2,2:ny+2,2:nz+2), order );
txx = permute( txx(2:nx+2,2:ny+2,2:nz+2), order );
tyy = permute( tyy(2:nx+2,2:ny+2,2:nz+2), order );
tzz = permute( tzz(2:nx+2,2:ny+2,2:nz+2), order );
tyz = permute( tyz(2:nx+2,2:ny+2,2:nz+2), order );
tzx = permute( tzx(2:nx+2,2:ny+2,2:nz+2), order );
txy = permute( txy(2:nx+2,2:ny+2,2:nz+2), order );

