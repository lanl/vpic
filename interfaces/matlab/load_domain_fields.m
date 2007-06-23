function [ g,                         ...
           ex, ey, ez,                ...
           bx, by, bz,                ...
           jfx, jfy, jfz,             ...
           rhob,                      ...
           material,                  ...
           tcax, tcay, tcaz,          ...
           rhof, div_e_err, div_b_err ...
         ] = load_domain_fields(filename,order);

% Usage: [ g, ...                         % Grid parameters
%          ex, ey, ez, ...                % Electric fields
%          bx, by, bz, ...                % Magnetic fields
%          jfx, jfy, jfz, ...             % Free currents density
%          rhob, ...                      % Bound charge density
%          material, ...                  % Mesh of material IDs
%          tcax, tcay, tcaz,              % DEBUGGING
%          rhof, div_e_err, div_b_err ... % DEBUGGING
%        ] = load_domain_fields(filename,order);
%
% g - Grid parameters:
%     g = [ nt nx ny nz ...
%           dt dx dy dz ...
%           cvac eps0 damp
%           x0 y0 z0
%           spid spqm
%           rank npro ]
%   where:
%     nt  nx, ny, nz - Time level and grid resolution
%     dt, dx, dy, dz - Time step and grid spacing
%     cvac           - Speed of light
%     eps0           - Permittivity of free space
%     damp           - Radiation damping parameter
%     x0, y0, z0     - Offset of this domain
%     spid spqm      - Irrevelant for domain fields dump
%     rank           - ID of this domain
%     npro           - Total number of domains (i.e. nproc)
%
% filename - Name of the fields dump file to load.
%
% order - (optional) Indicates the desired field indexing.
%   [2 1 3] (default) - YXZ indexing (compatible with "meshgrid"
%                       and MATLAB's 3d plotting routines)
%   [1 2 3] - XYZ indexing (compatible with "ndgrid")
%   Other orderings are easy to figure out.
%
% Notes:
% - All fields are time centered appropriately except jfx, jfy, jfz.
%   and divj. These are known half a time step back from the
% - Returned fields are located on the Yee mesh.
%   Ghost values are not included in the returned fields.
%   The materials mesh is a node centered mesh with a spacing of
%   dx/2,dy/2,dz/2
% - All fields are in the appropriate units for the problem.
%   To understand how the units of the problem are defined,
%   the following equations govern the fields:
%     partial_t B =       -curl E            ... Faraday's Law
%     partial_t E = cvac^2 curl B - Jf/eps0  ... Ampere's Law
%     div E = (rhof+rhob)/eps0               ... E boundary cond.
%     div B = 0                              ... B boundary cond.
% - tcax, tcay, tcaz, rhof, div_e_err and div_b_err are for
%   code debugging purposes only. They may not be current values!
%   (They are used for code debugging!)
%
% Written by:
%   Kevin J. Bowers, Ph.D.
%   Plasma Physics Group (X-1)
%   Los Alamos National Lab
%   March 2004 -  Adapted from V4PIC version 1 load_fields.m 

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
if type~=1,      fclose(handle); error('Invalid type'); end
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
if spid~=-1,     fclose(handle); error('Invalid spid'); end
if spqm~=0,      fclose(handle); error('Invalid spqm'); end
if elsz~=80,     fclose(handle); error('Invalid elsz'); end
if ndim~=3,      fclose(handle); error('Invalid ndim'); end
if dim0~=(nx+2), fclose(handle); error('Invalid dim0'); end
if dim1~=(ny+2), fclose(handle); error('Invalid dim1'); end
if dim2~=(nz+2), fclose(handle); error('Invalid dim2'); end

% Setup the grid parameters array

g = [ nt nx ny nz dt dx dy dz cvac eps0 damp x0 y0 z0 spid spqm rank npro ];

% Read the raw field data

data_start = ftell(handle);
shape = [nx+2 ny+2 nz+2];
nc = prod(shape);

data = fread(handle,[20,nc],'single');
ex        = reshape( data(1, :), shape );
ey        = reshape( data(2, :), shape );
ez        = reshape( data(3, :), shape );
div_e_err = reshape( data(4, :), shape );
bx        = reshape( data(5, :), shape ) / cvac;
by        = reshape( data(6, :), shape ) / cvac;
bz        = reshape( data(7, :), shape ) / cvac;
div_b_err = reshape( data(8, :), shape ) / cvac;
tcax      = reshape( data(9, :), shape );
tcay      = reshape( data(10,:), shape );
tcaz      = reshape( data(11,:), shape );
rhob      = reshape( data(12,:), shape );
jfx       = reshape( data(13,:), shape );
jfy       = reshape( data(14,:), shape );
jfz       = reshape( data(15,:), shape );
rhof      = reshape( data(16,:), shape );

% Handle material data

fseek(handle,data_start,'bof'); % Seek back to field data
fread(handle,16,'single');      % Line up with material data

data = fread(handle,[8,nc],'8*uint16',64);
ematx = reshape( data(1,:), shape );
ematy = reshape( data(2,:), shape );
ematz = reshape( data(3,:), shape );
nmat  = reshape( data(4,:), shape );
fmatx = reshape( data(5,:), shape );
fmaty = reshape( data(6,:), shape );
fmatz = reshape( data(7,:), shape );
cmat  = reshape( data(8,:), shape );

clear data
fclose(handle);

% Throw out ghost data for edge mesh quantities
ex    = permute( ex(  2:nx+1,2:ny+2,2:nz+2), order );
ey    = permute( ey(  2:nx+2,2:ny+1,2:nz+2), order );
ez    = permute( ez(  2:nx+2,2:ny+2,2:nz+1), order );
jfx   = permute( jfx( 2:nx+1,2:ny+2,2:nz+2), order );
jfy   = permute( jfy( 2:nx+2,2:ny+1,2:nz+2), order );
jfz   = permute( jfz( 2:nx+2,2:ny+2,2:nz+1), order );
tcax  = permute( tcax(2:nx+1,2:ny+2,2:nz+2), order );
tcay  = permute( tcay(2:nx+2,2:ny+1,2:nz+2), order );
tcaz  = permute( tcaz(2:nx+2,2:ny+2,2:nz+1), order );

% Throw out ghost data for face mesh quantities
bx    = permute( bx(2:nx+2,2:ny+1,2:nz+1), order );
by    = permute( by(2:nx+1,2:ny+2,2:nz+1), order );
bz    = permute( bz(2:nx+1,2:ny+1,2:nz+2), order );

% Throw out ghost data for node mesh quantities
rhof      = permute( rhof(     2:nx+2,2:ny+2,2:nz+2), order );
rhob      = permute( rhob(     2:nx+2,2:ny+2,2:nz+2), order );
div_e_err = permute( div_e_err(2:nx+2,2:ny+2,2:nz+2), order );

% Throw out ghost data for cell mesh quantities
div_b_err = permute( div_b_err(2:nx+1,2:ny+1,2:nz+1), order );

% Put material ID data into a unified array

ematx = ematx(2:nx+1,2:ny+2,2:nz+2);
ematy = ematy(2:nx+2,2:ny+1,2:nz+2);
ematz = ematz(2:nx+2,2:ny+2,2:nz+1);
nmat  = nmat(2:nx+2,2:ny+2,2:nz+2);
fmatx = fmatx(2:nx+2,2:ny+1,2:nz+1);
fmaty = fmaty(2:nx+1,2:ny+2,2:nz+1);
fmatz = fmatz(2:nx+1,2:ny+1,2:nz+2);
cmat  = cmat(2:nx+1,2:ny+1,2:nz+1);

% FIXME: PREALLOCATE material??
material(1:2:2*nx+1,1:2:2*ny+1,1:2:2*nz+1) = nmat;  clear nmat
material(2:2:2*nx,  1:2:2*ny+1,1:2:2*nz+1) = ematx; clear ematx
material(1:2:2*nx+1,2:2:2*ny,1:2:2*nz+1  ) = ematy; clear ematy
material(1:2:2*nx+1,1:2:2*ny+1,2:2:2*nz  ) = ematz; clear ematz
material(1:2:2*nx+1,2:2:2*ny,  2:2:2*nz  ) = fmatx; clear fmatx
material(2:2:2*nx,  1:2:2*ny+1,2:2:2*nz  ) = fmaty; clear fmaty
material(2:2:2*nx,  2:2:2*ny,  1:2:2*nz+1) = fmatz; clear fmatz
material(2:2:2*nx,  2:2:2*ny,  2:2:2*nz  ) = cmat;  clear cmat

material = permute( material, order );
