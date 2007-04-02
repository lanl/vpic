function [ g, bc, range, neighbor, gid ] = ...
           load_domain_grid( filename, order );

% Usage: [ g, ...        % Grid parameters
%          bc, ...       % Grid boundary conditions
%          range, ...    % Global ID the ranges of cells owned by each rank
%          neighbor, ... % Global ID of neighbors to cells on this grid
%          gid  ...      % Global ID of cells on this grid
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
%     spid spqm      - Irrevelant for domain grid dump
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
% Written by:
%   Kevin J. Bowers, Ph.D.
%   Plasma Physics Group (X-1)
%   Los Alamos National Lab
%   March 2004 - Original version written

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
if type~=0,      fclose(handle); error('Invalid type'); end
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

% Setup the grid parameters array

g = [ nt nx ny nz dt dx dy dz cvac eps0 damp x0 y0 z0 spid spqm rank npro ];


% Read the bc header

elsz = fread(handle,1,'int32');
ndim = fread(handle,1,'int32');
dim0 = fread(handle,1,'int32');
dim1 = fread(handle,1,'int32');
dim2 = fread(handle,1,'int32');

if elsz~=4, fclose(handle); error('Invalid bc elsz'); end
if ndim~=3, fclose(handle); error('Invalid bc ndim'); end
if dim0~=3, fclose(handle); error('Invalid bc dim0'); end
if dim1~=3, fclose(handle); error('Invalid bc dim1'); end
if dim2~=3, fclose(handle); error('Invalid bc dim2'); end

bc = reshape( fread(handle,[1,dim0*dim1*dim2],'int32'), ...
              [dim0 dim1 dim2] );

% Read the range array

elsz = fread(handle,1,'int32');
ndim = fread(handle,1,'int32');
dim0 = fread(handle,1,'int32');

if elsz~=4,      fclose(handle); error('Invalid range elsz'); end
if ndim~=1,      fclose(handle); error('Invalid range ndim'); end
if dim0~=npro+1, fclose(handle); error('Invalid range dim0'); end

range = fread(handle,[1,dim0],'int32');

% Read the neighbor array

elsz = fread(handle,1,'int32');
ndim = fread(handle,1,'int32');
dim0 = fread(handle,1,'int32');
dim1 = fread(handle,1,'int32');
dim2 = fread(handle,1,'int32');
dim3 = fread(handle,1,'int32');

if elsz~=4,    fclose(handle); error('Invalid elsz'); end
if ndim~=4,    fclose(handle); error('Invalid ndim'); end
if dim0~=6,    fclose(handle); error('Invalid dim0'); end
if dim1~=nx+2, fclose(handle); error('Invalid dim1'); end
if dim2~=ny+2, fclose(handle); error('Invalid dim2'); end
if dim3~=nz+2, fclose(handle); error('Invalid dim2'); end

neighbor = reshape( fread(handle,[1,dim0*dim1*dim2*dim3],'int32'), ...
                    [dim0 dim1 dim2 dim3] );

fclose(handle);

[x,y,z] = ndgrid(0:nx+1,0:ny+1,0:nz+1);
gid = range(rank+1) + x + (nx+2)*( y + (ny+2)*z );
clear x y z

neighbor = permute(neighbor,[1 order+1]);
gid      = permute(gid,order);
