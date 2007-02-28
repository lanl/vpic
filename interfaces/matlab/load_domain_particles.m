function [ g, px, py, pz, pux, puy, puz, pq ] = ...
           load_domain_particles(filename);

% Usage: [ g, px, py, pz, pux, puy, puz, pq ] = ...
%          load_domain_particles(filename);
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
% filename - Name of the particle dump file to load.
%
% Notes:
% - All values are time centered appropriately and in the appropriate units.
%
% Written by:
%   Kevin J. Bowers, Ph.D.
%   Plasma Physics Group (X-1)
%   Los Alamos National Lab
%   March 2004 -  Adapted from V4PIC version 1 load_species.m 

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
if type~=3,      fclose(handle); error('Invalid type'); end
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

% Setup the grid parameters array

g = [ nt nx ny nz dt dx dy dz cvac eps0 damp x0 y0 z0 spid spqm rank npro ];

% Read the raw particle data

elsz = fread(handle,1,'int32');
ndim = fread(handle,1,'int32');
dim0 = fread(handle,1,'int32');

if elsz~=32, fclose(handle); error('Invalid elsz'); end
if ndim~=1,  fclose(handle); error('Invalid ndim'); end
if dim0<0,   fclose(handle); error('Invalid dim0'); end

if dim0==0,
  px=[]; py=[]; pz=[]; pux=[]; puy=[]; puz=[]; pq=[];
  fclose(handle);
  return;
end

data_start = ftell(handle);
data = fread(handle,[8,dim0],'single');
pox = data(1,:);
poy = data(2,:);
poz = data(3,:);
pux = data(5,:);
puy = data(6,:);
puz = data(7,:);
pq  = data(8,:);
clear data

fseek(handle,data_start,'bof');
fread(handle,3,'single');
pix = fread(handle,[1,dim0],'int32',28);

fclose(handle);

piy = floor(pix/(nx+2));
pix = pix - piy*(nx+2);
piz = floor(piy/(ny+2));
piy = piy - piz*(ny+2);

px = x0 + (pix+0.5*pox-0.5)*dx;
py = y0 + (piy+0.5*poy-0.5)*dy;
pz = z0 + (piz+0.5*poz-0.5)*dz;


