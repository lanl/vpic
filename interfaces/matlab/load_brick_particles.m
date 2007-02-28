function [ G, PX, PY, PZ, PUX, PUY, PUZ, PQ ] = ...
           load_brick_particles(filebase,topo);

% Usage: [ g, px, py, pz, pux, puy, puz, pq ] = ...
%          load_brick_particles(filebase,topo);
%
% g - Grid parameters:
%     g = [ nt nx ny nz ...
%           dt dx dy dz ...
%           cvac eps0 damp ...
%           x0 y0 z0 ...
%           spid spqm ...
%           rank ndom ]
%
% topo - Topology of the brick (processors along each edge)
%
%   where:
%     nt  nx, ny, nz - Time level and grid resolution
%     dt, dx, dy, dz - Time step and grid spacing
%     cvac           - Speed of light
%     eps0           - Permittivity of free space
%     damp           - Radiation damping parameter
%     x0, y0, z0     - Offset of this domain
%     spid spqm      - Species ID and the charge to mass ratio
%     rank           - -1 (all domains)
%     ndom           - Total number of domains (i.e. nproc)
%
% filename - Name of the particle dump file to load.
%
% Notes:
% - All values are time centered appropriately and in the
%   appropriate units.
% - You should not use this routine unless you know you have
%   enough physical RAM to hold all the particles in memory at
%   once.
%
% Written by:
%   Kevin J. Bowers, Ph.D.
%   Plasma Physics Group (X-1)
%   Los Alamos National Lab
%   April 2004 - Original version

if length(topo)~=3, error( 'Bad topo' ); end
if nargin<3, order = [ 2 1 3 ]; end

for kk=0:topo(3)-1,
  for jj=0:topo(2)-1,
    for ii=0:topo(1)-1,

      rr = ii + topo(1)*( jj + topo(2)*kk );
      filename = sprintf('%s.%i',filebase,rr);
      [ g, px, py, pz, pux, puy, puz, pq ] = ...
        load_domain_particles( filename );

      nt = g(1); nx = g(2); ny = g(3); nz = g(4);
      dt = g(5); dx = g(6); dy = g(7); dz = g(8);
      cvac = g(9); eps0 = g(10); damp = g(11);
      x0 = g(12); y0 = g(13); z0 = g(14);
      spid = g(15); spqm = g(16);
      rank = g(17); npro = g(18);

      if rank~=rr, error( 'Inconsistent rank' ); end

      if rr==0,

        % First time, construct the output grid info
        NX = nx*topo(1);
        NY = ny*topo(2);
        NZ = nz*topo(3);
        G = [ nt nx*topo(1) ny*topo(2) nz*topo(3) ...
              dt dx dy dz ...
              cvac eps0 damp ...
              x0 y0, z0 ...
              spid spqm ...
              -1 npro ];

        PX  = [];
        PY  = [];
        PZ  = [];
        PUX = [];
        PUY = [];
        PUZ = [];
        PQ  = [];

      else
        % Make sure current domain is consistent
        if nt~=G(1),         error( 'Inconsistent nt'   ); end
        if nx~=G(2)/topo(1), error( 'Inconsistent nx'   ); end
        if ny~=G(3)/topo(2), error( 'Inconsistent ny'   ); end
        if nz~=G(4)/topo(3), error( 'Inconsistent nz'   ); end
        if dt~=G(5),         error( 'Inconsistent dt'   ); end
        if dx~=G(6),         error( 'Inconsistent dx'   ); end
        if dy~=G(7),         error( 'Inconsistent dy'   ); end
        if dz~=G(8),         error( 'Inconsistent dz'   ); end
        if cvac~=G(9),       error( 'Inconsistent cvac' ); end
        if eps0~=G(10),      error( 'Inconsistent eps0' ); end
        if damp~=G(11),      error( 'Inconsistent damp' ); end
        % FIXME: Add robust checks for x0, y0, z0 values
        if spid~=G(15),      error( 'Inconsistent spid' ); end
        if spqm~=G(16),      error( 'Inconsistent spqm' ); end
        % G(17) is -1 to indicate assembled values
        if npro~=G(18),      error( 'Inconsistent cvac' ); end
      end

      % Drop in the fields from the current domain
      PX  = [ PX,  px  ];
      PY  = [ PY,  py  ];
      PZ  = [ PZ,  pz  ];
      PUX = [ PUX, pux ];
      PUY = [ PUY, puy ];
      PUZ = [ PUZ, puz ];
      PQ  = [ PQ,  pq  ];

    end
  end
end
