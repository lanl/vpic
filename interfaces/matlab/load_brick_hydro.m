function [ G, ...
           RHO, JX, JY, JZ, ...
           KE, PX, PY, PZ, ...
           TXX, TYY, TZZ, ...
           TYZ, TZX, TXY ...
         ] = load_brick_hydro(filebase,topo,order);

% Usage: [ g, ...             % Grid and species parameters
%          RHO, ...           % Charge density
%          Jx, Jy, Jz, ...    % Current density (charge flux)
%          K, ...             % Kinetic energy density
%          Px, Py, Pz, ...    % Momentum density (energy flux)
%          Txx, Tyy, Tzz, ... % Stress tensor diagonal
%          Tyz, Tzx, Txy ...  % Stress tensor off-diagonal
%        ] = load_brick_hydro(filebase,topo,order);
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
%     rank           - -1 (all domains)
%     ndom           - Total number of domains (i.e. nproc)
%
% filename - Name of the hydro dump file to load.
%
% topo - Topology of the brick (processors along each edge)
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
% - You should not use this routine unless you know you have
%   enough physical RAM to hold all the hydro fields in memory
%   at once.
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
      [ g, ...
        rho,jx,jy,jz, ...
        ke,px,py,pz, ...
        txx,tyy,tzz, ...
        tyz,tzx,txy ...
      ] = load_domain_hydro(filename,[1 2 3]);

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

        % Allocate outputs
        RHO = zeros(NX+1,NY+1,NZ+1);
        JX  = zeros(NX+1,NY+1,NZ+1);
        JY  = zeros(NX+1,NY+1,NZ+1);
        JZ  = zeros(NX+1,NY+1,NZ+1);

        KE  = zeros(NX+1,NY+1,NZ+1);
        PX  = zeros(NX+1,NY+1,NZ+1);
        PY  = zeros(NX+1,NY+1,NZ+1);
        PZ  = zeros(NX+1,NY+1,NZ+1);

        TXX = zeros(NX+1,NY+1,NZ+1);
        TYY = zeros(NX+1,NY+1,NZ+1);
        TZZ = zeros(NX+1,NY+1,NZ+1);

        TYZ = zeros(NX+1,NY+1,NZ+1);
        TZX = zeros(NX+1,NY+1,NZ+1);
        TXY = zeros(NX+1,NY+1,NZ+1);

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

      IIs = 1 + ii*nx; IIe = 1 + (ii+1)*nx;
      JJs = 1 + jj*ny; JJe = 1 + (jj+1)*ny;
      KKs = 1 + kk*nz; KKe = 1 + (kk+1)*nz;
      
      RHO(IIs:IIe,JJs:JJe,KKs:KKe) = rho;
      JX( IIs:IIe,JJs:JJe,KKs:KKe) = jx;
      JY( IIs:IIe,JJs:JJe,KKs:KKe) = jy;
      JZ( IIs:IIe,JJs:JJe,KKs:KKe) = jz;

      KE( IIs:IIe,JJs:JJe,KKs:KKe) = ke;
      PX( IIs:IIe,JJs:JJe,KKs:KKe) = px;
      PY( IIs:IIe,JJs:JJe,KKs:KKe) = py;
      PZ( IIs:IIe,JJs:JJe,KKs:KKe) = pz;

      TXX(IIs:IIe,JJs:JJe,KKs:KKe) = txx;
      TYY(IIs:IIe,JJs:JJe,KKs:KKe) = tyy;
      TZZ(IIs:IIe,JJs:JJe,KKs:KKe) = tzz;

      TYZ(IIs:IIe,JJs:JJe,KKs:KKe) = tyz;
      TZX(IIs:IIe,JJs:JJe,KKs:KKe) = tzx;
      TXY(IIs:IIe,JJs:JJe,KKs:KKe) = txy;


    end
  end
end

if any(order~=[1 2 3]),

  RHO = permute( RHO, order );
  JX  = permute( JX,  order );
  JY  = permute( JY,  order );
  JZ  = permute( JZ,  order );

  KE  = permute( KE,  order );
  PX  = permute( PX,  order );
  PY  = permute( PY,  order );
  PZ  = permute( PZ,  order );

  TXX = permute( TXX, order );
  TYY = permute( TYY, order );
  TZZ = permute( TZZ, order );

  TYZ = permute( TYZ, order );
  TZX = permute( TZX, order );
  TXY = permute( TXY, order );

end
