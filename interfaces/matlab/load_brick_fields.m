function [ G,                         ...
           EX, EY, EZ,                ...
           BX, BY, BZ,                ...
           JX, JY, JZ,                ...
           RHOB,                      ...
           MATERIAL,                  ...
           TCAX, TCAY, TCAZ,          ...
           RHOF, DIV_E_ERR, DIV_B_ERR ...
         ] = load_brick_fields(filebase,topo,order);

% Usage: [ g,                         ... % Grid parameters
%          ex, ey, ez,                ... % Electric fields
%          bx, by, bz,                ... % Magnetic fields
%          jfx, jfy, jfz,             ... % Free currents density
%          rhob,                      ... % Bound charge density
%          material,                  ... % Mesh of material IDs
%          tcax, tcay, tcaz,          ... % DEBUGGING
%          rhof, div_e_err, div_b_err ... % DEBUGGING
%        ] = load_brick_fields(filebase,topo,order);
%
% g - Grid parameters:
%     g = [ nt nx ny nz    ...
%           dt dx dy dz    ...
%           cvac eps0 damp ...
%           x0 y0 z0       ...
%           spid spqm      ...
%           rank npro ]
%   where:
%     nt  nx, ny, nz - Time level and grid resolution
%     dt, dx, dy, dz - Time step and grid spacing
%     cvac           - Speed of light
%     eps0           - Permittivity of free space
%     damp           - Radiation damping parameter
%     x0, y0, z0     - Offset of this domain
%     spid spqm      - Irrevelant for domain fields dump
%     rank           - -1 (all domains)
%     npro           - Total number of domains (i.e. nproc)
%
% filebase - Base name of the fields dump file to load.
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
% - rhob may not be in sync on domain boundaries when this is
%   called.  It is doubtful this will cause any problems but
%   if so, let KJB know and this script can be updated.
% - tcax, tcay, tcaz, rhof, div_e_err and div_b_err are for
%   code debugging purposes only. They may not be current values
%   or in sync from node to node.  (They are used for debugging!)
% - You should not use this routine unless you know you have
%   enough physical RAM to hold all the fields in memory
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
        ex, ey, ez, ...
        bx, by, bz, ...
        jx, jy, jz, ...
        rhob, ...
        material, ...
        tcax, tcay, tcaz, ...
        rhof, div_e_err, div_b_err ...
      ] = load_domain_fields( filename, [1 2 3] );

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
        EX        = zeros(NX,  NY+1,NZ+1);
        EY        = zeros(NX+1,NY,  NZ+1);
        EZ        = zeros(NX+1,NY+1,NZ  );

        BX        = zeros(NX,  NY+1,NZ+1);
        BY        = zeros(NX+1,NY,  NZ+1);
        BZ        = zeros(NX+1,NY+1,NZ  );

        JX        = zeros(NX,  NY+1,NZ+1);
        JY        = zeros(NX+1,NY,  NZ+1);
        JZ        = zeros(NX+1,NY+1,NZ  );

        RHOB      = zeros(NX+1,NY+1,NZ+1);

        MATERIAL  = zeros(2*NX+1,2*NY+1,2*NZ+1);

        TCAX      = zeros(NX,  NY+1,NZ+1);
        TCAY      = zeros(NX+1,NY,  NZ+1);
        TCAZ      = zeros(NX+1,NY+1,NZ  );

        RHOF      = zeros(NX+1,NY+1,NZ+1);
        DIV_E_ERR = zeros(NX+1,NY+1,NZ+1);
        DIV_B_ERR = zeros(NX,  NY,  NZ  );

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

      IIs = 1 + ii*nx; IIe = (ii+1)*nx;
      JJs = 1 + jj*ny; JJe = (jj+1)*ny;
      KKs = 1 + kk*nz; KKe = (kk+1)*nz;

      IIms = 1 + ii*2*nx; IIme = 1 + (ii+1)*2*nx;
      JJms = 1 + jj*2*ny; JJme = 1 + (jj+1)*2*ny;
      KKms = 1 + kk*2*nz; KKme = 1 + (kk+1)*2*nz;

      EX(       IIs:IIe,  JJs:JJe+1,KKs:KKe+1) = ex;
      EY(       IIs:IIe+1,JJs:JJe,  KKs:KKe+1) = ey;
      EZ(       IIs:IIe+1,JJs:JJe+1,KKs:KKe  ) = ez;
      BX(       IIs:IIe+1,JJs:JJe,  KKs:KKe  ) = bx;
      BY(       IIs:IIe,  JJs:JJe+1,KKs:KKe  ) = by;
      BZ(       IIs:IIe,  JJs:JJe,  KKs:KKe+1) = bz;
      JX(       IIs:IIe,  JJs:JJe+1,KKs:KKe+1) = jx;
      JY(       IIs:IIe+1,JJs:JJe,  KKs:KKe+1) = jy;
      JZ(       IIs:IIe+1,JJs:JJe+1,KKs:KKe  ) = jz;
      RHOB(     IIs:IIe+1,JJs:JJe+1,KKs:KKe+1) = rhob;
      MATERIAL( IIms:IIme,JJms:JJme,KKms:KKme) = material;
      TCAX(     IIs:IIe,  JJs:JJe+1,KKs:KKe+1) = tcax;
      TCAY(     IIs:IIe+1,JJs:JJe,  KKs:KKe+1) = tcay;
      TCAZ(     IIs:IIe+1,JJs:JJe+1,KKs:KKe  ) = tcaz;
      RHOF(     IIs:IIe+1,JJs:JJe+1,KKs:KKe+1) = rhof;
      DIV_E_ERR(IIs:IIe+1,JJs:JJe+1,KKs:KKe+1) = div_e_err;
      DIV_B_ERR(IIs:IIe,  JJs:JJe,  KKs:KKe  ) = div_b_err;

    end
  end
end

if any(order~=[1 2 3]),

  EX        = permute( EX,        order );
  EY        = permute( EY,        order );
  EZ        = permute( EZ,        order );
  BX        = permute( BX,        order );
  BY        = permute( BY,        order );
  BZ        = permute( BZ,        order );
  JX        = permute( JX,        order );
  JY        = permute( JY,        order );
  JZ        = permute( JZ,        order );
  RHOB      = permute( RHOB,      order );
  MATERIAL  = permute( MATERIAL,  order );
  TCAX      = permute( TCAX,      order );
  TCAY      = permute( TCAY,      order );
  TCAZ      = permute( TCAZ,      order );
  RHOF      = permute( RHOF,      order );
  DIV_E_ERR = permute( DIV_E_ERR, order );
  DIV_B_ERR = permute( DIV_B_ERR, order );

end
