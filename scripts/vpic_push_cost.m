function flop = vpic_push_cost( vt, vd, cell, dt )

% Compute the expected number of flops used to push a particle from a
% given distribution with given simulation parameters in VPIC.
%
% A random distribution of particles is created and pushed one step.  The
% resulting spatial distribution is then analyzed to compute the probability a
% particle will cross a given number of boundaries.  Results are generally
% accurate to 3 sig figs (increase N if you want more).
%
% Flop model: 246 flop + 105 flop per boundary crossed
%
% Determined by op counting the inner loop.
%
% Particles are assumed uniformly distributed in space.
% Particles are assumed to have a drifting Maxwellian distribution.
%
% vt is the particle thermal velocity (in LENGTH/TIME).  This is generally
% a scalar but if you wish to model anisotropic distributions you can pass
% a 3x3 matrix here as well.
%
% vd is the particle drift velocity (in LENGTH/TIME).
% Use a 1x3 row vector here.
%
% cell gives the cell dimensions in length (in LENGTH).  Use either a scalar or
% a % 1x3 row vector here.
%
% dt gives the timestep.  In TIME. Use a scalar here.
%
% Distribution not completely accurate for relativistic temperatures

N = 10000000;
x = rand( N, 3 ) + ( repmat( vd, N, 1 ) + randn( N, 3 )*vt' )*diag( dt./cell );
n = sum( x<0 | x>1, 2 );
i = 0:max(n);
flop = 246 + 105*sum( i.*hist( n, i, 1 ) );
