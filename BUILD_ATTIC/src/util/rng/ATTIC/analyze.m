function analyze( file, type, range, integer, fno )
  disp( [ 'Analyzing ' file ] );

  f     = fopen( file, 'rb' );
  [X,N] = fread( f, type );
  fclose( file );

  % Sanity check the data

  X0 = min(X);
  X1 = max(X);
  if X0<range(1) || X1>range(2),
    disp( 'Out of range values detected' );
  end
  if integer,
    if any(round(X)!=X ),
      disp( 'Non-integral values detected' );
    end
  end

  % Estimate the density and distribution functions 

  M = floor( N^(1/3) );
  if integer,
    MI = range(2) - range(1) + 1;
    if MI<2048,
      M = MI;
      X0 = X0 - 0.5;
      X1 = X1 + 0.5;
    end
  end
  x_ = linspace( X0, X1, M+1 );
  x  = 0.5*( x_(2:M+1) + x_(1:M) ); % Density function at x
  x_ = x_(2:M+1);                   % Distrbution function at x_
  dx = (X1-X0)/M;
  f  = hist( X, x ) / (N*dx);
  F_ = cumsum( f ) * dx;

  % Octave doesn't handle subplots on multiple figures right
  figure( 2*fno-1 );
  plot( x, f );
  axis( [ X0, X1, 0, max(f)*1.05 ] );
  title( [ file ' density function' ] );
  xlabel( 'x' );
  ylabel( 'f_X(x)' );

  figure( 2*fno );
  plot( x_, F_ );
  axis( [ X0, X1, 0, 1 ] );
  title( [ file ' distribution function' ] );
  xlabel( 'x' );
  ylabel( 'F_X(x)' );

