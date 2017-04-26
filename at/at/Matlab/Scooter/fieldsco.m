function fieldsco( filename, PlotTitle, freq, atten, Pos, Gtemp, Opt, Rminkm, Rmaxkm, Nrr )

% Calculates the field using the Green's function produced by SCOOTER
% This version uses the trapezoidal rule directly to do a DFT, rather than an FFT
%
% usage:
%    fieldsco( filename )
%  or:
%    fieldsco( filename, PlotTitle, freq, atten, Pos, Gtemp, option, Rminkm, Rmaxkm, Nrr )
%
% You must include the file extension, if it exists
% Optionlly reads control info from a file fields.flp
% mbp, Dec. 2012 based on fieldsco

% phase speed limits for optional tapering
% this is for user play (at your own risk)
global cmin cmax

if ( isempty( cmin ) )
   cmin = 1e-10;
   cmax = 1e30;
end

% if user is just supplying a file name then read all the variables

if ( strcmp( filename( end - 3 : end ), '.mat' ) )
   fileroot = filename( 1 : end - 8 );
else
   fileroot = filename( 1 : end - 4 );
end

% disp( [ 'fieldsco ' filename ] )

if ( nargin == 1 )
   % read fields.flp data to select range limits
   fid = fopen( [ fileroot '.flp' ], 'r' );
   if ( fid < 3 )
      error( 'flp file missing' )
   end
   Opt = fgetl( fid );   % x or r coordinate; positive/negative/both sides of spectrum
   % Extract letters between the quotes
   nchars = strfind( Opt, '''' );   % find quotes
   Opt    = Opt( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
   if ( length( Opt ) <= 3 )
      Opt( 4 : 4 ) = 'O';   % default beampattern is omni
   end

   Rminkm = fscanf( fid, '%f', 1 );
   Rmaxkm = fscanf( fid, '%f', 1 );
   Nrr    = fscanf( fid, '%i', 1 );
   fclose( fid );
   
   % read in the Green's function
   [ PlotTitle, ~, freq, atten, Pos, Gtemp ] = read_shd( filename );
   k = Pos.r.range;
   deltak = ( k( end ) - k( 1 ) ) / length( k );
end

% optionally read in a source beam pattern
SBP = Opt( 4 : 4 );
SrcBmPat = readpat( 'shaded', SBP );

Nsd = length( Pos.s.depth );    % # of source depths
Nrd = length( Pos.r.depth );    % # of source depths
Nk  = length( k );              % # of wavenumbers

fprintf( '\n\n--- FieldSCO --- \n\nRminkm = %d, Rmaxkm = %d, Nrr = %i \n', Rminkm, Rmaxkm, Nrr )
if ( Rminkm < 0 )
   %warning( 'Rmin <0, G(-r) defined as G(r) if r<0' )
   error( 'Rmin must be nonnegative' )
end
Rmax    = 1000 * Rmaxkm;
Rmin    = 1000 * Rminkm;
rr      = linspace( Rmin, Rmax, Nrr );

ck = k + 1i * atten;
x  = ck * rr;

switch Opt( 1 : 1 )
   case 'X'
      X  = exp( -1i * x );   % e^( -i k r ) matrix
      X2 = exp( +1i * x );   % e^( +i k r ) matrix
      
      factor1 = 1;
      factor2 = deltak ./ sqrt( 2 * pi );
   case 'R'
      X  = exp( -1i * ( x - pi / 4 ) );   % e^( -i k r ) matrix
      X2 = exp( +1i * ( x - pi / 4 ) );   % e^( +i k r ) matrix
      
      factor1 = sqrt( ck );
      factor2 = deltak ./ sqrt( 2 * pi * rr );
end

omega  = 2.0 * pi * freq;
kleft  = omega / cmax; % left  limit for tapering
kright = omega / cmin; % right limit for tapering

pressure = zeros( 1, Nsd, Nrd, Nrr );

for isd = 1: Nsd
   G = squeeze( Gtemp( 1, isd, :, : ) );
   if size( G, 2 ) == 1 % if G is a vector, expand it into a matrix with one row
      G = reshape( G, 1, length( G ) );
   end
   
   % apply the source beam pattern
   if ( SBP == '*' )
      c = 1500;   % reference sound speed, should be speed at the source depth
      kz2 = omega^2 / c^2 - k.^2;   % vertical wavenumber squared
      kz2( kz2 < 0 ) = 0;           % remove negative values
      
      theta = atand( sqrt( kz2 ) ./ k );   % calculate the angle in degrees
      S = interp1( SrcBmPat( :, 1 ), SrcBmPat( : , 2 ), theta );
      G = scalecol( G, S );         % apply the shading
   end
   
   G = taper( G, k, Nk, kleft, kright );
   G( abs( G ) < 1e-40 ) = 0;   % avoid underflows--- they slow the following down a huge amount
   
   switch Opt( 2: 2 )
      case 'P'
         G = scalecol( G, factor1 );
         Y = -G * X;   % here's where the DFT is done
      case 'N'   % check this !!!
         G = scalecol( G, factor1 );
         Y = -G * X2;
      case 'B'   % check this !!!
         G = scalecol( G, factor1 );
         Y = -G * ( X + X2 );
         %Y = -G * cos(( x - pi / 4 ) );
   end
   
   pressure( 1, isd, :, : ) = scalecol( Y, factor2 );  % cylindrical spreading, etc.
   
   fprintf( 'Transform completed for source depth %f \n', Pos.s.depth( isd ) );
end

Pos.r.range = rr;
atten       = 0;
PlotType    = 'rectilin  ';

save( [ fileroot '.shd.mat' ], 'PlotTitle', 'PlotType', 'freq', 'atten', 'Pos', 'pressure' )

end

%%

function G = taper( G, k, Nk, kleft, kright )

% windowing to smooth out any discontinuities at the end of the spectrum

if ( kleft > k( 1 )  )
   Nwinleft = 2 * round( ( kleft - k( 1 ) ) / ( k(end) - k( 1 ) ) * Nk ) + 1;  % odd number covering 10% of the spectrum
   winleft  = hanning( Nwinleft )';
   Nwinleft = ( Nwinleft - 1 ) / 2;
   winleft = winleft( 1: Nwinleft ); % taking just the left half of the symmetric Hanning window
else
   Nwinleft = 0;
   winleft  = [];
end

if ( kright < k( end )  )
   Nwinright = 2 * round( ( k(end) - kright ) / ( k(end) - k( 1 ) ) * Nk ) + 1; % odd number covering 10% of the spectrum
   winright  = hanning( Nwinright )';
   Nwinright = ( Nwinright - 1 ) / 2;
   winright  = winright( end - Nwinright + 1: end ); % taking just the right half of the symmetric Hanning window
else
   Nwinright = 0;
   winright  = [];
end

if ( Nk < Nwinleft + Nwinright )
   error( [ mfilename, 'phase speed limits for windowing are invalid' ] );
end

window  = [ winleft ones( 1, Nk - Nwinleft - Nwinright ) winright ];

G = scalecol( G, window );
end % of taper
