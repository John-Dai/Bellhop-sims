clear
btyfil = 'Seamount.bty';
interp_type = 'R';
interp_type = 'C';

alpha0 = 15; % 
rmax   = 15;

dxy = 0.5;
dxy = 0.5;
xmin = -rmax;  xmax = rmax;
ymin = -rmax;  ymax = rmax;

x  = ( xmin : dxy : xmax );
nx = length( x );
y  = ( ymin : dxy : ymax );
ny = length( y );

% different construction ...
nx = 150;
ny = 150;
nx = 30;
ny = 30;
x = linspace( xmin, xmax, nx );
y = linspace( ymin, ymax, ny );

z = zeros( ny, nx );
[ X, Y ] = meshgrid( x, y );
R = sqrt( X.^2 + Y.^2 );
z = R * tand( alpha0 ) * 1000;

Bathy.X = x;
Bathy.Y = y;
Bathy.depth = z;

writebty3d( btyfil, Bathy )

plotbty3d( btyfil )

