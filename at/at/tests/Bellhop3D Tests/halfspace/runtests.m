% Run tests to verify the surface and bottom reflection coefficents are
% incorporated correctly

% free space

global units
units = 'km';

scooter( 'freeSPoint' )
plotshd( 'freeSPoint.shd.mat', 4, 1, 1 );
caxisrev( [ 60 80 ] )

bellhop3d( 'freeBhat' )
plotshd( 'freeBhat.shd', 4, 1, 2 );
caxisrev( [ 60 80 ] )

bellhop3d( 'freeBhat_raycen' )
plotshd( 'freeBhat_raycen.shd', 4, 1, 3 );
caxisrev( [ 60 80 ] )

bellhop3d( 'freeBgaussian' )
plotshd( 'freeBgaussian.shd', 4, 1, 4 );
caxisrev( [ 60 80 ] )

%%
% polar plots

bellhop3d( 'freeBhatpolar' )
figure
plotshdpol( 'freeBhatpolar.shd', 0.0, 0.0, 3000 );
caxisrev( [ 60 80 ] )

%%

bellhop3d( 'freeBgaussianpolar' )
figure
plotshdpol( 'freeBgaussianpolar.shd', 0.0, 0.0, 3000 );
caxisrev( [ 60 80 ] )
%%
% halfspaces

scooter( 'lower_halfS' )
plotshd( 'lower_halfS.shd.mat', 2, 2, 1 );
caxisrev( [ 60 80 ] )

bellhop3d( 'lower_halfB' )
plotshd( 'lower_halfB.shd', 2, 2, 2 );
caxisrev( [ 60 80 ] )

scooter( 'upper_halfS' )
plotshd( 'upper_halfS.shd.mat', 2, 2, 3 );
caxisrev( [ 60 80 ] )

bellhop3d( 'upper_halfB' )
plotshd( 'upper_halfB.shd', 2, 2, 4 );
caxisrev( [ 60 80 ] )
