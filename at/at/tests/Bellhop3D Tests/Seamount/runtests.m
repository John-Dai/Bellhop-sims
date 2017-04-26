% run the seamount test case
% p517 computational acoustics
clear

global units
units = 'km';
%%
makebty              % make the bathymetry
figure
plotbty3d Seamount
shading flat

%% ray trace
copyfile( 'Seamount.bty', 'Seamount_ray.bty' )   % copy over the bathymetry file
bellhop3d Seamount_ray     % run BELLHOP3D on the wedge3d test case

hold on
plotray3d Seamount_ray.ray

%%
copyfile( 'Seamount.bty', 'Seamount2D.bty' )   % copy over the bathymetry file

bellhop3d Seamount2D     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'Seamount2D.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%%
% 3d run (GeoHat Cartesian)

copyfile( 'Seamount.bty', 'Seamount3DHatcart.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DHatcart     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'Seamount3DHatcart.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%%
% 3d run (GeoHat Ray centered)

copyfile( 'Seamount.bty', 'Seamount3DHatRaycen.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DHatRaycen     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'Seamount3DHatRaycen.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%%
% 3d run (GeoGaussian)

copyfile( 'Seamount.bty', 'Seamount3DGaussian.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DGaussian     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'Seamount3DGaussian.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%print -depsc2 Seamount3DGaussian
%print -djpeg Seamount3DGaussian

