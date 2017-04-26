% run the ParaBot test case
% p517 computational acoustics

global units
units = 'km';
%%
makebty              % make the bathymetry
figure
plotbty3d ParaBot
shading flat

% ray trace
copyfile( 'ParaBot.bty', 'ParaBot_ray.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot_ray.ati' )   % copy over the bathymetry file

bellhop3d ParaBot_ray     % run BELLHOP3D on the wedge3d test case

hold on
plotray3d ParaBot_ray.ray

%%
copyfile( 'ParaBot.bty', 'ParaBot2D.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot2D.ati' )   % copy over the bathymetry file

bellhop3d ParaBot2D     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'ParaBot2D.shd' )
caxisrev( [ 60 100 ] )

%%
% 3d run (GeoHat Cartesian)

copyfile( 'ParaBot.bty', 'ParaBot3DHatcart.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DHatcart.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DHatcart     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'ParaBot3DHatcart.shd' )
caxisrev( [ 60 100 ] )

%%
% 3d run (GeoHat Ray centered)

copyfile( 'ParaBot.bty', 'ParaBot3DHatRaycen.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DHatRaycen.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DHatRaycen     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'ParaBot3DHatRaycen.shd' )
caxisrev( [ 60 100 ] )

%%
% 3d run (GeoGaussian)

copyfile( 'ParaBot.bty', 'ParaBot3DGaussian.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DGaussian.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DGaussian     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'ParaBot3DGaussian.shd' )
caxisrev( [ 60 100 ] )

%print -depsc2 ParaBot3DGaussian
%print -djpeg ParaBot3DGaussian

