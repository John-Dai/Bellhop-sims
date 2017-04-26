
% compares noise from the modal and spectral integral formulas

% run KRAKEN; done this way since we only want the modes, not the field
filename = 'Pekeris';
runkraken = which( 'kraken.exe' );

if ( isempty( runkraken ) )
   error( 'kraken.exe not found in your Matlab path' )
else
   eval( [ '! "' runkraken '" ' filename ] );
end

% plot_noise

Rmax_km = linspace( 1, 100, 50 );   % max range of integration of noise sources

ModeFile = 'pekeris.mod';

sd = 1.25;  % depth of noise sources
rd = 0:1:100;

SL       = 168.0;   % source level in dB
rho_ship = 3;       % ships per square degree

% following assumes a degree at the equator ...
km_per_degree = 107;   % (is this number right?)
rho_SL = 10^( SL / 10 ) * rho_ship / ( km_per_degree * 1000.0 )^2;  % q2 or source level per m^2
rho_SL_dB = 10 * log10( rho_SL );

rho_SL_dB = 54;   % monopole strength from Kewley for 800 Hz and 40 knots

Component = 'P';

NL = modal_noise_full( ModeFile, rho_SL_dB, sd, rd, Rmax_km, Component );

% line plots vs. disc size
% figure
% plot( Rmax_km, NL )
% xlabel( 'Radius of noise disc (km)' )
% ylabel( 'Noise Level (dB)' )

% image plot
figure
imagesc( Rmax_km, rd, NL )
colorbar
caxis( [ 60 80 ] )
xlabel( 'Range (km)' )
ylabel( 'Depth (m)' )
title( 'KRAKEN' )

%%
% Modal noise vs. depth using full matrix
% we just take the value from the matrix NL for the largest range

figure
plot( rd, NL( :, end ) )
xlabel( 'Receiver depth (m)' )
ylabel( 'NL (dB)' )
axis( [ 0 100 60 80 ] )
title( 'KRAKEN' )

%%
% Modal noise vs. depth using diagonal terms only

NL = modal_noise_diag( ModeFile, rho_SL_dB, sd, rd, Rmax_km, Component );
hold on
plot( rd, NL( :, end ), 'g' )


%%
% Spectral noise vs. depth
% Note that SCOOTER has to be run with stabilizing attenuation disabled
% ( TopOpt( 6 : 6 ) = '0' )

GreenFile = 'pekeris.grn';
scooter pekeris

NL = spectral_noise( GreenFile, rho_SL_dB, sd, rd );
hold on
plot( rd, NL, 'r' )

legend( 'NM Full', 'NM Diag', 'Spectral' )

%%

% tests for halfspace
% If the source is put at lambda/4 then
% the far field should have a level 13.2 dB higher than rho_SL_dB

GreenFile = 'halfspace.grn';

scooter halfspace

NL = spectral_noise( GreenFile, rho_SL_dB, sd, rd );

figure
plot( rd, NL )
xlabel( 'Receiver depth (m)' )
ylabel( 'NL (dB)' )
title( 'Halfspace, spectral formula' )

