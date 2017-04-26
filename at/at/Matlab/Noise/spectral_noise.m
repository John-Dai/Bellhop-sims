function NL = spectral_noise( GreenFile, rho_SL_dB, sd, rd )
%
% Calculates the NL dues to a uniform sheet of monopoles
%
% This version works off a depth-separated Green's function calculated by
% SCOOTER
%
% rho_SL_dB source level density (dB expressed per m^2 of the noise sheet)
% GreenFile name of the Green's function file
% sd        source depth (m)
% rd        receiver depth (m) (can be a vector)
%
% returns
% NL in dB with nrd rows and nrr columns
%
% The source level density assumes source strength reference to the usual
% standard (a source that produces 1 microPascal as 1 m) scaled by the
% density per meter

% The derivation of this noise formula in JKPS assumes the integration is
% done along the real axis. Therefore stabilizing attenuation should not be
% used

% mike porter, 2012

Pos.s.depth = sd;  % depth of noise sources
Pos.r.depth = rd;   % vector of receiver depths
Nrd         = length( Pos.r.depth );

rho_SL = 10^( rho_SL_dB / 10 ); % q2 or source level per m^2

% *** read in the Green's function
[ PlotTitle, ~, freq, atten, PosG, Gtemp ] = read_shd( GreenFile );
k = PosG.r.range;
atten = 0;   % stabilizing attenuation not allowed !!!

% Get Green's function for given source and receiver depth

% FIXME: will fail if there is only one rd as interp1 needs two

% interpolate Green's fn at source depth
% note: interp1 won't go if there is only one point in depth

if length( PosG.s.depth ) > 1
   Gtemp2 = interp1( PosG.s.depth, squeeze( Gtemp( 1, :, :, : ) ), Pos.s.depth );
   Gtemp2 = squeeze( Gtemp2 );
else
   zs   = Pos.s.depth;
   isd  = find( PosG.s.depth >= zs );    % index of source depth
   isd  = isd( 1 );
   Gtemp2 = squeeze( Gtemp( 1, isd, :, : ) );
end

size( Gtemp2 )

G = interp1( PosG.r.depth, squeeze( Gtemp2( :, : ) ), Pos.r.depth );


%%

omega = 2 * pi * freq;
c     = 1500;        % nominal sound speed
k0    = omega / c;   % nominal vertical wavenumber

% *** Construct Power ***
% See Eq. (9.17), p. 667 in JKPS

dk = ( k( end ) - k( 1 ) ) / length( k );
NL = zeros( Nrd, 1 );

for ird = 1 : Nrd
   
   gamma = sqrt( k0^2 - ( k - 1i * atten ).^2 );
   %    if ( ird == 51 )
   %       figure; plot( k, abs( G( ird, : ) ), 'b' )
   %       hold on
   %    end
   %
   %    % analytic Green's function
   %    % You get the best answer if you set atten = 0
   %    G( ird, : ) = ( exp( 1i * gamma * abs( rd( ird ) - sd ) )  - ...
   %                    exp( 1i * gamma * abs( rd( ird ) + sd ) ) )./ gamma;
   %    if ( ird == 51 )
   %       plot( k, abs( G( ird, : ) ), 'r' )
   %    end
   
   G2 = abs( G( ird, : ) ).^2 * ( k - 1i * atten );
   NL( ird ) = G2 * dk;   % trapezoidal rule for integral over k
end

% Here's the factor to convert a SL density to q2
q2 = k0^2 * 4 * pi * rho_SL;
NL = ( 8 * pi^2 * q2 / k0^2 ) * NL;   % scale

NL = NL / ( 4 * pi )^2;    % additional scaling to take the SCOOTER definition of g to the jkps one.

NL = 10 * log10( NL );   % convert to dB
