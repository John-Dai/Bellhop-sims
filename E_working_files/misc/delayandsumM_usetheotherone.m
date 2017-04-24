%% E's version
% Convolve a source timeseries with the channel impulse response to
% get the received waveform.
% This is a delay-and-sum process in which echoes of the source waveform
% are combined based on their amplitude and arrival time
%
% mbp 8/96, Universidade do Algarve
% 4/09 addition of Doppler effects
% editing this version to work with bellhopM

%clear all
% deals with motion
v  = [ 2, 0 ];   % Receiver motion vector (vr, vz ) in m/s
%v  = [ 0, 0 ];   % Receiver motion vector (vr, vz ) in m/s
c0 = 1500;       % reference speed to convert v to proportional Doppler

%ARRFIL = 'autecArr.arr'; 

fileroot = 'simple1';
ARRFIL   = fileroot;
pcmfile  = [ fileroot '.pcm' ];

% source timeseries
STSFIL = 'srcsignal.wav'; %!!!
[ stsTemp, sample_rate ] = audioread( STSFIL );
sts = stsTemp( 1 : 2 : 5*sample_rate*2 );  % sub-sample down to sample_rate

% generated by 'cans'
% sample_rate = 20000;
% TT          = linspace( -0.2, 0.8, 20000 )';
% sample_rate = 20000;
% omega       = 2 * pi * 4000
% Pulse       = 'T'
% [ sts, PulseTitle ] = cans( TT, omega, Pulse );

% read from Sandipa's mat file
% load rate2_48kHz
% sts         = y_pb48k;
% sample_rate = 48000;
% deltat      = 1 / sample_rate;
% nts         = length( sts );
% TT          = linspace( 0, ( nts - 1 ) * deltat, nts );

% normalize the source time series so that it has 0 dB source level
% (measured based on the peak power in the timeseries)
sts = sts / max( abs( sts ) );
nts = length( sts );

%%

% optional Doppler effects
% pre-calculate Dopplerized versions of the source time series

% following can be further optimized if we know that ray-angle limits
% further restrict possible Doppler factors
optD = 0; %use this to control whether Doppler is calculated or not
% Doppler stuff takes too long - what's wrong?
if optD
     if ( norm( v ) > 0 )
         disp( 'Setting up Dopplerized waveforms' )
         v_vec = linspace( min( v ), max( v ), 10 );   % vector of Doppler velocity bins
         %v_vec = linspace( 1.9, 2, 51 );   % vector of Doppler velocity bins
         %what are these doing?
         alpha_vec   = 1 - v_vec / c0;                     % Doppler factors
         nalpha      = length( alpha_vec );
         sts_hilbert = zeros( nts, nalpha );
         % loop over Doppler factors (should be further vectorized ...)
         for ialpha = 1 : length( alpha_vec )
             %disp( ialpha )
             sts_hilbert( :, ialpha ) = hilbert( arbitrary_alpha( sts', alpha_vec( ialpha ) )' ); % Dopplerize
         end
     else
         sts_hilbert = hilbert( sts );   % need the Hilbert transform for phase changes (imag. part of amp)
     end
end
%%

%**************************************************************************

c = 1537.0;  % reduction velocity (should exceed fastest possible arrival)
             % read this: http://cires1.colorado.edu/people/jones.craig/GSSH/004-Refract/065-Reduction.htm
T = 5.0;     % time window to capture

Narrmx = 100; 

[ Arr, Pos ] = read_arrivals_bin( ARRFIL, Narrmx );  % read the arrivals file
disp( 'Done reading arrivals' ) 

% select which source/receiver index to use
ir  = length( Pos.r.range );
isd = length( Pos.s.depth );

deltat = 1 / sample_rate;
nt     = round( T / deltat );	% number of time points in rts
rtsmat = zeros( nt, length( Pos.r.range ) );

for ird = 1 : length( Pos.r.depth );
    for ir = 1 : length( Pos.r.range )   % loop over receiver ranges
        disp( [ ird, ir ] )
        tstart = Pos.r.range( ir ) / c - 0.1;   % min( delay( ir, :, ird ) )
        tend   = tstart + T - deltat;
        time   = tstart : deltat : tend;
        
        % compute channel transfer function
        
        rts = zeros( nt, 1 );	% initialize the time series
        
        for iarr = 1 : Arr.Narr( ir, ird, isd )   % loop over arrivals
            
            Tarr = Arr.delay( ir, iarr, ird, isd ) - tstart;   % arrival time relative to tstart
            
            it1  = round( Tarr / deltat + 1 );             % starting time index in rts for that delay
            it2  = it1 + nts - 1;                          % ending   time index in rts
            its1 = 1;                                      % starting time index in sts
            its2 = nts;                                    % ending   time index in sts
            
            % clip to make sure [ it1, it2 ] is inside the limits of rts
            if ( it1 < 1 )
                its1 = its1 + ( 1 - it1 );  % shift right by 1 - it1 samples
                it1  = 1;
            end
            
            if ( it2 > nt )
                its2 = its2 - ( it2 - nt );  % shift left by it2 - nt samples
                it2  = nt;
            end
            
            if ( norm( v ) > 0 )
                % identify the Doppler bin
                theta_ray = Arr.RcvrAngle( ir, iarr, ird ) * pi / 180;    % convert ray angle to radians
                tan_ray   = [ cos( theta_ray ) sin( theta_ray ) ];        % form unit tangent
                alpha     = 1 - dot( v / c0, tan_ray );                   % project Doppler vector onto ray tangent
                ialpha    = 1 + round( ( alpha - alpha_vec( 1 ) ) / ( alpha_vec( end ) - alpha_vec( 1 ) ) * ( length( alpha_vec ) - 1 ) );
                
                % check alpha index within bounds?
                if ( ialpha < 1 || ialpha > nalpha )
                    disp( 'Doppler exceeds pre-tabulated values' )
                    ialpha = max( ialpha, 1 );
                    ialpha = min( ialpha, nalpha );
                end
                
                % load the weighted and Dopplerized waveform into the received time series
                rts( it1 : it2 ) = rts( it1 : it2 ) + real( Arr.A( ir, iarr, ird ) * sts_hilbert( its1 : its2, ialpha ) );
            else
                
                %  rts( it1 : it2 ) = rts( it1 : it2 ) + real( Arr.A( ir, iarr, ird ) ) * real( sts_hilbert( its1 : its2, 1 ) ) ...
                %                                      - imag( Arr.A( ir, iarr, ird ) ) * imag( sts_hilbert( its1 : its2, 1 ) );
                
                % following is math-equivalent to above, but runs faster in
                % Matlab even though it looks like more work ...
                rts( it1 : it2 ) = rts( it1 : it2 ) + real( Arr.A( ir, iarr, ird ) * sts_hilbert( its1 : its2, 1 ) );
            end
        end   % next arrival, iarr
        
        %       % write to file
        %       if ( ird == 1 && ir == 1 )
        %          fid_out = fopen( pcmfile, 'w', 'ieee-le' );
        %       else
        %          fid_out = fopen( pcmfile, 'a', 'ieee-le' );
        %       end
        %       if ( fid_out == -1 )
        %          error('Can''t open PCM file for output.');
        %       end
        %
        %       fwrite( fid_out, 2^15 * 0.95 * rts / max( abs( rts ) ), 'int16' );
        %       fclose( fid_out );
        
        %wavwrite( 0.95 * rts / max( abs( rts ) ) , sample_rate, [ 'rts_Rd_' num2str( ird ) '_Rr_' num2str( ir ) '.wav' ] );
        %rtsmat( :, ir ) = 0.95 * rts / max( abs( rts ) );
        rtsmat( :, ir ) = rts;
    end
    eval( [ ' save ' fileroot '_Rd_' num2str( ird ) ' rtsmat Pos sample_rate' ] );
end

% save autec
% foo = reshape( rtsmat, 1, nt * length( Pos.r.range ) );
%
% figure; pcolor( Pos.r.range, linspace( 0, T, nt ), 20 * log10 ( abs( hilbert( rtsmat' ) ) ) ); shading flat; colorbar
% xlabel( 'Time (s)' )
% ylabel( 'Range (m)' )
