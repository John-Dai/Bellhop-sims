
% Fourier synthesis to make a time series from the transfer function
% mbp 9/96
% This version updated 2014 to be compatible with the current file formats

% Need to set
%   Tstart = starting time
%   specfile (file containing spectrum)
%   root of field files
%   N = size of transform

clear all

Tstart   = 19.0;
specfile = 'munkts.fft';
shdfile  = 'sductB';
N        = 1000;

%%
% load source spectrum

%fid   = fopen( specfile );
%temp  = fscanf( fid, '%f', [ 3, inf ] );
%shat  = temp( 2, : ) + 1i * temp( 3, : );

'temporary force of a unit spectrum'
shat = ones( 91, 1 );

nfreq = length( shat );

nsd = 1;

for isd = 1 : nsd
   
   % read in the model transfer function
   
   clear rmodhat
   
   freq = zeros( nfreq, 1 );
   
   for ifreq = 1 : nfreq
      
      % read the shdfile
      
      filename = [ shdfile int2str( ifreq ) '.shd.mat' ];
      [ PlotTitle, PlotType, freqt, atten, Pos, pressure ] = read_shd( filename );
      
      nrd = length( Pos.r.depth );
      nrr = length( Pos.r.range );
      
      if ( ifreq == 1 )   % allocate space
         rmodhat = zeros( nfreq, nrd, nrr );
      end
      
      freq( ifreq ) = freqt;
      rmodhat( ifreq, 1 : nrd, 1 : nrr ) = squeeze( pressure );
      
      imagesc( abs( squeeze( rmodhat( ifreq, :, : ) ) ) )
      drawnow
   end   % next frequency
   
   %%
   % compute the received time series
   
   % remove the travel time delay
   
   for ir = 1 : nrr
      for ird = 1 :nrd
         rmodhat( :, ird, ir ) = rmodhat( :, ird, ir ) .* exp( -1i * 2 * pi * Tstart * freq );
      end
   end
   
   % weight transfer function by source spectrum
   
   for ifreq = 1 : nfreq
      rmodhat( ifreq, :, : ) = rmodhat( ifreq, :, : ) * shat( ifreq );
   end
   
   rmod = ifft( rmodhat, N );   % inverse FFT to compute received time-series
   
   % set up time vector based on usual FFT sampling rules
   
   deltaf = freq( 2 ) - freq( 1 );
   Tmax   = 1 / deltaf;
   deltat = Tmax / N;
   time   = linspace( 0.0, Tmax - deltat, N );
   
   %%
   % plot
   figure
   
   for itime = 1 : N
      imagesc( Pos.r.range / 1000, Pos.r.depth, 10 * log10( abs( squeeze( rmod( itime, :, : ) ) ) ) )
      xlabel( 'Range (km)' )
      ylabel( 'Depth (m)' )
      drawnow
   end
   
   %%
   % heterodyne with the base frequency
   % Be careful that rmod has adequate time-sampling ...
   
   for ir = 1 : nrr
      for ird = 1 : nrd
         rmod( :, ird, ir ) = rmod( :, ird, ir ) .* exp(  1i * 2 * pi * time' * freq( 1 ) );
      end
   end
   
   % spectrum is conjugate symmetric so:
   rmod = 2 * real( rmod );
   
   figure
   imagesc( time, Pos.r.depth, abs( rmod.' ) )
   xlabel( 'Time (s)' )
   ylabel( 'Depth (m)' )
   
   %%
   % plot a specific receiver
   
   %    irr = 1;
   %    rmod = squeeze( rmod( :, irr, : ) );	% can only plot 2D matrix
   %
   %    figure; plot( rmod( :, 26 ) )
   %
   %    figure
   %    peak = max( max( rmod ) );
   %    % imagesc( time, rd, squeeze( rmod )' );
   %    nint = 1000;
   %    imagesc( time, rd, rmod' ); ...
   %       caxis( [ -peak/5, peak/5 ] ); colorbar
   %    xlabel( 'Time (s)' )
   %    ylabel( 'Depth (m)' )
   %    title( 'KRAKEN impulse response (dB)' )
   %
   %rmod = abs( hilbert( rmod( :, : ) ) );
   
   % normalize
   
   %for ir = 1:nrr
   %  temp = 20 * log10( rmod( :, ir ) / max( rmod( :, ir ) ) ) + 30;
   %  I = find( temp < 0 );
   %  temp( I ) = zeros( size( I ) );
   %  rmod( :, ir ) = temp / norm( temp );
   %end
   %rmod = rmod';
   
   % save for future use ...
   
   %fname = [ 'model' int2str( isd ) ]
   %eval( [ 'save ' fname ' nrr rr rmod' ] );
   
end   % next source depth
