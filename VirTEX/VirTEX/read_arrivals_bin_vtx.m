function [Arr, Pos] = read_arrivals_bin(ARRFIL)

% Read the arrival time/amplitude data computed by BELLHOP
%
% Usage:
%
%   [ Arr, Pos ] = read_arrivals_bin( ARRFIL );
%
% Arr is a structure containing all the arrivals information
% Pos is a structure containing the positions of source and receivers
%
% ARRFIL is the name of the Arrivals File
% mbp 9/96

% open the arrivals file (read-only)

fid = fopen(ARRFIL, 'r');	% open the file

if fid < 0,
  error( [mfilename, ': error opening ARRFIL: ', ARRFIL] );
end;

% read the header information

fseek( fid, 4, -1 );
freq = fread( fid, 1, 'float' );
Nsd  = fread( fid, 1, 'long'  );
Nrd  = fread( fid, 1, 'long'  );
Nrr  = fread( fid, 1, 'long'  );

Pos.freq = freq;
Pos.Nsd = Nsd;
Pos.Nrd = Nrd;
Pos.Nrr = Nrr;

fseek( fid, 8, 0 );
Pos.s.depth   = fread( fid, Nsd, 'float' );

fseek(fid,8,0);
Pos.r.depth   = fread( fid, Nrd, 'float' );

fseek(fid,8,0);
Pos.r.range   = fread( fid, Nrr, 'float' );

fseek(fid,8,0);
Narrmax = fread(fid, 1, 'int', 8);

% loop to read all the arrival info (delay and amplitude)

Arr.Narr      = zeros( Nrr, Nrd, Nsd );
Arr.A         = zeros( Nrr, Narrmax, Nrd, Nsd );
Arr.delay     = zeros( Nrr, Narrmax, Nrd, Nsd );
Arr.SrcAngle  = zeros( Nrr, Narrmax, Nrd, Nsd );
Arr.RcvrAngle = zeros( Nrr, Narrmax, Nrd, Nsd );
Arr.NumTopBnc = zeros( Nrr, Narrmax, Nrd, Nsd );
Arr.NumBotBnc = zeros( Nrr, Narrmax, Nrd, Nsd );

for isd = 1:Nsd,
   for ird = 1:Nrd,
      for irr = 1:Nrr,
         Narr = fread( fid, 1, 'int', 8 );
         Arr.Narr( irr, ird, isd ) = Narr;
         if Narr > 0   % do we have any arrivals?
            da = fread( fid, [9, Narr], 'float');
            Arr.A(         irr, 1:Narr, ird, isd ) = da( 1, 1:Narr ) .* exp( 1i * da( 2, 1:Narr ) * pi/180 );
            Arr.delay(     irr, 1:Narr, ird, isd ) = da( 3, 1:Narr );
            Arr.SrcAngle(  irr, 1:Narr, ird, isd ) = da( 4, 1:Narr );
            Arr.RcvrAngle( irr, 1:Narr, ird, isd ) = da( 5, 1:Narr );
            Arr.NumTopBnc( irr, 1:Narr, ird, isd ) = da( 6, 1:Narr );
            Arr.NumBotBnc( irr, 1:Narr, ird, isd ) = da( 7, 1:Narr );
         end
      end		% next receiver range
   end		% next receiver depth
end	% next source depth

fclose( fid );
