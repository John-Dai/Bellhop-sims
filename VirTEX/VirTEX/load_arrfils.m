function arrival_params = load_arrfils(ARRFIL_fmt, bh_time);

% LOAD_ARRFILS - load a sequence of (binary) BELLHOP ray arrival files.
%
% USAGE:
%
%   arrival_params = load_arrfils(ARRFIL_fmt, bh_time);
%
% INPUTS
%
%   ARRFIL_fmt ==> An sprintf compatible format that specifies the names
%                  of the BELLHOP arrival file names. It is assumed that
%                  files with names sprintf(ARRFIL_fmt, n) exist and are
%                  readable for values of n = 1:length(bh_time)
%
%   bh_time    ==> Array of time values associated with the arrival files.
%
% OUTPUTS:
%
%   arrival_params ==> Structure containing the ray arrival data. It is
%                      compatible with the moving_environment function
%                      that implements the VIRTEX algorithm.
%
% USAGE EXAMPLE:
%
%   !bellhop < tstest_01.env > tstest_01.prt ; mv ARRFIL tstest_01.arr
%   !bellhop < tstest_02.env > tstest_02.prt ; mv ARRFIL tstest_02.arr
%   !bellhop < tstest_03.env > tstest_03.prt ; mv ARRFIL tstest_03.arr
%   bh_time = [0.0 0.1 0.2];
%   ARRFIL_fmt = 'tstest_%02d.arr';
%   arrival_params = load_arrfils(ARRFIL_fmt, bh_time);
%
% M. Siderius Aug. 2002
%
% $Id: load_arrfils.m,v 1.4 2011/05/12 01:31:31 jcp Exp $

i = sqrt(-1.0d0);		% everyone knows "j" is an index !!!

pi = acos(-1.0d0);

deg2rad = pi / 180.0d0;

% loop over all time values in bh_time

for itime = 1:length(bh_time),

  % open the file read-only

  filename = sprintf(ARRFIL_fmt, itime);

  fid = fopen(filename, 'r', 'native');

  if fid < 0,
    error( ['load_arrfils: unable to open: ', filename] );
  end;

  % read the header info

  fseek(fid, 4, -1);
  freq = fread(fid, 1, 'float');
  nsd  = fread(fid, 1, 'long');
  nrd  = fread(fid, 1, 'long');
  nrr  = fread(fid, 1, 'long');

  fseek(fid, 8, 0);
  sd   = fread(fid, nsd, 'float');

  fseek(fid, 8, 0);
  rd   = fread(fid, nrd, 'float');

  fseek(fid, 8, 0);
  rr   = fread(fid, nrr, 'float');

  fseek(fid, 8, 0);
  narrmax = fread(fid, 1, 'int', 8);

  % check for more than one source (not supported at this time)

  if nsd > 1,
    error('load_arrfils: more than one source detected');
  end;

  % sanity checks for consistency from frame to frame

  if itime > 1,

    % verify that the receiver grid size has not changed

    if (nsd ~= arrival_params.nsd),
      error('load_arrfils: BELLHOP runs have different number of sources');
    end;

    if (nrd ~= arrival_params.nrd),
      error('load_arrfils: BELLHOP runs have different number of rcv depths');
    end;

    if (nrr ~= arrival_params.nrr),
      error('load_arrfils: BELLHOP runs have different number of rcv ranges');
    end;

  else

    % store the source and receiver grid parameters

    arrival_params.sd = sd;
    arrival_params.rr = rr;
    arrival_params.rd = rd;
    arrival_params.nsd = nsd;
    arrival_params.nrd = nrd;
    arrival_params.nrr = nrr;

  end;

  % allocate memory for amp, delay, etc.

  narrmat   = zeros(nrr, nrd);

  amp       = zeros(narrmax, nrr, nrd);
  delay     = zeros(narrmax, nrr, nrd);
  SrcAngle  = zeros(narrmax, nrr, nrd);
  RcvAngle  = zeros(narrmax, nrr, nrd);
  NumTopBnc = zeros(narrmax, nrr, nrd);
  NumBotBnc = zeros(narrmax, nrr, nrd);

  % loop over all the receiver ranges, depths

  for ird = 1:nrd,

    for irr = 1:nrr,
  
      narr = fread(fid, 1, 'int', 8);

      narrmat(irr, ird) = narr;
  
      % do we have any arrivals?

      if narr > 0,
        data = fread( fid, [9, narr], 'float');
        phase = exp(i * data(2, :) * deg2rad);
        amp(      1:narr, irr, ird) = data(1, :) .* phase;	% complex
        delay(    1:narr, irr, ird) = data(3, :);
        SrcAngle( 1:narr, irr, ird) = data(4, :);
        RcvAngle( 1:narr, irr, ird) = data(5, :);
        NumTopBnc(1:narr, irr, ird) = data(6, :);
        NumBotBnc(1:narr, irr, ird) = data(7, :);
      end;

    end;	% next receiver range

  end;	% next receiver depth

  fclose(fid);

  % store the reorganized arrival data

  irr = 1:nrr;
  ird = 1:nrd;

  maxarr = max(max(narrmat));

  if maxarr ~= narrmax,
    error( [mfilename, ': internal error, maxarr is not equal to narrmax!'] );
  end;

  arrival_params.narrmat(irr,ird,itime) = narrmat;

  arrival_params.amp(      1:maxarr,irr,ird,itime) = amp(      1:maxarr,:,:);
  arrival_params.delay(    1:maxarr,irr,ird,itime) = delay(    1:maxarr,:,:);
  arrival_params.SrcAngle( 1:maxarr,irr,ird,itime) = SrcAngle( 1:maxarr,:,:);
  arrival_params.RcvAngle( 1:maxarr,irr,ird,itime) = RcvAngle( 1:maxarr,:,:);
  arrival_params.NumTopBnc(1:maxarr,irr,ird,itime) = NumTopBnc(1:maxarr,:,:);
  arrival_params.NumBotBnc(1:maxarr,irr,ird,itime) = NumBotBnc(1:maxarr,:,:);

  arrival_params.bh_time(itime) = bh_time(itime);

end;	% next time

return;

% end of load_arrfils.m
