function rtsmat = delay_sum(ARRFIL, src_wf, v_src, v_rcv, Trcv)

% DELAY_SUM - convolve a source time series with the channel impulse response
%
% This is VirTEX for Platform Motion which estimates a receiver time series 
% by convolving the source time series with the channel impulse response
% computed by BELLHOP. The case of a moving receiver is supported by
% using a library of "dopplerized" versions of the source waveform.
%
% Usage: rtsmat = delay_sum(ARRFIL, src_wf, v_src, v_rcv, Trcv)
%
%   ARRFIL - filename of the (ascii) BELLHOP arrivals file
%   src_wf - struct containing library of dopplerized src waveform
%    v_src - velocity (row) vector for source <v_r, v_z> (meters/second)
%    v_rcv - velocity (row) vector for receiver <v_r, v_z> (meters/second)
%     Trcv - total duration of receiver time series, optional (seconds)
%
% NOTES:
%
% 1) For the src and rcv velocity vectors, v_r > 0 corresponds to increasing
% radial distance, and v_z > 0 corresponds to motion towards deeper depths.
%
% $Id: delay_sum.m,v 1.5 2011/07/08 18:12:55 jcp Exp $

isd = 1;	% source depth to use, hard-wired for now

% check for missing or optional arguments

if nargin < 4,
  error( [mfilename, ': one or more required input arguments is missing'] );
end;

if nargin < 5,
  % we will compute a value later from all the BELLHOP arrivals...
  Trcv = 0.0;
end;

% read the (binary) arrivals file

[Arr, Pos] = read_arrivals_bin(ARRFIL);

rr = Pos.r.range;
rd = Pos.r.depth;
nrr = length(rr);
nrd = length(rd);

% unpack the source waveform struct produced by delay_sum_lib

c0 = src_wf.c0;
fs = src_wf.fs;
alpha_vec = src_wf.alpha_vec;
sts_hilbert = src_wf.sts_hilbert;
nsamples_src = size(sts_hilbert, 1);
nalpha = length(alpha_vec);
Tsrc = nsamples_src / fs;

% compute the maximum channel time spread (over all receivers)

ch_max = zeros(nrr, nrd);

for ird = 1:nrd
  for irr = 1:nrr
     %E. Oon - fixing indexing problems 
    %narr = Arr.Narr(irr, ird, isd);
    narr = Arr.Narr(ird,irr,isd);
    if (narr > 0)
      %ch_max(irr,ird) = max(squeeze( Arr.delay(irr, 1:narr, ird, isd) )) ...
                     % - min(squeeze( Arr.delay(irr, 1:narr, ird, isd) ));
        iDelay = Arr.delay(ird, irr,1:narr, isd);
        ch_max(irr,ird) = max(squeeze(iDelay)) - min(squeeze(iDelay));
    end
  end
end

ch_max = max(max(ch_max));

% verify the duration of the receiver time series is adequate (or compute one)

Tmin = Tsrc + ch_max;

if Trcv > 0.0,
  if Trcv < Tmin,
    Tmin_str = sprintf('%.5f', Tmin);
    error( [mfilename, ': duration of the receiver time series is too small, must be at least: ', Tmin_str] );
  end;
else
  % compute the smallest possible value
  Trcv = Tmin;
  % round up to the nearest 1 millisecond (1/1000 = 0.001 seconds)
  Trcv = ceil(1000.0*Trcv)/1000.0;
end;

% allocate memory for the receiver (channel spread) time series

nsamples_rcv = ceil(fs * Trcv);

rtsmat = zeros(nsamples_rcv, nrr, nrd);

% loop over receiver depths

for ird = 1:nrd,

  % loop over receiver ranges

  for irr = 1:nrr,

    % number of arrivals for this receiver
    
    %narr = Arr.Narr(irr, ird, isd); gotta fix indexing 
    narr = Arr.Narr(ird,irr,isd);

    % find the min of all arrival times at this receiver
    iDelay = Arr.delay(ird, irr,1:narr, isd); %fixing indexing 
    Arr_min = min(squeeze(iDelay));

    % initialize the time series for this receiver

    rts = zeros(nsamples_rcv, 1);

    % loop over the arrivals for this receiver

    for iarr = 1:narr

      % arrival time relative to start of receiver time series

      Tarr = Arr.delay(irr, iarr, ird, isd) - Arr_min;

      % compute the starting and ending time sample indices

      it1  = round(fs * Tarr) + 1;

      it2  = it1 + nsamples_src - 1;

      % identify which dopplerized source waveform to use from the library

      theta_ray = Arr.SrcAngle(irr, iarr, ird) * pi / 180;

      tan_ray = [cos(theta_ray), sin(theta_ray)];	% tangent vector to ray

      plroc_src = dot(v_src, tan_ray);	% path length rate of change (for src)

      theta_ray = Arr.RcvrAngle(irr, iarr, ird) * pi / 180;

      tan_ray = [cos(theta_ray), sin(theta_ray)];	% tangent vector to ray

      plroc_rcv = dot(v_rcv, tan_ray);	% path length rate of change (for rcv)

      alpha = 1.0 - (plroc_rcv - plroc_src)/c0;

      phi = (alpha - alpha_vec(1)) / (alpha_vec(end) - alpha_vec(1));

      ialpha = 1 + round(phi * (nalpha-1));

      % verify that the alpha index is within bounds

      if (ialpha < 1 || ialpha > nalpha),
        error('Doppler exceeds the pre-tabulated values')
      end

      % add the contribution from this arrival

      rts(it1:it2) = rts(it1:it2) ...
           + real( Arr.A(irr, iarr, ird) * sts_hilbert(:, ialpha) );

    end   % next arrival, iarr

    % save the receiver time series

    rtsmat(:, irr, ird) = rts;

  end

end

return;

%
% end of delay_sum.m
