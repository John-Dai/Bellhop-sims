function [rcv_ts, Tstart] = delay_sum_surface(ARRFIL, RAYFIL, ...
  src_huswf, src_usfs, surf_t, surf_r, surf_h, rcv_fs, v_src, v_rcv, chan_T)


% DELAY_SUM_SURFACE - delay and sum arrivals (with surface doppler)
%
% Usage: [rcv_ts, Tstart] = delay_sum_surface(ARRFIL, RAYFIL, ...
%                 src_huswf, src_usfs, surf_t, surf_r, surf_h, rcv_fs, chan_T)
%  Inputs:
%    ARRFIL    - string containing name of the BELLHOP arrivals file (binary)
%    RAYFIL    - string containing name of the BELLHOP eigenray file
%    src_huswf - amplitude samples of Hilbert transform of upsampled src packet
%    src_usfs  - sample rate of Hilbert transform of upsampled src packet (Hz)
%    surf_t    - time values corresponding to the surface height grid (seconds)
%    surf_r    - range values corresponding to the surface height grid (meters)
%    surf_h    - grid of surface wave heights (meters, 2D array, range by time)
%                or function handle f(time, range) (pass [] for surf_t, surf_r)
%    rcv_fs    - sample rate of the received time series (Hz)
%    v_src     - velocity (row) vector for source <v_r, v_z> (meters/second)
%    v_rcv     - velocity (row) vector for receiver <v_r, v_z> (meters/second)
%    chan_T    - channel time spread (seconds)
%
%  Output:
%    rcv_ts    - amplitude samples of the received time series
%    Tstart    - start time of the received time series (seconds)
%
% NOTES:
%
% 1) Arguments v_src, v_rcv, and chan_T are optional. If v_src and/or v_rcv
% are missing or empty, they are assumed to be zero. If chan_T is missing or
% empty, a value is computed using the BELLHOP arrivals.
%
% 2) For the src and rcv velocity vectors, v_r > 0 corresponds to increasing
% radial distance, and v_z > 0 corresponds to motion towards deeper depths.
%
% $Id: delay_sum_surface.m,v 1.7 2011/05/22 03:27:44 jcp Exp $

% check for missing or optional arguments

if nargin < 8,
  error( [mfilename, ': one or more required input arguments is missing'] );
end;

if nargin < 9,
  % if v_src is missing, it is assumed to be zero
  v_src = [0.0, 0.0];
else
  % if v_src is present, but empty, it is assumed to be zero
  if isempty(v_src),
    v_src = [0.0, 0.0];
  end;
end;

if nargin < 10,
  % if v_rcv is missing, it is assumed to be zero
  v_rcv = [0.0, 0.0];
else
  % if v_rcv is present, but empty, it is assumed to be zero
  if isempty(v_rcv),
    v_rcv = [0.0, 0.0];
  end;
end;

if nargin < 11,
  % if chan_T is missing we will compute a value later from BELLHOP...
  chan_T = 0.0;
else
  % also compute a value later if chan_T is present, but empty
  if isempty(chan_T),
    chan_T = 0.0;
  end;
end;

% derived values

src_vel = sqrt(dot(v_src, v_src));
rcv_vel = sqrt(dot(v_rcv, v_rcv));

h_isfloat = isfloat(surf_h);	% is surf_h a matrix? else assumed function

% assign values to "hard-wired" parameters

depthL = 1.0;

pad_T = 0.050;	% This is the time duration of extra samples added before
		% the start and after the end of the received (the channel
		% spread) waveform. This is to allow for the possibility of
		% both early and late arrivals due to small changes in the
		% arrival time due to surface motion. (Namely, the channel
		% spread can be slightly smaller and/or larger than what
		% BELLHOP predicts). Users DON'T need to worry about this,
		% as the code will check if this value was adequate or not...

% read the BELLHOP arrivals file

[Arr, Pos] = read_arrivals_bin(ARRFIL);

Nsd = Pos.Nsd;
Nrd = Pos.Nrd;
Nrr = Pos.Nrr;

rd = Pos.r.depth;
rr = Pos.r.range;

% verify there is only one source, one receiver (fix this restriction later...)

if Nsd ~= 1,
  err_str = sprintf('%s: number of sources (%d) in arrivals file: %s > 1\n', mfilename, Nsd, ARRFIL);
  error(err_str);
end;

if Nrd ~= 1,
  err_str = sprintf('%s: number of rcv depths (%d) in arrivals file: %s > 1\n', mfilename, Nrd, ARRFIL);
  error(err_str);
end;

if Nrr ~= 1,
  err_str = sprintf('%s: number of rcv ranges (%d) in arrivals file: %s > 1\n', mfilename, Nrr, ARRFIL);
  error(err_str);
end;

% get the number of arrivals
isd = 1;
ird = 1;	% for sake of clarity ...
irr = 1;
narr = Arr.Narr(irr, ird, isd);

% determine the min, max of all arrival times at this receiver
Arr_min = min( squeeze( Arr.delay(irr, 1:narr, ird, isd) ) );
Arr_max = max( squeeze( Arr.delay(irr, 1:narr, ird, isd) ) );

% compute a start time for the rcv time series based on BELLHOP arrival times
Tstart = max( [ Arr_min - pad_T, 0.0 ] );
Tstart = floor(Tstart * rcv_fs)/rcv_fs;

% verify the provided channel time spread is adequate (or compute one)
if chan_T > 0.0,
  if chan_T < (Arr_max - Arr_min + 2*pad_T),
    chan_str = sprintf('%.5f', Arr_max - Arr_min + 2*pad_T);
    error( [mfilename, ': channel spread is too small, must be at least: ', chan_str] );
  end;
else
  % compute a conservative value for the channel time spread (allow for doppler)
  chan_T = Arr_max - Arr_min + 2*pad_T;
  % round up to the nearest 1 millisecond (1/1000 = 0.001 seconds)
  chan_T = ceil(1000.0*chan_T)/1000.0;
end;

% compute the duration of the receiver time series
src_T = length(src_huswf) / src_usfs;
rcv_T = src_T + chan_T;	% add channel time spread to duration of the source

% allocate memory for the computation of arrival times at the receiver
nsamples_src = round(src_T * rcv_fs);
time_src = (0:(nsamples_src-1))/rcv_fs;

% allocate memory for the receiver (channel spread) time series
nsamples_rcv = round(rcv_T * rcv_fs);
rcv_ts = zeros(nsamples_rcv, 1);

% verify that the surface heights are defined over time interval of interest
if (h_isfloat),
  if ((Tstart + rcv_T) > surf_t(end)),
    error([ mfilename, ': surface is not defined over requested time span' ]);
  end;
end;

% read the BELLHOP eigenray file
[Hdr, Rays] = read_rayfil(RAYFIL);

depthT = Hdr.depthT;
depthB = Hdr.depthB;
Nrays = length(Rays);
for iray = 1:Nrays,
  Alpha0(iray) = Rays(iray).Alpha0;
end;

% loop over arrivals at this receiver

for iarr = 1:narr,

  % arrival time at receiver for this ray
  Tarr = Arr.delay(irr, iarr, ird, isd);

  % complex amplitude of this arrival
  CAmp = Arr.A(irr, iarr, ird, isd);

  % beam/ray angle at the source for this arrival
  SrcAngle = Arr.SrcAngle(irr, iarr, ird, isd);

  % arrival angle at the receiver for this arrival
  RcvAngle = Arr.RcvrAngle(irr, iarr, ird, isd);

  % number of surface interactions for this arrival
  NumTopBnc = Arr.NumTopBnc(irr, iarr, ird, isd);

  % number of bottom interactions for this arrival
  NumBotBnc = Arr.NumBotBnc(irr, iarr, ird, isd);

  % arrival/travel time of source samples at the receiver (at rcv sample rate)
  time_rcv = time_src;

  % adjust arrival time to account for steady source motion (if present)
  if (src_vel > 0.0),
    theta_ray = SrcAngle * pi / 180;
    tan_ray = [cos(theta_ray), sin(theta_ray)];	% tangent vector to ray
    dr_dt = -dot(v_src, tan_ray);	% time rate of change of path length
    time_rcv = time_rcv + dr_dt * (time_rcv - min(time_rcv)) / 1500.0;
  end;

  % check if this ray interacts with the surface

  if NumTopBnc > 0,

    % locate the ray associated with this arrival
    [min_v, min_i] = min( abs(Alpha0 - SrcAngle) );
    Nsteps = Rays(min_i).Nsteps;
    r = Rays(min_i).r;
    d = Rays(min_i).d;

    % compute the (total) path length along the ray (which goes "past" rcv)
    del_r = r(2:Nsteps) - r(1:(Nsteps-1));
    del_d = d(2:Nsteps) - d(1:(Nsteps-1));
    del_path_len = sqrt(del_r.*del_r + del_d.*del_d);
    path_len = [ 0.0, cumsum(del_path_len) ];

    % compute the ray path lengths between surface interactions (or src/rcv)
    num_segs = 0;
    iNumTopBnc = 0;
    last_r = r(1);
    last_i = 1;

    % find the next interaction of the ray with the surface
    while iNumTopBnc < NumTopBnc,
      ibeg = last_i;
      % ignore first segment starting at the source
      if iNumTopBnc > 0,
        while d(ibeg) < (depthT+depthL),
          % still near the surface from the previous interaction
          ibeg = ibeg + 1;
        end;
      end;
      % march along the ray until we get near the surface again
      while d(ibeg) >= (depthT+depthL),
        % still too far below the surface
        ibeg = ibeg + 1;
      end;
      % march along the ray until we get below the surface again
      iend = ibeg;
      while d(iend+1) < (depthT+depthL),
        iend = iend + 1;
      end;
      % choose the index with minimum depth
      [min_v, min_i] = min(d(ibeg:iend));
      istep = ibeg + min_i - 1;
      % save all the attributes of the interaction that we will need later
      num_segs = num_segs + 1;
      seg_beg_r(num_segs) = last_r;
      seg_beg_cos(num_segs) = cos_ray(r, d, last_i, +1);
      seg_end_r(num_segs) = r(istep);
      seg_end_cos(num_segs) = cos_ray(r, d, istep, -1);
      seg_len(num_segs) = path_len(istep) - path_len(last_i);
      % next surface interaction
      iNumTopBnc = iNumTopBnc + 1;
      last_r = r(istep);
      last_i = istep;
    end;

    % compute the path length from the last surface interaction to the receiver
    del_r = r(last_i:Nsteps) - Pos.r.range(irr);
    del_d = d(last_i:Nsteps) - Pos.r.depth(ird);
    dist_to_rcv_sqr = del_r.*del_r + del_d.*del_d;
    [min_v, min_i] = min(dist_to_rcv_sqr);
    istep = last_i + min_i - 1;

    num_segs = num_segs + 1;
    seg_beg_r(num_segs) = last_r;
    seg_beg_cos(num_segs) = cos_ray(r, d, last_i, +1);
    seg_end_r(num_segs) = Pos.r.range(irr);
    seg_end_cos(num_segs) = cos_ray(r, d, istep, -1);
    seg_len(num_segs) = path_len(istep) - path_len(last_i);

    % compute the total path length along the ray from src -> rcv
    total_len = sum(seg_len(1:num_segs));

    % compute the mean sound speed over the path from src -> rcv
    cbar = total_len / Tarr;

    % estimated arrival time (will adjust later to force it to be Tarr)
    est_Tarr = 0.0;

    % work forward in time from the src and compute the arrival times

    for iseg = 1:num_segs,
      % beg of segment
      if (iseg > 1),
        beg_r = seg_beg_r(iseg);
        beg_cos = seg_beg_cos(iseg);
        if h_isfloat,
          % surf_h is a matrix, estimate h using 2-D interpolation in t, r
          del_seg = interp2(surf_t, surf_r, surf_h, time_rcv, beg_r) * beg_cos;
        else
          % surf_h is a function handle, call it with desired t, r values
          del_seg = surf_h(time_rcv, beg_r) * beg_cos;
        end;
        % add the extra (estimated) travel time for this surface interaction
        time_rcv = time_rcv + del_seg/cbar;
      end;
      % add the (estimated) travel time for this segment of the total path
      time_rcv = time_rcv + seg_len(iseg)/cbar;
      % keep track of the estimated travel time without surface motion
      est_Tarr = est_Tarr + seg_len(iseg)/cbar;
      % end of segment
      if (iseg < num_segs),
        end_r = seg_end_r(iseg);
        end_cos = seg_end_cos(iseg);
        if h_isfloat,
          % surf_h is a matrix, estimate h using 2-D interpolation in t, r
          del_seg = interp2(surf_t, surf_r, surf_h, time_rcv, end_r) * end_cos;
        else
          % surf_h is a function handle, call it with desired t, r values
          del_seg = surf_h(time_rcv, end_r) * end_cos;
        end;
        % add the extra (estimated) travel time for this surface interaction
        time_rcv = time_rcv + del_seg/cbar;
      end;
    end;

    % tweak the arrival time to match BELLHOP when there is no surface motion
    time_rcv = time_rcv + (Tarr - est_Tarr);

  else

    % handle the case where the ray does not interact with the surface
    time_rcv = time_rcv + Tarr;

  end;

  % adjust arrival time to account for steady receiver motion (if present)
  if (rcv_vel > 0.0),
    theta_ray = RcvAngle * pi / 180;
    tan_ray = [cos(theta_ray), sin(theta_ray)];	% tangent vector to ray
    dr_dt = dot(v_rcv, tan_ray);	% time rate of change of path length
    time_rcv = time_rcv + dr_dt * (time_rcv - min(time_rcv)) / 1500.0;
  end;

  % determine the indices where we need to compute the rcv time series
  itime_beg =  ceil((min(time_rcv)-Tstart) * rcv_fs) + 1;
  itime_end = floor((max(time_rcv)-Tstart) * rcv_fs) + 1;

  % verify the arrival times aren't outside of the rcv array (is pad_T valid?)
  if (itime_beg < 1 || itime_end > nsamples_rcv),
    error( [mfilename, ': internal error, increase pad_T and try again'] );
  end;

  % time values where the rcv time series will be sampled (at rcv_fs rate)
  times_rcv = Tstart + ((itime_beg-1):(itime_end-1))/rcv_fs;

  % use interpolation to determine the time when these samples left the source
  times_src = interp1(time_rcv, time_src, times_rcv, 'linear');

  % compute the indices of the source waveform at the delayed time values
  isamples = round(times_src * src_usfs) + 1;	% Matlab arrays start at 1

  % accumulate the contribution from this arrival

  rcv_ts(itime_beg:itime_end) = rcv_ts(itime_beg:itime_end) ...
    + ( real(CAmp)*real( src_huswf(isamples) ) ) ...
    - ( imag(CAmp)*imag( src_huswf(isamples) ) );

end;	% next arrival, iarr
      
% that's all folks!!!

return;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% function to estimate the cosine of the incoming/outgoing ray angle
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function cos_th = cos_ray(range, depth, ibeg, flag)

fudge = 2.0;

if (flag > 0),

  % outgoing ray

  iend = ibeg;
  path_len = 0.0;
  while (path_len < fudge),
    iend = iend + 1;
    d_range = range(iend) - range(ibeg);
    d_depth = depth(iend) - depth(ibeg);
    path_len = sqrt(d_range^2 + d_depth^2);
    cos_th = abs(d_depth) / path_len;
  end;

else

  % incoming ray

  iend = ibeg;
  path_len = 0.0;
  while (path_len < fudge),
    ibeg = ibeg - 1;
    d_range = range(iend) - range(ibeg);
    d_depth = depth(iend) - depth(ibeg);
    path_len = sqrt(d_range^2 + d_depth^2);
    cos_th = abs(d_depth) / path_len;
  end;

end;

return;

%
% end of delay_sum_surface.m
