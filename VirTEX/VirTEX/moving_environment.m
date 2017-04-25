function [rcv_ts, Tstart] = moving_environment(rcv_pos, ...
                  src_huswf, src_usfs, arrival_params, rcv_fs, chan_T)

% MOVING_ENVIRONMENT - implements VIRTEX for a moving environment
%
% A simulator that computes the received time series after propagating
% through a sound channel for a given source waveform. Any of the source,
% receiver, or sea surface can be either static or moving.
%
% USAGE:
%
%  [rcv_ts, Tstart] = moving_environment(rcv_pos, ...
%           src_huswf, src_usfs, arrival_params, rcs_fs, chan_T)
%
% INPUTS:
%
%  rcv_pos        - Nx3 array that describes the position of the receiver
%                   versus time. The first column specifies the receiver
%                   range (meters), the second column specifies the
%                   receiver depth (meters), and the third column
%                   specifies the time (seconds).
%  src_huswf      - hilbert transform of upsampled src packet
%  src_usfs       - sample rate of upsampled src packet (Hz)
%  arrival_params - arrival amplitudes and delays (from load_arrfils)
%  rcv_fs         - sample rate of the received time series (Hz)
%  chan_T         - channel time spread, optional (seconds)
%
% OUTPUTS:
%  rcv_ts         - amplitude samples of the received time series
%  Tstart         - start time of the received time series (seconds)
%
% NOTES:
%  Make sure to compute a time-range-depth grid in Bellhop that contains
%  all of the possible time-range-depth positions of the receiver. That
%  is, make sure the receiver doesn't move outside of the pre-computed
%  time-range-depth grid.
%
%  M. Siderius August, 2002
%
% $Id: moving_environment.m,v 1.7 2011/04/13 18:42:44 jcp Exp $

% check for the required arguments

if (nargin < 5),
  error('moving_environment: one or more required arguments is missing?');
end

% check for the optional chan_T argument

if (nargin < 6),
  % we will compute a value later from BELLHOP...
  chan_T = 0.0;
end

% extract required parameters from arrival_params

rr = arrival_params.rr;
rd = arrival_params.rd;

nrr = length(rr);
nrd = length(rd);

if isfield(arrival_params, 'bh_time'),
  bh_time = arrival_params.bh_time;
else
  bh_time = 0.0;
end

if isfield(arrival_params, 'c0'),
  c0 = arrival_params.c0;
else
  c0 = 1500.0;
end

amp = arrival_params.amp;
delay = arrival_params.delay;
narrmat = arrival_params.narrmat;
src_angle = (pi/180.0) * arrival_params.SrcAngle;
rcv_angle = (pi/180.0) * arrival_params.RcvAngle;
numtopbnc = arrival_params.NumTopBnc;
numbotbnc = arrival_params.NumBotBnc;

clear arrival_params

% determine the min, max of all arrival times at the receiver(s)

Arr_min = min(min(min(min(delay(delay > 0.0)))));
Arr_max = max(max(max(max(delay(delay > 0.0)))));

Arr_min = floor(1000.0*Arr_min)/1000.0;	% round down to nearest 1 millisecond
Arr_max =  ceil(1000.0*Arr_max)/1000.0;	% round   up to nearest 1 millisecond

% verify the user provided channel time spread is adequate (or compute one)

if chan_T > 0.0,
  if chan_T < (Arr_max - Arr_min),
    chan_str = sprintf('%.6f', Arr_max - Arr_min);
    error( ['moving_environment: channel spread is too small, must be at least: ', chan_str] );
  end;
else
  % compute a value for the channel time spread
  chan_T = Arr_max - Arr_min;
end;

% compute a start time for the rcv time series based on BELLHOP arrival times

Tstart = Arr_min;

Tstart = floor(Tstart * rcv_fs)/rcv_fs;

% compute the duration of the receiver time series

nsamples_src = length(src_huswf);

src_T = nsamples_src / src_usfs;

rcv_T = src_T + chan_T;

% verify the received waveform is within the time spanned by the BELLHOP runs

if ((Tstart < min(bh_time)) || ((Tstart+rcv_T) >= max(bh_time))),
  error('moving_environment: received time is outside of BELLHOP time grid')
end

% compute grid spacing in (BELLHOP) time dimension (uniform spacing is assumed)

nframes = length(bh_time);

if (length(bh_time) > 1);
  dbh_time = (bh_time(end) - bh_time(1))./(length(bh_time) - 1);
else
  amp(:,:,:,2) = amp(:,:,:,1);
  delay(:,:,:,2) = delay(:,:,:,1);
  dbh_time = 1;
end

% verify that the reciever position stays within of the BELLHOP reciever grid

if (min(rcv_pos(:,1)) < min(rr) || max(rcv_pos(:,1)) > max(rr)),
  error('moving_environment: track range moves outside of BELLHOP range grid')
end

if (min(rcv_pos(:,2)) < min(rd) || max(rcv_pos(:,2)) > max(rd)),
  error('moving_environment: track depth moves outside of BELLHOP depth grid')
end

% compute the time values where the received time series is to be computed

nsamples = round(rcv_T * rcv_fs);

t_rcv = Tstart + (0:(nsamples-1)) / rcv_fs;

rcv_ts = zeros(nsamples, 1);

% compute the receiver position at each time value (using linear interpolation)

r_rcv = interp1(rcv_pos(:,3), rcv_pos(:,1), t_rcv, 'linear', 'extrap');
d_rcv = interp1(rcv_pos(:,3), rcv_pos(:,2), t_rcv, 'linear', 'extrap');

% For the sake of efficiency, the case of a single receiver, and the case of
% a grid of receivers (at least two ranges, two depths) are treated separately.
% While there is much duplication of code, this eliminates "if" statements
% inside "for" loops, which tends to slow things down quite a bit.

hilb_uswf_r = real(src_huswf);
hilb_uswf_i = imag(src_huswf);

if (nrr == 1 && nrd == 1),

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Case of a single receiver
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% verify there is no receiver motion

if ( (min(rcv_pos(:,1)) < max(rcv_pos(:,1))) ...
  || (min(rcv_pos(:,2)) < max(rcv_pos(:,2))) ),
  error('moving_environment: no motion allowed for a single BELLHOP receiver')
end

% allocate memory for the interpolated arrival amp, delay values

maxarr = max(max(max( narrmat )));

ampi   = zeros(maxarr, 2);
delayi = zeros(maxarr, 2);

% loop over the time values where the received time series is to be computed

for isample = 1:nsamples,

  % find the index of the bh_time interval that contains this time sample

  it_1 = 1 + floor((t_rcv(isample) - bh_time(1)) / dbh_time);

  it_2 = it_1 + 1;

  % determine the maximum number of arrivals (over both ends of this interval)

  narrivals = [ narrmat(1, 1, it_1), ...
                narrmat(1, 1, it_2) ];

  maxarr = max(narrivals);

  iarr = 1:maxarr;

  % compute the interpolation weight with respect to bh_time

  twgt = (t_rcv(isample) - bh_time(it_1)) / dbh_time;

  % compute the (linear interpolation) ray amplitude weights

  ampi(iarr, 1) = (1.0-twgt)*amp(iarr, 1, 1, it_1);

  ampi(iarr, 2) = (    twgt)*amp(iarr, 1, 1, it_2);

  % get the BELLHOP time delays that correspond to each of the end points

  delayi(iarr, 1) = delay(iarr, 1, 1, it_1);	% time 1

  delayi(iarr, 2) = delay(iarr, 1, 1, it_2);	% time 2

  % loop over the end points of this bh_time interval

  ts_value = 0.0;

  for inode = 1:2,

    % delay, sum each of the arrivals for this vertex of the cube

    iarr = 1:narrivals(inode);

    % compute the time when these arrivals left the source

    t_src = t_rcv(isample) - delayi(iarr, inode);

    % compute the associated sample indices in the upsampled source waveform

    isamples = round(t_src * src_usfs) + 1;	% Matlab arrays start at 1

    indx = find((isamples > 0) & (isamples <= nsamples_src));

    % compute the contributions from these arrivals to the sample value

    ts_values = real(ampi(iarr(indx), inode)) .* hilb_uswf_r(isamples(indx)) ...
              - imag(ampi(iarr(indx), inode)) .* hilb_uswf_i(isamples(indx));

    % sum all of the contributions from these arrivals

    ts_value = ts_value + sum(ts_values);

    % next vertex of the cube

  end		% for inode = 1:2

  % saved the computed time series sample value

  rcv_ts(isample) = ts_value;

  % next sample ...

end		% for isample = 1:nsamples

else		% (nrr == 1 && nrd == 1) is TRUE

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Case of a (uniformly spaced) grid of receivers (at least two ranges, depths)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% verify that there are at least two ranges, and at least two depths

if (nrr > 1),
  dr = (rr(end) - rr(1)) / (nrr - 1);
else
  error('moving_environment: receiver grid must have at least two ranges');
end

if (nrd > 1),
  dz = (rd(end) - rd(1)) / (nrd - 1);
else
  error('moving_environment: receiver grid must have at least two depths');
end

% for a grid of receivers, the upper limit of the cubes are open intervals [L,H)

if max(rcv_pos(:,1)) >= max(rr),
  error('moving_environment: track range moves outside of BELLHOP range grid')
end

if max(rcv_pos(:,2)) >= max(rd),
  error('moving_environment: track depth moves outside of BELLHOP depth grid')
end

% allocate memory for the arrival amp, delay, and (cube) ray coordinates

maxarr = max(max(max( narrmat )));

ampi   = zeros(maxarr, 8);
delayi = zeros(maxarr, 8);
ray_r  = zeros(maxarr, 8);
ray_d  = zeros(maxarr, 8);

% loop over the time values where the received time series is to be computed

for isample = 1:nsamples,

  % find indices of the (range, depth, bh_time) cube that contains the receiver

  ir_1 = 1 + floor((r_rcv(isample) - rr(1)) / dr);
  id_1 = 1 + floor((d_rcv(isample) - rd(1)) / dz);
  it_1 = 1 + floor((t_rcv(isample) - bh_time(1)) / dbh_time);

  ir_2 = ir_1 + 1;
  id_2 = id_1 + 1;
  it_2 = it_1 + 1;

  % determine the maximum number of arrivals (over all 8 corners of the cube)

  narrivals = [ narrmat(ir_1, id_1, it_1), ...
                narrmat(ir_2, id_1, it_1), ...
                narrmat(ir_2, id_2, it_1), ...
                narrmat(ir_1, id_2, it_1), ...
                narrmat(ir_1, id_1, it_2), ...
                narrmat(ir_2, id_1, it_2), ...
                narrmat(ir_2, id_2, it_2), ...
                narrmat(ir_1, id_2, it_2) ];

  maxarr = max(narrivals);

  iarr = 1:maxarr;

  % compute the weights for interpolation in range, depth, and bh_time

  rwgt = (r_rcv(isample) - rr(ir_1)) / dr;
  dwgt = (d_rcv(isample) - rd(id_1)) / dz;
  twgt = (t_rcv(isample) - bh_time(it_1)) / dbh_time;

  % compute the (trilinear interpolation) ray amplitude weights

  ampi(iarr, 1) = (1-rwgt)*(1-dwgt)*(1-twgt)*amp(iarr, ir_1, id_1, it_1);
  ampi(iarr, 2) = (  rwgt)*(1-dwgt)*(1-twgt)*amp(iarr, ir_2, id_1, it_1);
  ampi(iarr, 3) = (  rwgt)*(  dwgt)*(1-twgt)*amp(iarr, ir_2, id_2, it_1);
  ampi(iarr, 4) = (1-rwgt)*(  dwgt)*(1-twgt)*amp(iarr, ir_1, id_2, it_1);

  ampi(iarr, 5) = (1-rwgt)*(1-dwgt)*(  twgt)*amp(iarr, ir_1, id_1, it_2);
  ampi(iarr, 6) = (  rwgt)*(1-dwgt)*(  twgt)*amp(iarr, ir_2, id_1, it_2);
  ampi(iarr, 7) = (  rwgt)*(  dwgt)*(  twgt)*amp(iarr, ir_2, id_2, it_2);
  ampi(iarr, 8) = (1-rwgt)*(  dwgt)*(  twgt)*amp(iarr, ir_1, id_2, it_2);

  % get the BELLHOP time delays that correspond to each corner of the cube

  delayi(iarr, 1) = delay(iarr, ir_1, id_1, it_1);	% time 1
  delayi(iarr, 2) = delay(iarr, ir_2, id_1, it_1);
  delayi(iarr, 3) = delay(iarr, ir_2, id_2, it_1);
  delayi(iarr, 4) = delay(iarr, ir_1, id_2, it_1);

  delayi(iarr, 5) = delay(iarr, ir_1, id_1, it_2);	% time 2
  delayi(iarr, 6) = delay(iarr, ir_2, id_1, it_2);
  delayi(iarr, 7) = delay(iarr, ir_2, id_2, it_2);
  delayi(iarr, 8) = delay(iarr, ir_1, id_2, it_2);

  % find the range, depth coordinates relative to the associated cube vertex

  rcv_cube_r(1) = r_rcv(isample) - rr(ir_1);	% time 1
  rcv_cube_r(2) = r_rcv(isample) - rr(ir_2);
  rcv_cube_r(3) = r_rcv(isample) - rr(ir_2);
  rcv_cube_r(4) = r_rcv(isample) - rr(ir_1);

  rcv_cube_r(5) = rcv_cube_r(1);		% time 2 is the same
  rcv_cube_r(6) = rcv_cube_r(2);
  rcv_cube_r(7) = rcv_cube_r(3);
  rcv_cube_r(8) = rcv_cube_r(4);

  rcv_cube_d(1) = d_rcv(isample) - rd(id_1);	% time 1
  rcv_cube_d(2) = d_rcv(isample) - rd(id_1);
  rcv_cube_d(3) = d_rcv(isample) - rd(id_2);
  rcv_cube_d(4) = d_rcv(isample) - rd(id_2);

  rcv_cube_d(5) = rcv_cube_d(1);		% time 2 is the same
  rcv_cube_d(6) = rcv_cube_d(2);
  rcv_cube_d(7) = rcv_cube_d(3);
  rcv_cube_d(8) = rcv_cube_d(4);

  % get the components of the tangent vectors to the rays at the grid points

  ray_r(iarr, 1) = cos( rcv_angle(iarr, ir_1, id_1, it_1) ); % time 1
  ray_r(iarr, 2) = cos( rcv_angle(iarr, ir_2, id_1, it_1) );
  ray_r(iarr, 3) = cos( rcv_angle(iarr, ir_2, id_2, it_1) );
  ray_r(iarr, 4) = cos( rcv_angle(iarr, ir_1, id_2, it_1) );

  ray_r(iarr, 5) = cos( rcv_angle(iarr, ir_1, id_1, it_2) ); % time 2
  ray_r(iarr, 6) = cos( rcv_angle(iarr, ir_2, id_1, it_2) );
  ray_r(iarr, 7) = cos( rcv_angle(iarr, ir_2, id_2, it_2) );
  ray_r(iarr, 8) = cos( rcv_angle(iarr, ir_1, id_2, it_2) );

  ray_d(iarr, 1) = sin( rcv_angle(iarr, ir_1, id_1, it_1) ); % time 1
  ray_d(iarr, 2) = sin( rcv_angle(iarr, ir_2, id_1, it_1) );
  ray_d(iarr, 3) = sin( rcv_angle(iarr, ir_2, id_2, it_1) );
  ray_d(iarr, 4) = sin( rcv_angle(iarr, ir_1, id_2, it_1) );

  ray_d(iarr, 5) = sin( rcv_angle(iarr, ir_1, id_1, it_2) ); % time 2
  ray_d(iarr, 6) = sin( rcv_angle(iarr, ir_2, id_1, it_2) );
  ray_d(iarr, 7) = sin( rcv_angle(iarr, ir_2, id_2, it_2) );
  ray_d(iarr, 8) = sin( rcv_angle(iarr, ir_1, id_2, it_2) );

  % loop over the 8 nodes, vertices of the (range, depth, bh_time) cube

  ts_value = 0.0;

  for inode = 1:8,

    % compute the time delays from this cube corner to current receiver position

    pos_delays = ( rcv_cube_r(inode) * ray_r(:, inode) ...
               +   rcv_cube_d(inode) * ray_d(:, inode) ) / c0;

    % delay, sum each of the arrivals for this vertex of the cube

    iarr = 1:narrivals(inode);

    % compute the time when these arrivals left the source

    t_src = t_rcv(isample) - (delayi(iarr, inode) + pos_delays(iarr));

    % compute the associated sample indices in the upsampled source waveform

    isamples = round(t_src * src_usfs) + 1;	% Matlab arrays start at 1

    indx = find((isamples > 0) & (isamples <= nsamples_src));

    % compute the contributions from these arrivals to the sample value

    ts_values = real(ampi(iarr(indx), inode)) .* hilb_uswf_r(isamples(indx)) ...
              - imag(ampi(iarr(indx), inode)) .* hilb_uswf_i(isamples(indx));

    % sum all of the contributions from these arrivals

    ts_value = ts_value + sum(ts_values);

    % next vertex of the cube

  end		% for inode = 1:8

  % saved the computed time series sample value

  rcv_ts(isample) = ts_value;

  % next sample ...

end		% for isample = 1:nsamples

end		% (nrr == 1 && nrd == 1) is FALSE

% that's all folks !!!

return;

% last line of moving_environment.m
