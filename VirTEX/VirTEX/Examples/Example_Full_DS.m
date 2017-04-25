% Example_Full_Down_Swell
%
% Example script to run Full VirTEX for down swell case
%
% $Id: Example_Full_DS.m,v 1.2 2011/07/08 22:43:17 jcp Exp $

% configure low and high "water marks" for execution of BELLHOP in parallel
jobLoLevel = 12;
jobHiLevel = 16;

% Source waveform time series parameters
src_wf_filename = 'cw_5kHz_1s.wav';
src_upsample_ratio = 16;

% Receive waveform time series parameters
chan_T = 0.0;	% let VirTEX figure it out
rcv_fs = 48000;

% Construct spatial-temporal snapshots of sea surface
Tsurf = 6.5;		% end time for surface (seconds)
Rsurf = 2000.0;		% max range for surface (meters)

swell_amp = 5.0;	% caution: not the same as significant wave height
swell_period = 8.0;
no_range_pnts = 1024;
no_time_pnts = 256;
surf_r = Rsurf*(0:no_range_pnts-1)/(no_range_pnts-1);
surf_t = Tsurf*(0:no_time_pnts-1)/(no_time_pnts-1);
surf_h = zeros(no_range_pnts, no_time_pnts);

grav = 9.81;
omega = 2.0*pi/swell_period;
k = omega^2 / grav;

for j = 1:no_range_pnts,
  % generate a simple swell wave
  surf_h(j,:) = swell_amp * sin(omega*surf_t - k*surf_r(j));
end;

% source depth (meters)
src_depth = 50.0;

% receiver motion array; range(meters), depth(meters), time(sec)
rcv_pos = [1500.0, 50.0,   0.0; ...
           1500.0, 50.0, Tsurf;];

% spacing of range, depth grid in BELLHOP calculations
delta_r = 5.0;		% range spacing of BELLHOP grid (meters)
delta_d = 5.0;		% depth spacing of BELLHOP grid (meters)

% configure all of the BELLHOP run parameters by populating the bparams struct

bparams.bellhop_path = '/usr/local/bin/bellhop.exe';	% you will VERY LIKELY need to edit this!
bparams.calc_top_dir = 'BellhopRuns';
bparams.calc_sub_dir = 'frame_%04d';
bparams.arrfil_fmt = 'BELLHOP_%04d.arr';

bparams.template_env = 'Example_Full.env';
bparams.btyfil = 'Example_Full.bty';
bparams.freqs = 5000.0;
bparams.nbeams = 1000;
bparams.angles = [-85 85];
bparams.runtype = 'aB R';

% Now construct the source waveform. The VIRTEX processor will assume
% that the source waveform amplitudes vanish for time values less than
% zero or greater than the duration of the waveform. Hence it is NOT
% necessary to surround the waveform with any kind of zero padding.

% Read the discrete sampled source waveform (pcm/wav format)
[src_wf, src_fs, src_nbits] = wavread(src_wf_filename);

% Upsample the waveform by zero padding in the frequency domain
src_uswf = fourier_upsample(src_wf, src_upsample_ratio);
src_usfs = src_upsample_ratio * src_fs;

% Compute the Hilbert transform of the source waveform
src_huswf = hilbert(src_uswf);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Run BELLHOP to compute the ray arrivals at each temporal time step          %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

disp('Working...');

% verify the requested time is within the time spanned by the BELLHOP runs

Tmax = (length(src_wf)/src_fs) + chan_T;

if (Tmax >= Tsurf),
  error('Tmax is outside of the time interval spanned by the BELLHOP runs!');
end;

% check for receiver motion and setup an appropriate BELLHOP range, depth grid

if ( max(rcv_pos(:,1)) > min(rcv_pos(:,1)) ...
  || max(rcv_pos(:,2)) > min(rcv_pos(:,2)) ),

% receiver is moving, at least 2 ranges, and 2 depths are required in the grid

ngrid_pnts_range = 2 + floor((max(rcv_pos(:,1)) - min(rcv_pos(:,1))) / delta_r);
ngrid_pnts_depth = 2 + floor((max(rcv_pos(:,2)) - min(rcv_pos(:,2))) / delta_d);

grid_ranges = (0:ngrid_pnts_range-1) * delta_r;
grid_depths = (0:ngrid_pnts_depth-1) * delta_d;

grid_ranges = mean(rcv_pos(:,1)) + grid_ranges - mean(grid_ranges);
grid_depths = mean(rcv_pos(:,2)) + grid_depths - mean(grid_depths);

else

% no receiver motion, only one point is needed

grid_ranges = rcv_pos(1,1);
grid_depths = rcv_pos(1,2);

end;

% finish populating the bparams struct with the BELLHOP run parameters

bparams.src_depth = src_depth;
bparams.rec_ranges = grid_ranges/1000.0; % two or more receiver ranges
bparams.rec_depths = grid_depths; % two or more receiver depths (to form a box)

% create ALL of the BELLHOP input files

disp('Generating BELLHOP input files ...');

Bellhop_inputs(bparams, surf_t, surf_r, surf_h);

% execute ALL of the BELLHOP input files

disp('Computing BELLHOP phase delays over range, depth, time grid');

ntime = length(surf_t);

Bellhop_execute(bparams, ntime, jobLoLevel, jobHiLevel);

% load the BELLHOP arrival files
clear arrival_params
arrival_params = load_arrfils('BELLHOP_%04d.arr', surf_t);

% compute the received time series using Full VirTEX
disp('Convolving BELLHOP phase delays with the source waveform');
[rcv_ts, Tstart] = moving_environment(rcv_pos, ...
         src_huswf, src_usfs, arrival_params, rcv_fs, chan_T);

% compute the sonogram, spectrogram and plot the results
npts_wind = 8192;
npts_overlap = 7680;
npts_fft = 16384;
spectrogram(rcv_ts, npts_wind, npts_overlap, npts_fft, rcv_fs);
set(gca, 'XLim', [4900 5100]);
set(gca, 'Fontweight', 'bold');
set(gca, 'Tickdir', 'out');
colorbar();

%
% End of Example_Full_Down_Swell.m
