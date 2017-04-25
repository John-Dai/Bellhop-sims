% Example_Lite_DS
%
% Example script to run VirTEX for Sea Surface Motion for down swell case
%
% NOTE: to run this example, the WAFO toolkit must be installed
%
% $Id: Example_Lite_DS.m,v 1.4 2011/07/08 19:16:09 jcp Exp $

% Source waveform time series parameters
src_ts_filename = 'cw_5kHz_1s.wav';	% 5 kHz CW tone
src_upsample_ratio = 16;

% Receive waveform time series parameters
chan_T = 0.0;	% let VirTEX figure it out
rcv_fs = 48000;

% BELLHOP output filenames
ARRFIL = 'Example_Lite.arr';
RAYFIL = 'Example_Lite.ray';

% Now construct the source waveform. VIRTEX-lite (delay_sum_surface) will
% assume that the source waveform amplitudes vanish for time values less
% than zero or greater than the duration of the waveform. Hence it is NOT
% necessary to surround the waveform with any kind of zero padding.

% Read the discrete sampled source waveform (pcm/wav format)
[src_ts, src_fs, src_nbits] = wavread(src_ts_filename);

% Upsample the waveform by zero padding in the frequency domain
src_usts = fourier_upsample(src_ts, src_upsample_ratio);
src_usfs = src_upsample_ratio * src_fs;

% Compute the Hilbert transform of the source waveform
src_husts = hilbert(src_usts);

% Parameters for spatial-temporal sea surface model
Tsurf = 5.0;		% end time for surface (seconds)
Rsurf = 1500.0;		% max range for surface (meters)

sparams.surface_type = 'wafo-pmspec';
sparams.cutoff_freq = 0.40;
sparams.peak_period = 8.0;
sparams.height = 10.0;
sparams.os_time = 32;
sparams.os_space = 4;

% Construct spatial-temporal snapshots of sea surface height using WAFO toolbox
[surf_t, surf_r, surf_h, sample_tf, sample_sf] = ...
                                    surface_height(sparams, Tsurf, Rsurf);

% Compute the received time series using VirTEX for Sea Surface Motion
[rcv_ts, tstart] = ...
   delay_sum_surface(ARRFIL, RAYFIL, ...
      src_husts, src_usfs, surf_t, surf_r, surf_h, rcv_fs, [], [], chan_T);

% Compute the sonogram, spectrogram and plot the results
npts_wind = 8192;
npts_overlap = 7680;
npts_fft = 16384;
spectrogram(rcv_ts, npts_wind, npts_overlap, npts_fft, rcv_fs);
set(gca, 'XLim', [4900 5100]);
set(gca, 'Fontweight', 'bold');
set(gca, 'Tickdir', 'out');
colorbar();

%
% End of Example_Lite_DS.m
