% Example_Lite_CS
%
% Example script to run VirTEX for Sea Surface Motion for cross swell case
%
% $Id: Example_Lite_CS.m,v 1.6 2011/07/08 19:04:32 jcp Exp $

% Source time series parameters
src_ts_filename = 'cw_5kHz_1s.wav';	% 5 kHz CW tone
src_upsample_ratio = 16;	% controls interpolation accuracy

% Receive time series parameters
chan_T = 0.0;	% let VirTEX figure it out
rcv_fs = 48000;

% BELLHOP output filenames
ARRFIL = 'Example_Lite.arr';
RAYFIL = 'Example_Lite.ray';

% Parameters for spatial-temporal sea surface model
tau_s = 8.0;		% period of swell wave (seconds)
sparams.h0 = 10.0;	% wave amplitude (not same as significant wave height)
sparams.k = 0.0;	% cross swell case has no spatial range dependence
sparams.omega = 2.0*pi/tau_s;	% temporal angular frequency
sparams.phi_0 = 0.0;	% initial phase angle

% Initialize the function that computes the surface wave height (this
% example uses the option to pass a function handle to delay_sum_surface
% that computes the surface heights "on the fly" instead of passing arrays
% with precomputed heights).
gravity_wave([], [], sparams);

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

% Compute the received time series using VirTEX for Sea Surface Motion

[rcv_ts, tstart] = ...
   delay_sum_surface(ARRFIL, RAYFIL, ...
      src_husts, src_usfs, [], [], @gravity_wave, rcv_fs, [], [], chan_T);

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
% End of Example_Lite_CS.m
