function src_wf = delay_sum_lib(wf_filename, v_min, v_max, numvel, c0)

% DELAY_SUM_LIB - construct library of Dopplerized src waveforms for delay_sum
%
% Usage: src_wf = delay_sum_lib(wf_filename, v_min, v_max, numvel, c0)
%
%   wf_filename - filename of the source waveform (PCM/WAV format)
%         v_min - minimum signed velocity (meters/second)
%         v_max - maximum signed velocity (meters/second)
%        numvel - number of velocity bins (integer)
%            c0 - reference sound speed to convert velocity to Doppler
%
%        src_wf - struct containing the library of dopplerized src waveforms
%
% NOTES:
%
% 1) The source waveform should have a little bit of zero padding both
% before and after the actual waveform. If the non-zero portion of the
% waveform consists of N samples, then there should be *at least*
% N*v_max/c0 additional samples at each end. This is to insure that
% the stretched doppler replicas do not get trucated.
%
% 2) Be sure that the PCM/WAV file has the proper sampling rate.
%
% 3) The arbitrary_alpha() function from the Acoustics Toolbox is used
% to compute the "stretched" doppler replicas. It internally calls the
% standard Matlab fftfilt() routine. Unfortunately, fftfilt has a quirk
% which can result in *very* long run times when the number of samples
% in the source waveform exceeds roughly 2^18 = 262144. We hope to fix
% this issue at a future date.
%
% $Id: delay_sum_lib.m,v 1.2 2011/05/18 19:25:28 jcp Exp $
% 23-Apr-17 1:09PM changed out of wavread to audioread 
% read the source timeseries (assumed to be a PCM/WAV format file for now)

[sts, fs] = audioread(wf_filename);

% normalize the source timeseries so that it has 0 dB (peak) source level

sts = sts(:) / max(abs(sts));

% pre-calculate Dopplerized versions of the source timeseries

nsamples_src = length(sts);

sts_hilbert = zeros(nsamples_src, numvel);

alpha_vec = 1.0 - linspace(v_min, v_max, numvel) / c0;

for ialpha = 1:numvel,
  % resample to get the steady motion stretched waveform
  sts_strectched = arbitrary_alpha(sts.', alpha_vec(ialpha));
  % compute the hilbert transform of the stretched waveform
  sts_hilbert(:, ialpha) = hilbert(sts_strectched(1:nsamples_src).');
end

% pack everything into a single structure for convenience

src_wf.c0 = c0;
src_wf.fs = fs;
src_wf.sts = sts;
src_wf.alpha_vec = alpha_vec;
src_wf.sts_hilbert = sts_hilbert;

% all done

return;

%
% end of delay_sum_lib.m
