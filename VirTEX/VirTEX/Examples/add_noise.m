%
% add_noise - assemble received simulation packets and add ambient noise
%
% $Id: add_noise.m,v 1.3 2011/05/23 03:51:49 jcp Exp $

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Configure the case study "matrix" parameters
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Number of "top level" cases in the study (e.g. surface wave heights)
Ncase = 40;

% Number of packets to transmit for each case (add_noise will add 1 sync pkt)
Nxmit = 45;

% Time duration of transmitted packets (zero to make it as short as possible)
xmit_T = 6.0;

% Format for the receive time series saveset filenames (icase, ixmit)
rcv_saveset_fmt = 'ts_case_%03d_%03d.mat';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Configure the output PCM file parameters
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

pcm_filename = '2011-05-23-A.pcm';

max_int = 2^15 - 2;	% abs value of maximum integer value for PCM file

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Configure the ambient noise parameters
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

AN_db = 35.0;	% ambient noise level (power spectral density, dB)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Configure the source level and modem packet parameters
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

SL_db = 185.0;	% source level (total power, dB)

src_fc = 10000.0;	% modem packet center frequency (Hz)

src_BW = 10000.0;	% modem packet bandwidth (Hz)

swf_saveset_name = 'src_wf.mat';	% source waveform saveset name

fade_beg = 0.025;	% duration of tukey transition at waveform start (sec)
fade_end = 0.025;	% duration of tukey transition at waveform end   (sec)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Nothing more to edit below this line ...
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

i = sqrt(-1.0d0);	% everyone knows "j" is an index

pi = acos(-1.0d0);

% Read the source time series saveset

saveset = load(swf_saveset_name, '-mat');
nsamples_src = length(saveset.src_wf);
src_wf = saveset.src_wf;
src_fs = saveset.src_fs;

% Read the first receive time series saveset

rcv_saveset_name = sprintf(rcv_saveset_fmt, 1, 1);
saveset = load(rcv_saveset_name, '-mat');
nsamples_rcv = length(saveset.rcv_ts);
rcv_fs = saveset.rcv_fs;
clear saveset;

% Down sample the source waveform to the receive sample rate (if needed)

if src_fs > rcv_fs,
  src_wf = resample(src_wf, rcv_fs, src_fs);
  nsamples_src = length(src_wf);
end;

% Compute the number of samples in each transmission

if xmit_T > 0.0,

  % transmission length was specified by the user

  nsamples_xmit = round( xmit_T * rcv_fs );

  if ((nsamples_xmit < nsamples_rcv) || (nsamples_xmit < nsamples_src)),
    error( [mfilename, ': xmit_T is shorter than src or rcv waveform'] );
  end;

else

  % compute the smallest possible transmission length (may cause modem glitches)

  if nsamples_rcv < nsamples_src,
    error( [mfilename, ': nsamples_rcv < nsamples_src is not supported'] );
  end;

  nsamples_xmit = max( [nsamples_src, nsamples_rcv] );

end;

% Construct the transition regions of the time domain tukey taper

nsamples_beg = round(fade_beg * rcv_fs);
nsamples_end = round(fade_end * rcv_fs);

tukey_beg = [ sin( 0.5 * pi * (0:(nsamples_beg-1)).'/nsamples_beg ) ].^2;
tukey_end = [ sin( 0.5 * pi * (1:(nsamples_end  )).'/nsamples_end ) ].^2;
tukey_end = flipud( tukey_end );

% Parameters for construction of ambient (white) noise

SL_amp = 10.0 ^ (SL_db/20.0);
AN_amp = 10.0 ^ (AN_db/20.0);

freq_max = src_fc + 0.5*src_BW;
for j = floor(rcv_fs/(2.0*freq_max)):-1:1,
  if (rem(rcv_fs,j) == 0),
    us_ratio = j;
    break;
  end;
end;

nsamples_row = (Nxmit+1) * nsamples_xmit;

nsamples_amb = nsamples_row / us_ratio;

% Open the output PCM file

fid = fopen(pcm_filename, 'w', 'ieee-le');

% Display a summary of the processing to stdout

total_duration = Ncase*(Nxmit + 1)*nsamples_xmit/rcv_fs;

fprintf('Output filename ............................ %s\n', pcm_filename);
fprintf('Sample rate (Hz) ........................... %.2f\n', rcv_fs);
fprintf('Total duration (sec) ....................... %.2f\n', total_duration);
fprintf('Number of Packets (including sync) ......... %.2f\n', (Nxmit+1)*Ncase);
fprintf('Source level (total power, dB) ............. %.2f\n', SL_db);
fprintf('Noise level (power spectral density, dB) ... %.2f\n', AN_db);
fprintf('Noise cutoff frequency (Hz) ................ %.2f\n', 0.5*rcv_fs/us_ratio);

% Loop over cases

for icase = 1:Ncase,

  % Storage for one "row" of signal

  ts_sig = zeros(nsamples_row, 1);

  % Loop over transmissions

  max_amp = 0.0;

  ixmit_offset = nsamples_xmit;

  for ixmit = 1:Nxmit,

    % load this receive saveset

    rcv_saveset_name = sprintf(rcv_saveset_fmt, icase, ixmit);

    saveset = load(rcv_saveset_name, '-mat');

    rcv_ts = saveset.rcv_ts;

    % apply the time domain tukey taper

    rcv_ts(1:nsamples_beg) = tukey_beg.*rcv_ts(1:nsamples_beg);

    rcv_ts(end-nsamples_end+1:end) = tukey_end.*rcv_ts(end-nsamples_end+1:end);

    % update the maximum amplitude (used later for scaling the reference packet)

    max_amp = max( [max_amp, max(abs(rcv_ts))] );

    % store this packet in the "row" vector

    ts_sig(ixmit_offset+1:ixmit_offset+nsamples_rcv) = rcv_ts;

    % next packet

    ixmit_offset = ixmit_offset + nsamples_xmit;

  end;

  % scale and insert the un-distorted, reference packet for this "row"

  ts_sig(1:nsamples_src) = 0.5 * max_amp * src_wf;	% 1-st packet of row

  % generate band limited white Gaussian noise

  ts_amb = fourier_upsample(randn(nsamples_amb, 1), us_ratio);

  ts_amb = sqrt(2.0*rcv_fs/us_ratio) * AN_amp * ts_amb;

  % scale the signal by the source level and add in the white noise

  ts_sig = (SL_amp * ts_sig) + ts_amb;

  % convert to signed 16 bit integer, and write out this "row" to the PCM file

  fwrite(fid, round((max_int / max(abs(ts_sig))) * ts_sig), 'int16');

end;

fclose(fid);

%
% end of add_noise.m
