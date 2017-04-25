%
% Example script to run VirTEX for Platform Motion for moving receiver
%
% $Id: Example_MP.m,v 1.1 2011/07/08 18:43:39 jcp Exp $
% 23-Apr-17 1:09PM E. Oon: this is very outdated and has functions that no
% longer work in MATLAB 
% Identify the BELLHOP arrivals file (generate it by running BELLHOP)

ARRFIL = 'Example_MP';

% Identify the source timeseries/waveform (5 kHz CW tone, 1 sec duration)

wf_file = 'cw_5kHz_1s.wav';

% Generate the library of pre-stretched source waveforms. Note that in
% practice this can be done once, saved as a MATLAB saveset, and read
% from disk storage when needed.

v_min = -1.0;		% minimum net velocity
v_max =  1.0;		% maximum net velocity
numvel = 51;		% number of velocity "bins"
c0 = 1539.0;

src_wf = delay_sum_lib(wf_file, v_min, v_max, numvel, c0);

% Compute the timeseries observed at the receivers defined in the BELLHOP
% arrivals file. The output is a 3 dimensional array (number of samples by
% number of rcv ranges by number of rcv depths). The sample rate of the rcv
% timeseries is identical to the source.

v_src = [ 0.000, 0.000 ];		% source position is fixed
v_rcv = [ 0.866, 0.500 ];		% receiver velocity vector (m/s)

rcv_wf = delay_sum(ARRFIL, src_wf, v_src, v_rcv);

% End of Example_MP.m
