% HLS VIRTEX Toolbox 
%
% VirTEX for Platform Motion (VirTEX Extra Lite)
% delay_sum_lib      - pre-compute "dopplerized" source waveforms for delay_sum
% delay_sum          - compute receiver timeseries (for steady src/rcv motion)
%
% VirTEX for Sea Surface Motion (VirTEX Lite)
% delay_sum_surface  - compute receiver timeseries (with surface motion)
%
% VirTEX (Full VirTEX)
% Bellhop_inputs     - construct input files for BELLHOP for given parameters
% Bellhop_execute    - execute a list of BELLHOP runs (in parallel, UNIX only)
% moving_environment - compute receiver timeseries (VIRTEX core)
%
% Helper functions
% fourier_upsample   - upsample a timeseries (zero pads in frequency domain)
% load_arrfils       - load and collate a list of BELLHOP binary arrival files
% read_arrivals_bin  - read the contents of a BELLHOP binary arrival file
% read_rayfil        - read the contents of a BELLHOP eigenray data file
% surface_height     - compute space-time snapshots of sea surface waves
% wafo_wrapper       - wrapper function to access selected WAFO functions
%
% $Id: Contents.m,v 1.2 2011/07/06 19:03:51 jcp Exp $
