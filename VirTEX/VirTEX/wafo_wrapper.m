function [surf_t, surf_r, surf_h, sample_tf, sample_sf] = wafo_wrapper(params, Tsurf, Rsurf);

% WAFO_WRAPPER - interface WAFO wave spectra functions to VIRTEX driver
%
%   wafo-jonswap        JONSWAP spectral density
%   wafo-mccormick      McCormick spectral density
%   wafo-ohspec         Ochi-Hubble spectral density
%   wafo-ohspec2        Ochi-Hubble spectral density (bimodal)
%   wafo-ohspec3        Ochi-Hubble spectral density (mixed sea state)
%   wafo-oscspec        Harmonic Oscillator spectral density
%   wafo-pmspec         Pierson-Moskowitz spectral density
%   wafo-tmaspec        JONSWAP spectral density (finite water depth)
%   wafo-torsethaugen   Double peaked spectral density (swell + wind)
%   wafo-wallop         Wallop spectral density
%
% $Id: wafo_wrapper.m,v 1.2 2011/06/13 18:40:04 jcp Exp $

% verify the WAFO toolbox has been initialized before proceeding

which_initwafo = which('initwafo(''minimum'',1)');

if which_initwafo(end-length('initwafo.m')+1:end) ~= 'initwafo.m',
  error('wafo_wrapper: you must run initwafo before using WAFO sea spectra');
end;

% extract surface spectrum type field from params (required)

if isfield(params, 'surface_type'),
  surface_type = params.surface_type;
else
  error('wafo_wrapper: required field surface_type is missing from params');
end;

% extract optional fields from params that are common to all models

if isfield(params, 'cutoff_freq'),
  cutoff_freq = params.cutoff_freq;
else
  cutoff_freq = 1.0/2.0;	% default cutoff frequency of spectrum (0.5 Hz)
end;

if isfield(params, 'height'),
  Hm0 = params.height;
else
  Hm0 = 1.0;	% default significant wave height (1 meter peak to trough)
end;

if isfield(params, 'os_time'),
  os_time = params.os_time;
  if os_time < 1,
    error('surface_height: temporal over sampling factor cannot be < 1');
  end;
else
  os_time = 1.0;
end;

if isfield(params, 'os_space'),
  os_space = params.os_space;
  if os_space < 1,
    error('surface_height: spatial over sampling factor cannot be < 1');
  end;
else
  os_space = 1.0;
end;

if isfield(params, 'peak_period'),
  Tp = params.peak_period;
else
  Tp = 10.0;	% default wave period associated with spectral peak (10 seconds)
end;

if isfield(params, 'seed'),
  seed = params.seed;
  s = RandStream.create('mt19937ar', 'seed', seed);
  RandStream.setDefaultStream(s);
end;

% Compute the Nyquist spatial and temporal sample rates

g = gravity();

sample_tf = 2.0 * cutoff_freq;				% cycles / second

sample_sf = 2.0 * (2.0 * pi * cutoff_freq^2 / g);	% cycles / meter

% Compute the oversampled spatial and temporal sample rates

sample_tf = sample_tf * os_time;

sample_sf = sample_sf * os_space;

% Compute the sample period, sample length

dt = 1.0 / sample_tf;

dx = 1.0 / sample_sf;

% Compute the transform lengths in space and time

Nt = 1 + ceil(sample_tf*Tsurf);

Nx = 2 ^ nextpow2( 1 + ceil(sample_sf*Rsurf) );		% seasim uses FFTs

% Compute a one dimensional "frequency' wave spectrum

plotflag = 0;		% please - no obnoxious plots !!!

omega = 2.0 * pi * linspace(0.0, cutoff_freq, 257)';

switch lower(surface_type)

  case 'wafo-jonswap',

    S1d = jonswap(omega, [Hm0 Tp], plotflag);

  case 'wafo-mccormick',

    S1d = mccormick(omega, [Hm0 Tp], plotflag);

  case 'wafo-ohspec',

    S1d = ohspec(omega, [Hm0 Tp], plotflag);

  case 'wafo-ohspec2',

    S1d = ohspec2(omega, [Hm0 1], plotflag);

  case 'wafo-ohspec3',

    error(['wafo_wrapper: ', surface_type, 'is not implemented yet']);

  case 'wafo-oscspec',

    error(['wafo_wrapper: ', surface_type, 'is not implemented yet']);

  case 'wafo-pmspec',

    S1d = pmspec(omega, [Hm0 Tp], plotflag);

  case 'wafo-tmaspec',

    error(['wafo_wrapper: ', surface_type, 'is not implemented yet']);

  case 'wafo-torsethaugen',

    S1d = torsethaugen(omega, [Hm0 Tp], plotflag);

  case 'wafo-wallop',

    S1d = wallop(omega, [Hm0 Tp], plotflag);

  otherwise

    error(['wafo_wrapper: unknown surface model type: ', surface_type]);

end;	% switch lower(surface_type)

% Compute a directional wave spectrum

Dspread = spreading([], 'cos2s',[] ,[], omega, []);	% spreading function

S2d = mkdspec(S1d, Dspread, plotflag);	% multiply by spreading function

% Simulate the sea surface

fftdim = 2;

Y = seasim(S2d, Nx, 0, Nt, dx, dx, dt, fftdim, plotflag);

% Extract what we need from the output struct Y

surf_t = Y.t;
surf_r = Y.x;
surf_h = Y.Z;

% That's all folks...

return;

