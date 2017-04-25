function [surf_t, surf_r, surf_h, sample_tf, sample_sf] = surface_height(params, Tsurf, Rsurf)

% SURFACE_HEIGHT - construct space-time snapshots of sea surface
%
% $Id: surface_height.m,v 1.2 2010/05/28 21:06:05 jcp Exp $

% Variable glossary:
%
% -----------------------------------+-------------------------------------+
% Temporal variables                 | Spatial variables                   |
% -----------------------------------+-------------------------------------+
% tau        wave period (sec)       | lambda     wave length (meters)     |
% tfreq      frequency (1/sec)       | sfreq      frequency (1/meters)     |
% omega      angular frequency       | k          wave number              |
% sample_tf  temporal sampling rate  | sample_sf  spatial sampling rate    |
% -----------------------------------+-------------------------------------+

% declare mathematical and physical constants

pi = acos(-1.0);

g = 9.81;	% acceleration of gravity (meters/second^2)

% extract surface spectrum type field from params (required)

if isfield(params, 'surface_type'),
  surface_type = params.surface_type;
else
  error('surface_height: required field surface_type is missing from params');
end;

% for WAFO sea surface spectra, call the dedicated wrapper function

if strncmp(lower(surface_type), 'wafo-', 5),
  % call the WAFO interface wrapper function to handle everything ...
  [ surf_t, surf_r, surf_h, sample_tf, sample_sf ] = ...
                                      wafo_wrapper(params, Tsurf, Rsurf);
  % that's all folks !!!
  return;
end;

% extract optional fields from params that are common to remaining models

if isfield(params, 'os_time'),
  os_time = params.os_time;
  if os_time < 1,
    error('surface_height: temporal over sampling factor cannot be < 1');
  end;
else
  os_time = 4.0;
end;

if isfield(params, 'os_space'),
  os_space = params.os_space;
  if os_space < 1,
    error('surface_height: spatial over sampling factor cannot be < 1');
  end;
else
  os_space = 4.0;
end;

if isfield(params, 'seed'),
  seed = params.seed;
  s = RandStream.create('mt19937ar', 'seed', seed);
  RandStream.setDefaultStream(s);
end;

% perform the model dependent calculations

switch lower(surface_type)

  case 'flat'

    % period (seconds)
    if isfield(params, 'peak_period'),
      tau = params.peak_period;
    else
      error('surface_height: required field peak_period is missing from params');
    end;
    % compute temporal sample rate
    tfreq = 1.0 / tau;
    sample_tf = 2.0 * tfreq * os_time;
    % compute temporal sample points
    npnts_t = 1 + ceil(sample_tf*Tsurf);
    surf_t = (0:(npnts_t-1)) / sample_tf;
    % compute spatial sample rate (no need to worry about aliasing here)
    sfreq = 1.0 / 500.0;
    sample_sf = sfreq;
    % compute spatial sample points
    npnts_r = max([1 + floor(sample_sf*Rsurf), 2]);
    surf_r = (0:(npnts_r-1)) / sample_sf;
    % compute the surface wave heights
    if isfield(params, 'height'),
      % wave height means peak to trough distance
      h = 0.5* params.height;
    else
      error('surface_height: required field height is missing from params');
    end;
    omega = 2.0*pi*tfreq;
    k = 2.0*pi*sfreq;
    random_phase = 2.0*pi*rand(1);
    surf_h = zeros(npnts_r, npnts_t);
    for itime = 1:npnts_t,
      time = surf_t(itime);
      surf_h(:,itime) = h * cos(omega*time + random_phase) * ones(npnts_r,1);
    end;

  case {'sinusoidal', 'swell'}

    % period (seconds)
    if isfield(params, 'peak_period'),
      tau = params.peak_period;
    else
      error('surface_height: required field peak_period is missing from params');
    end;
    % compute temporal sample rate
    tfreq = 1.0 / tau;
    sample_tf = (2.0 * tfreq) * os_time;
    % compute temporal sample points
    npnts_t = 1 + ceil(sample_tf*Tsurf);
    surf_t = (0:(npnts_t-1)) / sample_tf;
    % wavelength (meters)
    lambda = g * tau^2 / (2.0*pi);
    % compute spatial sample rate
    sfreq = 1.0 / lambda;
    sample_sf = (2.0 * sfreq) * os_space;
    % compute spatial sample points
    npnts_r = 1 + ceil(sample_sf*Rsurf);
    surf_r = (0:(npnts_r-1)) / sample_sf;
    % compute the surface wave heights
    if isfield(params, 'height'),
      % wave height means peak to trough distance
      h = 0.5* params.height;
    else
      error('surface_height: required field height is missing from params');
    end;
    omega = 2.0*pi*tfreq;
    k = 2.0*pi*sfreq;
    random_phase = 2.0*pi*rand(1);
    surf_h = zeros(npnts_r, npnts_t);
    for itime = 1:npnts_t,
      time = surf_t(itime);
      surf_h(:,itime) = h * cos(omega*time - k*surf_r + random_phase);
    end;

  otherwise

    error(['surface_height: unknown surface model type: ', surface_type]);

end;	% switch surface_type

% that's all folks ...

return;

