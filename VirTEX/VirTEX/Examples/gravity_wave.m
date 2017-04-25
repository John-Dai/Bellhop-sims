function h = gravity_wave(t, r, sparams)

% GRAVITY_WAVE - compute the height of a simple swell (gravity) surface wave
%
% Usage: h = gravity_wave(t, r, sparams)
%
%  Inputs:
%    t       - time values (real values array of length N)
%    r       - range value (real values array of length 1)
%    sparams - structure for setting the wave parameters (first call)
%
%  Output:
%    h       - surface height values (real values array of length N)
%
% $Id: gravity_wave.m,v 1.1 2011/03/16 23:13:00 jcp Exp $

persistent h0 k omega phi_0;

% check for required arguments

if nargin < 2,
  error( [mfilename, ': one or more required input arguments is missing'] );
end;

% check for optional arguments

if nargin > 2,

  % store wave parameters in local (persistent) variables for future use
  h0    = sparams.h0;
  k     = sparams.k;
  omega = sparams.omega;
  phi_0 = sparams.phi_0;

else

  % compute the surface height for given time, range values

  h = h0 * sin(omega*t - k*r + phi_0);

end;

return;

% end of gravity_wave.m
