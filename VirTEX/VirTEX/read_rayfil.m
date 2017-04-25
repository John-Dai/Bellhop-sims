function [Hdr, Rays] = read_rayfil(RAYFIL)

% read_rayfil - read the ASCII RAYFIL written by BELLHOP
%
% $Id: read_rayfil.m,v 1.2 2011/05/22 02:56:00 jcp Exp $

% attempt to open the RAYFIL

fid = fopen(RAYFIL, 'r');

if (fid == -1)
  error( [ mfilename, ': error opening the RAYFIL: ', RAYFIL ] );
end;

% read RAYFIL header

title  = fgetl(  fid );
freq   = fscanf( fid, '%f', 1 );
NbeamA = fscanf( fid, '%i', 2 );
depthT = fscanf( fid, '%f', 1 );
depthB = fscanf( fid, '%f', 1 );

% remove single quotes and extraneous blanks from title string

nchars = strfind(title, '''');
title = title(nchars(1)+1:nchars(2)-1);

Nbeams = NbeamA(1);
Hdr.title  = title;
Hdr.freq   = freq;
Hdr.Nbeams = Nbeams;
Hdr.depthT = depthT;
Hdr.depthB = depthB;

% read rays

for ibeam = 1:Nbeams,

  Alpha0    = fscanf(fid, '%f', 1);
  Nsteps    = fscanf(fid, '%i', 1);
  NumTopBnc = fscanf(fid, '%i', 1);
  NumBotBnc = fscanf(fid, '%i', 1);

  if isempty(Nsteps),
    break;
  end;

  Rays(ibeam).Alpha0    = Alpha0;
  Rays(ibeam).Nsteps    = Nsteps;
  Rays(ibeam).NumTopBnc = NumTopBnc;
  Rays(ibeam).NumBotBnc = NumBotBnc;

  points = fscanf(fid, '%f', [2 Nsteps]);

  Rays(ibeam).r = points(1,:);
  Rays(ibeam).d = points(2,:);
   
end;	% next beam

fclose( fid );

