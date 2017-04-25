function Bellhop_inputs(params, surf_t, surf_r, surf_h)

% Bellhop_inputs - construct all of the needed BELLHOP input files
%
% $Id: Bellhop_inputs.m,v 1.4 2011/05/12 00:15:52 jcp Exp $

% Load the user's input file template for BELLHOP frame runs

if (isfield(params, 'template_env')),
  template_filename = params.template_env;
else
  % Default is 'Template.env'
  template_filename = 'Template.env';
end;

[ TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, fid ] = ...
            read_env(template_filename, 'BELLHOP');
fclose(fid);

% Extract the needed run parameters from the user specified params struct

if (isfield(params, 'bellhop_path')),
  bellhopexe = params.bellhop_path;
else
  % Default
  bellhopexe = 'bellhop.exe';
end;

if (isfield(params, 'calc_top_dir')),
  top_dir = params.calc_top_dir;
else
  % Default
  top_dir = 'BellhopRuns';
end;

if (isfield(params, 'calc_sub_dir')),
  sub_dir_fmt = params.calc_sub_dir;
else
  % Default
  sub_dir_fmt = 'frame_%04d';
end;

if (isfield(params, 'arrfil_fmt')),
  arrfil_fmt = params.arrfil_fmt;
else
  % Default
  arrfil_fmt = 'BELLHOP_%04d.arr';
end;

% Set implicitly valued params (in case of errors in the user's Template.env)

Opt = Bdry.Top.Opt;
Opt(5) = '*';
Bdry.Top.Opt = Opt;

% Copy any optional user specified params to the BELLHOP environment structures

if (isfield(params, 'freq')),
  freq = params.freq;
end;

if (isfield(params, 'src_depth')),
  Pos.s.depth = params.src_depth;
  Pos.Nsd = 1;
end;

if (isfield(params, 'rec_ranges')),
  rec_ranges = params.rec_ranges;
  Pos.Nrr = length(rec_ranges);
  Pos.r.range = rec_ranges;
end;

if (isfield(params, 'rec_depths')),
  rec_depths = params.rec_depths;
  Pos.Nrd = length(rec_depths);
  Pos.r.depth = rec_depths;
end;

if (isfield(params, 'nbeams')),
  Beam.Nbeams = params.nbeams;
end;

if (isfield(params, 'angles')),
  Beam.alpha = params.angles;
end;

if (isfield(params, 'runtype')),
  Beam.RunType = params.runtype;
end;

if (isfield(params, 'btyfil')),
  btyfil = params.btyfil;
else
  btyfil = [];
end;

% Save the current working directory (where VIRTEX is run from)

virtex_cwd = pwd;

% Change to top of the directory tree where the BELLHOP calculations are done

if (~exist(top_dir, 'dir')),
  % huh? directory does not exist, create it...
  mkdir(top_dir);
end;

cd(top_dir);

% Create the master BELLHOP environment file

write_env('BELLHOP.env', 'BELLHOP', TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax);

% Loop over all the time steps

ntime = length(surf_t);

for itime = 1:ntime,

  sub_dir_name = sprintf(sub_dir_fmt, itime);

  % Check for the existance of the working directory for this frame

  if (~exist(sub_dir_name, 'dir')),
    % huh? sub-directory does not exist, create it...
    mkdir(sub_dir_name);
  end;

  % Change to the working directory for this frame

  cd(sub_dir_name);

  % Write the BELLHOP altimetry for this frame

  system('rm -f BELLHOP.ati');

  fid = fopen('BELLHOP.ati','w');
  fprintf(fid,'%s\n', '''C''');
  fprintf(fid,'%i\n', length(surf_r));
  for irange = 1:length(surf_r),
    % Note that altimetry files specify the depth, (or -1 * surface height)
    fprintf(fid,'%.8f %.8f\n', 0.001*surf_r(irange), -surf_h(irange,itime));
  end;
  fclose(fid);

  % Copy the bathymetry template file (if one was specified)

  if (length(btyfil) > 0),
    cp_btyfil_cmd = sprintf('cp -p %s/%s ./BELLHOP.bty', virtex_cwd, btyfil);
    system(cp_btyfil_cmd);
  end;

  % Write the (Bourne) shell script to execute BELLHOP

  arrfil_name = sprintf(arrfil_fmt, itime);

  fid = fopen('bellhop.sh', 'w');

  fprintf(fid, '#!/bin/sh\n' );
  fprintf(fid, '#\n' );
  fprintf(fid, '# first, be certain we are in the correct dir\n');
  fprintf(fid, 'cd %s\n', pwd());
  fprintf(fid, '# remove the BELLHOP status file for a clean start\n');
  fprintf(fid, 'rm -f bellhop.rc\n');
  fprintf(fid, '# copy the "master" env file (up one dir) to current dir\n');
  fprintf(fid, 'cp -p ../BELLHOP.env ./\n');
  fprintf(fid, '# execute BELLHOP using the copy of the "master" env file\n');
  fprintf(fid, '%s BELLHOP > bellhop.log 2>&1\n', bellhopexe);
  fprintf(fid, '# save the exit status to the BELLHOP status file\n');
  fprintf(fid, 'echo $? > bellhop.rc\n');
  fprintf(fid, '# verify the arrivals file exists before trying to move it\n');
  fprintf(fid, 'if [ -s BELLHOP.arr ]; then\n');
  fprintf(fid, '  mv BELLHOP.arr %s/%s\n', virtex_cwd, arrfil_name);
  fprintf(fid, 'fi\n');

  fclose(fid);

  % Remove any BELLHOP status files left over from previous runs

  system('rm -f bellhop.rc');

  % Change back to the top of the directory tree

  cd('..');

  % Next frame

end;

% Change back to the current working directory (where VIRTEX is run from)

cd(virtex_cwd);

% That's all folks...

return;

