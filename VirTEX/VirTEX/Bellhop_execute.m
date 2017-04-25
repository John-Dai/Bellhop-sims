function Bellhop_execute(params, ntime, jobLoLevel, jobHiLevel)

% Bellhop_execute - execute all of the BELLHOP input files
%
% $Id: Bellhop_execute.m,v 1.2 2010/05/26 20:11:43 jcp Exp $

% Extract the needed run parameters from the user specified params struct

pause_sec = 0.10;

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

% Save the current working directory (where VIRTEX is run from)

virtex_cwd = pwd;

% Change to top of the directory tree where the BELLHOP calculations are done

cd(top_dir);

% Loop over all the time steps, generate the frame sub-directory names

frame_dirs = cell(ntime, 1);

for itime = 1:ntime,
  sub_dir_name = sprintf(sub_dir_fmt, itime);
  frame_dirs(itime) = cellstr(sub_dir_name);
end;

% Initialize the job status array

jobStatus = -1 * ones(1, ntime);	% flags to indicate frame job status
					% -1 Job has not been started yet
					%  0 Job was started, status unknown
					%  1 Job completed successfully
					% >1 Job completed with errors

% Start the first jobHiLevel jobs (or ntime jobs, if that is smaller)

for itime = 1:min([ntime, jobHiLevel]),
  % get the name of the sub-directory for this frame
  frame_dir = char(frame_dirs(itime));
  % start this job (in the background)
  system(['/usr/bin/env -i sh ', frame_dir, '/bellhop.sh &']);
  % tag this job as started
  jobStatus(itime) = 0;
  % feedback that job was started
  fprintf('.');
end;

% Loop until all frame jobs have been launched and completed

while any(jobStatus < 1),

   % get the indices of calculations that were started, but are not finished

   jobList = find(jobStatus == 0);

   % check each frame to see if the associated BELLHOP job has completed

   njobs = length(jobList);

   for itime = jobList,

     % get the name of the sub-directory for this frame

     frame_dir = char(frame_dirs(itime));

     % check for the existance of the BELLHOP status file first

     bellhop_status_file = [ frame_dir '/bellhop.rc' ];

     if ~exist(bellhop_status_file, 'file'),
        % nothing to see here, move along...
        continue;
     end;

     % check exit status of the BELLHOP calculation

     fid = fopen(bellhop_status_file, 'r');
     rc = fscanf(fid, '%d');
     fclose(fid);

     fprintf('(%d)', rc);

     if (rc ~= 0),
        % DOH !!!
        fprintf('Bellhop_execute: run for frame %03d had errors!', itime);
        jobStatus(itime) = 2;
        njobs = njobs - 1;
        continue;
     else
        % Ehhhxcellent ...
        jobStatus(itime) = 1;
        njobs = njobs - 1;
     end;

   end;   % next itime

   % decide if we should start more jobs

   jobList = find(jobStatus == 0);

   njobs = length(jobList);

   if (njobs < jobLoLevel),

     % get the indices of BELLHOP runs that have NOT been started yet

     jobList = find(jobStatus < 0);

     nleft = length(jobList);

     if (nleft > 0),

       % number of jobs to start

       nstart = min([ nleft, jobHiLevel-njobs ]);

       for ijob = 1:nstart,

         itime = jobList(ijob);

         % get the name of the sub-directory for this frame

         frame_dir = char(frame_dirs(itime));

         % start this job (in the background)

         system(['/usr/bin/env -i sh ', frame_dir, '/bellhop.sh &']);

         jobStatus(itime) = 0;	% tag this job as started

         fprintf('.');	% feedback that job was started

       end;	% if

     end;	% if nleft > 0

   else

     % wait a moment before checking again

     pause(pause_sec);

   end;		% if njobs < jobLoLevel

end;   % while any(jobStatus < 1)

fprintf('\n');

% Change back to the current working directory (where VIRTEX is run from)

cd(virtex_cwd);

% That's all folks...

return;

