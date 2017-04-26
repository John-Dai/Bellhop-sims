function [ Modes ] = read_modes( filename, modes )

% Read the modes produced by KRAKEN
% usage:
%    [ Modes ] = read_modes( filename, modes )
% filename should include the extension
% iProf is the profile number
% modes is an optional vector of mode indices

% mbp, May 2001

% identify the file type

% !!! put the following line in your calling routine, for the first call
% clear read_modes_bin % to force rewind to beginning of mode file

switch filename(1 : 6)
    case 'MODFIL'
        FileType = 'mod';
    case 'MOAFIL'
        FileType = 'moa';
    otherwise
        endchar = length( filename );
        if ( endchar >= 4 )
            FileType = lower( filename( endchar-2 : endchar ) );
        end
end

% read the modal data

switch FileType
    case 'mod' % binary format
        if nargin == 1
            Modes = read_modes_bin( filename );
        else
            Modes = read_modes_bin( filename, modes );
        end
    case 'mat' % Matlab format
        load( filename );
    case 'moa' % ascii format
        if nargin == 1
            Modes = read_modes_asc( filename );
        else
            Modes = read_modes_asc( filename, modes );
        end
        
    otherwise
        warndlg( 'Unrecognized file extension', 'Warning' )
end

% calculate wavenumbers in halfspaces (if there are any modes)

if ( Modes.M ~= 0 )
    
    if ( Modes.Top.bc == 'A' )   % top
        Modes.Top.k2     = ( 2 * pi * Modes.freq / Modes.Top.cp )^2;
        gamma2           = Modes.k .^ 2 - Modes.Top.k2;
        Modes.Top.gamma  = PekerisRoot( gamma2 );   % vertical wavenumber
        Modes.Top.phi    = Modes.phi( 1, : );       % mode value at halfspace
    else
        Modes.Top.rho   = 1.0;
        Modes.Top.gamma = zeros( size( Modes.k ) );
        Modes.Top.phi   = zeros( size( Modes.phi( 1, : ) ) );
    end
    
    if ( Modes.Bot.bc == 'A' )   % bottom
        Modes.Bot.k2    = ( 2 * pi * Modes.freq / Modes.Bot.cp )^2;
        gamma2          = Modes.k .^ 2 - Modes.Bot.k2;
        Modes.Bot.gamma = PekerisRoot( gamma2 );    % vertical wavenumber
        Modes.Bot.phi   = Modes.phi( end, : );      % mode value at halfspace
    else
        Modes.Bot.rho   = 1.0;
        Modes.Bot.gamma = zeros( size( Modes.k ) );
        Modes.Bot.phi   = zeros( size( Modes.phi( end, : ) ) );
    end
end