
cases = [ 'halfspace      '; ...
          'ParaBot        '; ...
          'KoreanSeas     '; ...
          'Weymouth       '; ...
          'KermitRoosevelt'; ...
          'Munk           '; ...
          'MunkRot        '; ...
          'PenetrableWedge'; ...
          'Seamount       '; ...
          'Taiwan         '; ...
          'TruncatedWedge '; ...
          'PerfectWedge   '; ...
              ];


for icase = 1: size( cases, 1 )
    directory = deblank( cases( icase, : ) )
    eval( [ 'cd ' directory ] );
    delete *.prt
    delete *.shd
    delete *.grn
    delete *.mod
    delete *.shd.mat
    delete *.arr
    delete *.asv
    delete *.rts
    delete *.ray
    delete *.irc
    delete *.brc
    delete ._*
    cd ..
end

