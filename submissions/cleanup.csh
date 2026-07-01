#!/bin/tcsh

set nonomatch
set list = ( */*.root */*.corrupted */*.recovered *_Chunk/*.gz */*.txt */core* */jobid  */LSFJOB*/ */log/* */output/* */error/* */*.DAT */*.cc *_Chunk/*.exe)

foreach f ( ${list} )
    if ( -e $f ) then
    rm -r $f
    endif
end
