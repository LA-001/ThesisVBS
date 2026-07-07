#!/bin/tcsh
# Search for missing root files
# parameter to be specified: mf = move failed jobs
#                            lf = link failed jobs (then you can run cleanup.csh; resubmit_Condor.csh in AAAFAIL)

#set echo

set opt=$1

@ runcounter = 0
@ failcounter = 0

if ( $#argv >= 2 ) then
  set gooddir=$2
else
  set gooddir=AAAOK
endif

if ( $#argv >= 3 ) then
  set baddir=$3
else
  set baddir=AAAFAIL
endif

foreach chunk ( *Chunk* )
 set fail="false"


 # Check that root file is existing and not empty
 set filename=/eos/user/l/lallasia/JobOutput/${chunk}/testoutput.root
 if ( ! -e $filename ) then
#   echo "Missing root file in " ${chunk}
   set fail="true"
 else if ( -z $filename ) then
#   echo "Empty file: " $filename
   set fail="true"
 endif

 # Check job exit status. Cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/JobExitCodes , https://twiki.cern.ch/twiki/bin/view/CMSPublic/StandardExitCodes
 set exitStatus = 0
 if(-e ${chunk}/exitStatus.txt) then
   set filesize = `ls -l ${chunk}/exitStatus.txt | awk '{print $5}'`
   if ( $filesize > 1 ) then
     set exitStatus=`cat ${chunk}/exitStatus.txt`
     set fail="true"
   endif
 else if(! -e ${chunk}/exitStatus.txt) then
   set fail="true"
 endif

 # Check for failures reported in the Condor log
 if ( $exitStatus == 0 || $exitStatus == 1) then 
  set nonomatch
  set logFile = ( ${chunk}/log/*.log )
  if ( -e $logFile[1] ) then
    if ( `grep -c -e "Job removed.*due to wall time exceeded" $logFile[$#logFile]` != 0 ) then
      set exitStatus=-152
      set fail="true"
    else if ( `grep -c -e "The job attribute PeriodicRemove expression.*evaluated to TRUE" $logFile[$#logFile]` != 0 ) then
      set exitStatus=-153
      set fail="true"
    else if ( `grep -c -e "Job was aborted by the user" $logFile[$#logFile]` != 0 ) then
      set exitStatus=-154
      set fail="true"
    endif
  endif
  unset nonomatch
 endif

 # Archive succesful jobs, or report failure
 if ( $fail == "false" ) then
  mkdir -p $gooddir
  mv $chunk $gooddir/
  echo $chunk ": DONE"
 else
  set description=""
   if ( $exitStatus == 0 || $exitStatus == 1) then
     # is the job terminated?
     set nonomatch
     set logFile = ( ${chunk}/log/*.log )
     if ( -e $logFile[1] ) then
	if ( `grep -c -e "Job terminated" $logFile[$#logFile]` != 0 ) then
	    echo $chunk ": terminated, unknown failure"
	    @ failcounter = $failcounter + 1
	else
	    @ runcounter = $runcounter + 1
	    echo $chunk ": still running (or unknown failure)"
	    endif
     else
	echo $chunk ": still running (or unknown failure)"
	@ runcounter = $runcounter + 1
     endif
     unset nonomatch
   else
     @ failcounter = $failcounter + 1
     if ( $exitStatus == 84 ) set description="(missing input file)"
     if ( $exitStatus == 85 ) set description="(error reading file)"
     if ( $exitStatus == 92 ) set description="(failed to open file)"
     if ( $exitStatus == 127 ) set description="(probably an error transfering files)"
     if ( $exitStatus == 134 ) set description="(Crashed)"
     if ( $exitStatus == 137 ) set description="(killed - probably Condor problem)"
     if ( $exitStatus == 139 ) set description="(error reading file)"
     if ( $exitStatus == 152 || $exitStatus == -152) set description="(Exceeded CPU time)"
     if ( $exitStatus == 153 || $exitStatus == -153) set description="(Condor crashed, see log file)"
     if ( $exitStatus == 154 || $exitStatus == -154) set description="(You cancelled the job)"
    echo $chunk ": failed, exit status = " $exitStatus $description
   endif
   if ( $opt == "mf" && $exitStatus != 0 ) then
     mkdir -p $baddir
     mv $chunk $baddir/
   else if ( $opt == "lf" && $exitStatus != 0 ) then
     mkdir -p $baddir
     ln -s ../$chunk $baddir/
     ln -sf ../condor.sub $baddir
   endif
 endif
end

set donecounter=`ls -1 AAAOK | wc -l`

set bridgetext = ""
set runtext = ""
if $runcounter != 0 then
    set runtext = "${bridgetext}\033[1;33m${runcounter} jobs running\033[0m"
    set bridgetext = ", "
    endif

set failtext = ""
if $failcounter != 0 then
    set failtext = "${bridgetext}\033[1;31m${failcounter} jobs failed\033[0m"
    set bridgetext = ", "
    endif

set donetext = ""
if $donecounter != 0 then
    set donetext = "${bridgetext}\033[1;32m${donecounter} jobs DONE\033[0m"
    endif

echo "\033[1;4mSummary\033[0m\033[1m: \033[0m"$runtext$failtext$donetext
