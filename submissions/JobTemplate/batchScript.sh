#!/bin/bash
set -euo pipefail

haha1=$3
haha2=$4

if [ -z ${_CONDOR_SCRATCH_DIR+x} ]; then
  #running locally
  runninglocally=true
  _CONDOR_SCRATCH_DIR=$(mktemp -d)
  SUBMIT_DIR=$(pwd)
else
  runninglocally=false
  SUBMIT_DIR=$1
fi

cd $SUBMIT_DIR
dir_name=$(basename $(pwd))

eval $(scram ru -sh)


cp ../macros/MACRO.exe $_CONDOR_SCRATCH_DIR
mkdir $_CONDOR_SCRATCH_DIR/jsons
cp -r ../jsons/* $_CONDOR_SCRATCH_DIR/jsons
cd $_CONDOR_SCRATCH_DIR

export HOME=$_CONDOR_SCRATCH_DIR

echo 'Running at:' $(date)
echo path: `pwd`

cmsRunStatus=   #default for successful completion is an empty file
./MACRO.exe "FILELOC" "IDENTIFIER" "XSEC" "NUMEVENTS" $haha1 $haha2 > log.txt || cmsRunStatus=$?

echo -n $cmsRunStatus > exitStatus.txt
echo 'cmsRun done at: ' $(date) with exit status: ${cmsRunStatus+0}
#gzip log.txt

export ROOT_HIST=0
#if [ -s testoutput.root ]; then
# root -q -b '${CMSSW_BASE}/src/Run3VBSPol/Run3VBSPol/test/rootFileIntegrity.r("testoutput.root")'
#else
# echo moving empty file
# mv testoutput.root testoutput.root.empty
#fi

echo "Files on node:"
ls -la

#delete mela stuff and $USER.cc
rm -f br.sm1 br.sm2 ffwarn.dat input.DAT process.DAT "$USER.cc"

#delete submission scripts, so that they are not copied back (which fails sometimes)
#rm -rf MACRO.exe batchScript.sh jsons/

echo '...done at' $(date)

#note cping back is handled automatically by condor
if $runninglocally; then
  cp testoutput.root* *.txt *.gz $SUBMIT_DIR
fi

echo "Copy to /eos/user/l/lallasia/JobOutput/"$dir_name
mv testoutput.root /eos/user/l/lallasia/JobOutput/$dir_name

exit $cmsRunStatus
