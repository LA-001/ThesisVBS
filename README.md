Repository for VBS analysis with taus in the final state.


# Setup

```
cmsrel CMSSW_14_0_9
cd CMSSW_14_0_9/src
git clone https://github.com/gmarozzo/Run3VBSWithTaus
cmsenv
scram b	-j 8
```

You will need to replace my local path with yours, in multiple locations. The ones I can think of (there might be others): submissions/checkProd.csh submissions/CreateJobs.py submissions/hadd.py 


# Basic usage

To compile the analysis macro (you can save this command as an alias in .bashrc if you want):

```
cd submissions/macros
g++ $(correction config --cflags --ldflags --rpath) -o "VBS.exe" "VBS.C" $(root-config --cflags --glibs)
cd ..
```

To create and submit jobs (from inside the "submissions" folder):

```

voms-proxy-info -e || voms-proxy-init --voms cms
python3 CreateJobs.py macros/VBS samples/samplelist_2024.csv
csh resubmit_Condor.csh
```

You can check the status of the jobs with condor_q. After they are done, check the status:

```
csh checkProd.csh
```

If some failed:

```
csh cleanup.csh && csh resubmit_Condor.csh
```

When all jobs are done:

```
python3 hadd.py
```