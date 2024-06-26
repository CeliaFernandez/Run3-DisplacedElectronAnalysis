# Run 3 superclusters analysis

## How to install

Recommended release for this analyzer is CMSSW_12_6_0 or later. Commands to setup the analyzer are:

```
cmsrel CMSSW_12_6_0

cd CMSSW_12_6_0/src

cmsenv

mkdir Analysis

cd Analysis

git clone git@github.com:CeliaFernandez/Run3-SuperclusterAnalysis.git

scram b -j 8
```

## Ntuple preparation
[...]

## Run the notebook

Within ```notebooks/```:

(0) Setup the environment for the first time (do it only one time):
```
source setup.sh
```

(1) Then in order to activate the environment:
```
source env/bin/activate
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc10/cms/cmssw/CMSSW_12_6_0/src ; eval `scramv1 runtime -sh` ; cd -
```

(2) To run the notebook:
```
jupyter notebook --no-browser --port=8893
```

(3) Outside, in your local computer, push the port to the remote machine:
```
ssh -N -f -L localhost:8893:localhost:8893 username@uaf-10.t2.ucsd.edu
```

(4) Work from the browser


