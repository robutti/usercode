#!/bin/bash

alias cmsenv='eval `scramv1 runtime -sh`'
cd /afs/cern.ch/user/r/robutti/Workspace
pwd
source cmsSetup_620.sh
cd ERobutti/HoughTransChecks
cmsRun $PYCONFIGFILE
