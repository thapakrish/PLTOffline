#!/bin/tcsh
setenv ROOTSYS /nfs/cern/root_v5.34.28_x64
setenv PATH $ROOTSYS/bin:${PATH}
setenv LD_LIBRARY_PATH {$ROOTSYS}/lib:${LD_LIBRARY_PATH}
echo "source /nfs/cern/root_v5.34.28_x64/bin/thisroot.csh"
#export ROOTSYS=/data/CernRoot/root_v5.34.28/root
#export PATH=$PATH:/data/CernRoot/root_v5.34.28/root/bin/
#export LD_LIBRARY_PATH=/data/CernRoot/root_v5.34.28/root/lib/

#export ROOTSYS=/nfs/cern/root_v5.26_x64
#export PATH=$PATH:/nfs/cern/root5.26_x64/root/bin/
#export LD_LIBRARY_PATH=/nfs/cern/root5.26_x64/root/lib/

