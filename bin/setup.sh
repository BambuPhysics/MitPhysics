#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Setup the MitPhysics package adjusting things that are needed for it to compile and run properly.
#
#                                                                   Jan 16, 2014 - V0 Christoph Paus
#                                                                   Apr 26, 2015 - V1 Yutaro Iiyama
#---------------------------------------------------------------------------------------------------

echo "*************************"
echo " MitPhysics/bin/setup.sh"
echo "*************************"

if [[ $HOSTNAME =~ t[23].*\.mit\.edu ]]
then
  if ! [ -L $CMSSW_BASE/src/MitPhysics/data ]
  then
    ln -s /cvmfs/cvmfs.cmsaf.mit.edu/hidsk0001/cmsprod/cms/MitPhysics/data $CMSSW_BASE/src/MitPhysics/data
  fi
else
  # download the MitPhysics/data directory
  # At MIT T2/T3, data files are available at MIT CVMFS
  $CMSSW_BASE/src/MitPhysics/bin/updateData.sh
fi

# Generate ROOT dictionaries for classes defined in this module
$CMSSW_BASE/src/MitCommon/bin/genDict.sh MitPhysics/{FakeMods,Init,Mods,SelMods,Skim,Utils,Validation}

exit 0
