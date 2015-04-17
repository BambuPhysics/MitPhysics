#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Setup the MitPhysics package adjusting things that are needed for it to compile and run properly.
#
#                                                                       Apr 17, 2015 - V1 Y. Iiyama
#---------------------------------------------------------------------------------------------------

if ! [[ $HOSTNAME =~ t[23].*\.mit\.edu ]]
then
  # download the MitPhysics/data directory
  # At MIT T2/T3, data files are available at MIT CVMFS
  $CMSSW_BASE/src/MitPhysics/bin/updateData.sh
fi

# Generate ROOT dictionaries for classes defined in this module
$CMSSW_BASE/src/MitCommon/bin/genDict.sh MitPhysics/{FakeMods,Init,Mods,SelMods,Skim,Utils,Validation}

# check for existing fastjet+contribution directory or install it
$CMSSW_BASE/src/MitPhysics/bin/installFastjetAndContrib.sh

# check for existing qjets directory or install it  
$CMSSW_BASE/src/MitPhysics/bin/installQjets.sh

exit 0
