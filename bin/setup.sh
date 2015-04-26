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

if ! [[ $HOSTNAME =~ t[23].*\.mit\.edu ]]
then
  # download the MitPhysics/data directory
  # At MIT T2/T3, data files are available at MIT CVMFS
  $CMSSW_BASE/src/MitPhysics/bin/updateData.sh
fi

# Generate ROOT dictionaries for classes defined in this module
$CMSSW_BASE/src/MitCommon/bin/genDict.sh MitPhysics/{FakeMods,Init,Mods,SelMods,Skim,Utils,Validation}

# check for existing qjets directory or install it  
$CMSSW_BASE/src/MitPhysics/bin/installQjets.sh

exit 0
