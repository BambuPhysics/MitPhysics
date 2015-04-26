#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Setup the MitPhysics package adjusting things that are needed for it to compile and run properly.
#
#                                                                   Jan 16, 2014 - V0 Christoph Paus
#                                                                   Apr 26, 2015 - V1 Yutaro Iiyama
#---------------------------------------------------------------------------------------------------

# check for existing qjets directory or install it  
$CMSSW_BASE/src/MitPhysics/bin/installQjets.sh

exit 0
