#!/bin/bash
#---------------------------------------------------------------------------------------------------------------------
#
# When CMSSW is setup in the condor area the external libraries are also copied but need to be setup with scram so
# they appear in the right places for linking. We can assume that when this script this called CMSSW is already setup.
#
#                                                                                            C.Paus, V0 (Mar 05, 2014)
#                                                                                          Y.Iiyama, V1 (Jul 13, 2015)
#---------------------------------------------------------------------------------------------------------------------

#---------------------
# SET PARAMETERS HERE
#
# PACKAGES TO INSTALL
PACKAGES="pwhg_cphto_reweight"
#---------------------

if ! [ "$CMSSW_BASE" ]
then
  echo "CMSSW_BASE not set"
  exit 1
fi

# Utility function to find the installation of the external
find-external() {
  local PACKAGE=$1
  local EXTERNAL

  # look in CMSSW default first
  for EXTERNAL in \
  /cvmfs/cms.cern.ch/$SCRAM_ARCH/external/$PACKAGE \
  /cvmfs/cvmfs.cmsaf.mit.edu/hidsk0001/cmsprod/cms/external/$PACKAGE
  do
    if [ -d $EXTERNAL ]
    then
      echo $EXTERNAL
      exit 0
    fi
  done

  EXTERNAL="/home/$USER/cms/external/$PACKAGE"
  echo " INFO - making directory at: $EXTERNAL"
  echo ""

  mkdir -p $EXTERNAL
  echo $EXTERNAL
}

### Installer functions

# Install fastjet
install-pwhg_cphto_reweight() {
  # add local fastjet external to scarm config

  local BASE=$(find-external pwhg_cphto_reweight)
  local TARGET=$CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/pwhg_cphto_reweight.xml

  echo ""
  echo " INFO - using pwhg_cphto_reweight at: $BASE"
  echo ""

  mv $TARGET ${TARGET}-last.$(date +%s) 2> /dev/null

  echo \
'  <tool name="pwhg_cphto_reweight" version="1">
    <lib name="cphto"/>
    <client>
      <environment name="LIBDIR" default="'$BASE'"/>
    </client>
  </tool>
' > $TARGET

  # commit the scram config changes
  cd $CMSSW_BASE/src
  scram setup pwhg_cphto_reweight
}

### Loop over PACKAGES

for PACKAGE in $PACKAGES
do
  case $PACKAGE in
    pwhg_cphto_reweight)
      install-pwhg_cphto_reweight
      ;;
    *)
      echo "Unknown external $PACKAGE"
      exit 1
      ;;
  esac
done
