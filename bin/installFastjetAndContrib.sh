#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Install fastjet and all relevant contributions to be able to fully compile MitPhysics.
#
# Use this script when the default CMSSW version of fastjet is outdated.
#
#                                                                   Jan 13, 2014 - V0 Christoph Paus
#                                                                   Apr 26, 2015 - V1 Yutaro Iiyama
#---------------------------------------------------------------------------------------------------
function configureScram {
  local EXTERNAL=$1
  local VERSION=$2

  if ! [ $VERSION ]
  then
    echo " Could not determine fastjet version. Please check installation at"
    echo " $EXTERNAL"
    exit 1
  fi

  # add local fastjet external to scarm config
  mv $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml \
     $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml-last.$(date +%s)

  echo \
'  <tool name="fastjet" version="'$VERSION'">
    <info url="http://www.lpthe.jussieu.fr/~salam/fastjet/"/>
    <lib name="fastjetplugins"/>
    <lib name="fastjettools"/>
    <lib name="siscone"/>
    <lib name="siscone_spherical"/>
    <lib name="fastjet"/>
    <lib name="fastjetcontrib"/>
    <client>
      <environment name="FASTJET_BASE" default="'$EXTERNAL'"/>
      <environment name="LIBDIR" default="$FASTJET_BASE/lib"/>
      <environment name="INCLUDE" default="$FASTJET_BASE/include"/>
    </client>
  </tool>
' > $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml

  # commit the scram config changes
  cd $CMSSW_BASE/src
  scram setup fastjet
}

####################################################
### EDIT THE LINES BELOW TO SPECIFY THE VERSIONS ###
####################################################

FASTJET="fastjet-3.0.6"
FJCONTRIB="fjcontrib-1.011_nsub-2.0.0-rc3"

####################################################

# test whether environment is clean and we can start
if [ -z "$CMSSW_BASE" ]
then
  echo ""
  echo " ERROR - cmssw is not setup. Please, setup release directory and checkout MitPhysics."
  echo ""
  exit 1
fi
if [ -d "$CMSSW_BASE/src/MitPhysics" ]
then
  echo ""
  echo " INFO - found MitPhysics location at: $CMSSW_BASE/src/MitPhysics"
  echo ""
else
  echo ""
  echo " ERROR - MitPhysics is not in your release. Please check it from GITHUB/CVS."
  echo ""
  exit 1
fi

FASTJET_VERSION=`echo $FASTJET | cut -d '-' -f2`

# the default location
EXTERNAL=/cvmfs/cvmfs.cmsaf.mit.edu/hidsk0001/cmsprod/cms/external
if [ -d $EXTERNAL/$FASTJET ]
then
  configureScram $EXTERNAL $FASTJET_VERSION
  exit 0
fi

EXTERNAL=/home/cmsprod/cms/external
if [ -d $EXTERNAL/$FASTJET ]
then
  configureScram $EXTERNAL $FASTJET_VERSION
  exit 0
fi

EXTERNAL="/home/$USER/cms/external"
echo " INFO - make own repository at: $EXTERNAL/$FASTJET"
echo ""

# Here the real work starts (fastjet, fjcontrib, tweaks)

mkdir -p $EXTERNAL

# install best fastjet version
#-----------------------------

# set all relevant variables
FASTJET_URL="http://fastjet.fr/repo"

# in the right location
cd $EXTERNAL

# now do the download
echo " INFO - download starting"
echo ""

wget "$FASTJET_URL/$FASTJET.tar.gz" -O $FASTJET.tar.gz

# unpack it (no debug, could add 't' option)
echo " INFO - unpacking"
echo ""

tar fzx $FASTJET.tar.gz

# cleanup to avoid junk
echo " INFO - cleaning up"
echo ""

rm -rf $FASTJET.tar.gz

# installing fastjet
echo " INFO - installing"
echo ""

cd $EXTERNAL/$FASTJET

./configure --prefix=$EXTERNAL
make
make check
make install


# install best fastjet contribution version
#------------------------------------------

# define all relevant variables
FJCONTRIB_URL="http://fastjet.hepforge.org/contrib/downloads"

# in the right location
cd $EXTERNAL

# now do the download
echo " INFO - download starting"
echo ""

wget "$FJCONTRIB_URL/$FJCONTRIB.tar.gz" -O $FJCONTRIB.tar.gz

# unpack it (no debug, could add 't' option)
echo " INFO - unpacking"
echo ""

tar fzx $FJCONTRIB.tar.gz

# cleanup to avoid junk
echo " INFO - cleaning up"
echo ""

rm -rf $FJCONTRIB.tar.gz

# installing fastjet
echo " INFO - installing"
echo ""

cd $EXTERNAL/$FJCONTRIB

./configure --fastjet-config=$EXTERNAL/$FASTJET/fastjet-config --prefix=$EXTERNAL CXXFLAGS="-O3 -Wall -Woverloaded-virtual -g -fPIC -I$EXTERNAL/include"
make
make check
make install

# make shared libraries for the contributions
g++ -shared -fPIC -o $EXTERNAL/lib/libfastjetcontrib.so -Wl,-soname,libfastjetcontrib.so $EXTERNAL/$FJCONTRIB/*/[A-Z]*.o 

# final adjustment to scram configuration
configureScram $EXTERNAL $FASTJET_VERSION

exit 0
