#!/bin/sh
#
# ksh/bash users source this file to install and build Starlab from scratch:
#
# Usage:	.  ./AAA_SOURCE_ME  [options]
#
# Note: All command-line arguments are passed to 'make config_new'

#
# It is not foolproof, but still intended as such. If it does not work, keep 
# the install.log and study it (or email the authors). Most of the time you
# will need to edit the local/cshrc.starlab file, source it, and try the
# 'make all' again
#

stat=0

trap 'if [ $stat -ne 0 ]; then
       echo "### $stat errors occured during execution of the AAA_SOURCE_ME.sh script"
      fi
      exit $stat
    ' 0 1 2 15

#	make sure the script runs from the current directory
if [ ! -f AAA_SOURCE_ME.sh -o ! -f starlab_start.in ]; then
   echo "Must . AAA_SOURCE_ME.sh  from the top-level starlab directory!"
   stat=1
   exit 1
fi

rm -f install.log 2> /dev/null
touch install.log


echo "(Sending output to install.log)"

if [ ! -f starlab_start ]; then
    echo Running configure
    ./configure >> install.log  2>&1
    if [ $? -ne 0 ]; then
	stat=`expr $stat + 1`
    fi
    echo '++++++++++++++++++++++++++++++++++++++++' >> install.log
    echo Installing configuration files
    make config_new "$@" >> install.log 2>&1
    if [ $? -ne 0 ]; then
	stat=`expr $stat + 1`
    fi
    echo '++++++++++++++++++++++++++++++++++++++++' >> install.log
fi

echo Initializing Starlab environment
. ./starlab_start.sh >> install.log 2>&1
if [ $? -ne 0 ]; then
    stat=`expr $stat + 1`
fi
echo '++++++++++++++++++++++++++++++++++++++++' >> install.log

echo "Building Starlab (this may take a few minutes)..."
make all >> install.log 2>&1
if [ $? -ne 0 ]; then
    stat=`expr $stat + 1`
fi
echo '++++++++++++++++++++++++++++++++++++++++' >> install.log

if [ $stat -eq 0 ]; then
    echo -n Testing kira...
    kira --help 2> /dev/null
    if [ $? -eq 0 ]; then
        echo 'OK :-)'
    else
        echo 'failed :-('
    fi

    echo "Done...   Run 'scripts/kira_test' to check a sample kira run output."
    echo ""
    echo "Note: Although your current shell now has the Starlab environment loaded,"
    echo "      new shells will not have Starlab pre-loaded. You would need to add"
    echo "      the following command/alias to your .cshrc (or equivalent) file:"
    echo "                   . $STARLAB_PATH/starlab_start.sh"

else
    echo Install failed: stat = $stat
fi
