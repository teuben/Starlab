#! /bin/sh
if [ "x$STARLAB_PATH" = "x/home/slevy/starlab" ]; then
  _do_echo=
else
  _do_echo=1
  STARLAB_PATH="/home/slevy/starlab"
  export STARLAB_PATH
fi

. $STARLAB_PATH/local/starlab_setup.sh

# for now, set STARLAB_VERSION, a few too many things in old .cshrc type files
# would break with this variable not set

STARLAB_VERSION=`cat $STARLAB_PATH/VERSION` 

#  always set a convenient shorter name for the last STARLAB_PATH loaded
STARLAB=$STARLAB_PATH
export STARLAB_VERSION STARLAB

if [ -n "$_do_echo" -a -n "$PS1" ]; then
  echo "Starlab version $STARLAB_VERSION loaded with STARLAB_PATH=$STARLAB_PATH"
fi

unset _do_echo
