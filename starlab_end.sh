#   sourcing this script will remove all STARLAB filth from your bash/ksh environment

#				remove the usual bin/sbin from the path

if [ -n "$STARLAB_PATH" ]; then
    PATH=`echo :$PATH: | sed -e "s=:$STARLAB_PATH/s*bin:=:=g" -e 's=^:==' -e 's=:$=='`
fi

#   				remove any STARLAB_ variables 

_scrap=`env | sed -n -e '/^STARLAB_/s/=.*//p'`
if [ -n "$scrap" ]; then
  _nscrap=`echo $_scrap | wc -w`
  echo "Removing all $_nscrap STARLAB_* environment variables"
  unset $_scrap
fi
unset _scrap _nscrap

