#!/bin/csh -f

# Starlab installation script for "developers" (who get Starlab from CVS).
# The extra layers are necessary to ensure that the correct autotools
# are always available when needed.

# Require that STARLAB_PATH be defined.

if ($?STARLAB_PATH == 0) then
    setenv STARLAB_PATH $PWD
endif
echo STARLAB_PATH = $STARLAB_PATH

set install_dir = $STARLAB_PATH/usr
echo install_dir \ = $install_dir

set log = install.log
echo Writing detailed output to $log

# Install autotools (place in install_dir):

echo Installing local versions of autotools in install_dir
(chdir src/packages; ./install-autotools $install_dir)	>&! $log

# Make sure that install_dir precedes system directories in search path

set path = ($install_dir/bin $path)
rehash

# First time through we need to do autoreconf

if (! -e configure) then
    echo Running autoreconf to build configure
    echo Running autoreconf to build configure		>&! $log
    autoreconf						>>& $log
endif

set OK = 0

if (-e configure) then

    # Perform the usual configuration process:

    echo Running configure
    echo 						>>& $log
    echo Running configure				>>& $log
    ./configure --prefix=$install_dir			>>& $log

    # Build the package.

    echo Running make
    echo 						>>& $log
    echo Running make					>>& $log
    make						>>& $log

    echo Running make install
    echo 						>>& $log
    echo Running make install				>>& $log
    make install					>>& $log

    # Then test whatever...

    rehash
    echo Testing...

    echo Looking for kira
    set which_kira = (`which kira | grep -v 'Command not found.'`)
    if ($#which_kira > 0) then
	echo -n kira':  '; kira --help |& grep created
	set OK = 1
    else
	echo Couldn\'t find kira...
    endif

endif

# Congratulate the user and explain some details.

if (X$OK == "X1") then
cat <<EOF

Congratulations, Starlab appears to have been successfully installed!

To use the Starlab tools, please make sure that the directory

        $install_dir/bin

is in your search path and *precedes* any system files.

For programmers, the Starlab header files are installed in

        $install_dir/include/starlab

end the libraries are in

        $install_dir/lib/starlab

Developers should find that any changes made to the automake "Makefile.am"
files will cause the relevant Makefiles ot be automagically regenerated.

                                          Enjoy!
                                          The Starlab Team

EOF
endif