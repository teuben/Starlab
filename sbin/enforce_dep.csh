#!/bin/csh -f

# Arguments are divided into "sources" and "dependencies",
# separated by a "@".

# Enforce dependencies not easily handled by make by touching
# any source older than one of the dependencies.

set src = ()
set dep = ()
set which = 0

foreach arg ($argv[*])

    if ("X$arg" == "X@") then

	which = 1

    else

	if ($which == 1) then
	    set src = ($src $arg)
	else
	    set dep = ($dep $arg)
	endif

    endif

end

if ($#src == 0) exit
if ($#dep == 0) exit

foreach s ($src)
    if (`filedate $s` < `filedate $dep`) then
	touch $s
    endif
end
