
c
c	Copyright (c) 1986,1987,1988,1989,1990,1991,1992,1993,
c       by Steve McMillan, Drexel University, Philadelphia, PA.
c
c       All rights reserved.
c
c	Redistribution and use in source and binary forms are permitted
c	provided that the above copyright notice and this paragraph are
c	duplicated in all such forms and that any documentation,
c	advertising materials, and other materials related to such
c	distribution and use acknowledge that the software was developed
c	by the author named above.
c
c	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
c	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
c	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
c

c	UNIX.F: UNIX lookalikes or replacements, for non-UNIX (or
c       	brain-dead UNIX) systems.
c
c	(Ultimately) allow access to the following (nonstandard)
c	UNIX functions:
c
c		flush		as	myflush			x
c		getlog		as	mygetlog
c		hostnm		as	myhostname
c		idate		as	mydate
c		itime		as	mytime
c		system		as	mysystem
c		getenv		as	mygetenv		x

c==========================================================================

	integer function mysystem(command)
	character*(*) command

c	Pass the specified command to the system.

#if defined(sun)

	integer system
	mysystem = system(command)

#elif defined(hpux) || defined(linux) || defined(decunix)

	integer csystem
	mysystem = csystem(command)

#endif

	end

c--------------------------------------------------------------------------

	subroutine myhostname(name)
	character*(*) name
c
c	Return the current hostname as a character string.
c
	name = 'unknown'

#ifdef sun
	call hostnm(name)
#endif

#ifdef hpux
        name = 'hpux'
#endif

#ifdef linux
        name = 'linux'
#endif

#ifdef decunix
        name = 'decunix'
#endif

	end

c--------------------------------------------------------------------------

	subroutine mygetlog(name)
	character*(*) name
c
c	Return the current username as a character string.
c
	name = 'unknown'

#if defined(sun)

	call getlog(name)

#elif defined(hpux) || defined(linux) || defined(decunix)

        call mygetenv('LOGNAME',name)
        if (name(1:1).le.' ') name = 'unknown'

#endif

	end

c--------------------------------------------------------------------------

	subroutine mydate(dat)
        character*(*) dat
c
c	Return the current date as a 10-character string (D/M/Y).
c
#ifdef sun
	integer jdate(3)
#endif

#ifdef hpux
$NOSTANDARD SYSTEM
#endif

        dat = ' '

#ifdef sun
	jdate(1) = 0
	jdate(2) = 0
	jdate(3) = 0
	call idate(jdate)
c
c	Returned order is D:M:Y.
c
	write(dat,'(i2,''/'',i2,''/'',i4)')
     &        jdate(2),jdate(1),jdate(3)
#endif

#ifdef hpux
        call date(dat)
        dat = ' '//dat
#endif
	end

c--------------------------------------------------------------------------

	subroutine mytime(tim)
        character*(*) tim
c
c	Return the current time as a 10-character string (H:M:S).
c
#ifdef sun
	integer jtime(3)
#endif

#ifdef hpux
$NOSTANDARD SYSTEM
#endif

        tim = ' '

#ifdef sun
	jtime(1) = 0
	jtime(2) = 0
	jtime(3) = 0
	call itime(jtime)
c
c	Returned order is H:M:S.
c
	write(tim,'(2x,i2.2,'':'',i2.2,'':'',i2.2)')
     &        jtime(1),jtime(2),jtime(3)
#endif

#ifdef hpux
        call time(tim)
        tim = '  '//tim
#endif

	end

c--------------------------------------------------------------------------

	subroutine myflush(unit)
	integer unit
c
c	Flush the specified logical unit
c
#if defined(sun)
	call flush(unit)

#elif defined(hpux) || defined(linux) || defined(decunix)

c
c	The best we can do?
c
	if (unit.eq.6) call cflush()

#endif

	end

c--------------------------------------------------------------------------

	subroutine mygetenv(name, value)
	character*(*) name,value

c	Return the value associated with the specified environment
c	variable, or a blank string.
c
c	Problem: some systems have only a C version, others only a
c 	FORTRAN version!

	value = ' '

#if defined(sun)

c
c	Use the FORTRAN routine:
c
	call getenv(name,value)

#elif defined(hpux) || defined(linux) || defined(decunix)

c
c	Use the C version:
c
	call cgetenv(name,value)

#endif
	end
