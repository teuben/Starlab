
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

	subroutine getline(inunit,line,nl,iprompt,nhist)
        save
c
	character*(*)line
	character*80 scratch,format
c
        character*80 device
        common/plot device/device,aspect,idev
        common /x input/ interact
c
        character*20 version
        common /mcd version/ nversion,version
c
	integer eof
c
	eof = 0
c
12	if (inunit.eq.5.and.iprompt.eq.1) then
	    write(scratch(1:6),'(i6)')nhist+1
  	    do 10012 i=1,6
      	        if (scratch(i:i).gt.' ')go to 10112
10012       continue
c
10112	    nfmt = 24-i+nversion
            format(1:nfmt) = '(''mcdraw-'//version(1:nversion)
     &                                  //'('//scratch(i:6)//'): ''$)'
  	    write(6,fmt=format(1:nfmt))
	end if
c
	n_xopen = num_win(17)
c
c       If n_xopen = 0, there aren't any X-windows to get input from, so
c       just default to standard FORTRAN input.  Note that we can get
c       input from X even if another device is current.
c

        if (inunit.ne.5.or.n_xopen.le.0.or.interact.eq.0) then
c
            call readin(inunit,line,nl,*20)
c
        else
c
c           Get input via X, keeping screen up to date.
c
            call myflush(6)
            line = '\0'
c
            call win_read_line(line)
c
            do 15 nl=len(line),1,-1
                if (line(nl:nl).gt.' ') go to 16
15          continue
            nl = 0
16          continue
        end if
c
c	The entire line is line(1:nl).
c	Comments start with a '#', all other lines must end in a blank.
c
  	if (nl.le.0.or.line(1:1).eq.'#') go to 12
c
  	nl=nl+1
  	line(nl:nl)=' '
	return
c
20	if (inunit.ne.5) then
	    close(inunit)
	    inunit = 5
	else
c
c	    For now, an interactive ^D will terminate the program.
c
	    if (eof.ge.0) then
		line = 'q'
		nl = 1
		return
	    end if
	    eof = eof + 1
	end if
	go to 12
c
	end
