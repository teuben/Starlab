c
c       Copyright (c) 1986,1987,1988,1989,1990,1991,1992,1993,
c       by Steve McMillan, Drexel University, Philadelphia, PA.
c
c       All rights reserved.
c
c       Redistribution and use in source and binary forms are permitted
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
c
c				      **********************************
c				      *****  PostScript routines:  *****
c				      **********************************
c
        subroutine psinit(psfile in, iorient, isparc)
        save
	character *(*) psfile in

c	*****************************************************************
c	*								*
c	*	PSINIT:   Initialize PostScript output file.		*
c	*								*
c	*****************************************************************

	character*80 device
	common /plot device/ device,aspect,idev
        common /dev init/ init

	common /dev status/ idevon,idevpen,idevwt
	common /ps strokes/ nstroke,nstrpage,nstroketot
	common /ps keep/ ikeepps,iprintps
	common /ps pages/ npage
	common /ps copies/ ncopies
	common /ps bounding box/ ixleft,iybot,ixright,iytop
	common /ps write/ iwrite

	character psfile*200,psdir*20,id*30,time*20,status*7,
     &		  lpr*50,print*100,newfile*50,line*132

        integer check_open

	character*40 font
	common /ps font/ font,ifsize
	save /ps font/

	common /ps output file/ psfile,psdir
	common /psfile open/ ipsopen
	common /ps orient/ ior
        common /sparcbug/ isp

	character*200 temp

        parameter (MAX_LINES = 20000, MAX_PS_FILES = 1000)

	save /ps output file/,npsf,npsdir,
     &       /ps orient/,/sparcbug/,npsfiles,npsfin

	data npsfiles/-1/ipsopen/0/isp/0/

	if (iorient.gt.0) ior = iorient
        if (isparc.eq.0.or.isparc.eq.1) isp = isparc

c-----------------------------------------------------------------------

	npsf = len(psfile in)
	npsfin = npsf
	if (npsfin.eq.1.and.psfile in(1:1).eq.' ') npsfin = 0
	psfile(1:npsf) = psfile in(1:npsf)

        psdir = './'
        npsdir = 2

	status = 'unknown'
        iflag = 0

c       (Re)open an output file.

10	if (npsfin.le.0) then

c           Create a standard file name.  On return, if iret = 0, the
c	    file psfile will be legal (i.e. we can open it) in the
c	    directory psdir.  Otherwise, it means that we cannot find
c	    a legal filename/directory pair to use.

            call get_std_name(npsfiles, MAX_PS_FILES,
     &              	      psdir, npsdir, psfile, npsf,
     &              	      status, iret)
            if (iret.eq.1) return

	end if

        if (iflag.gt.0) then
            call make_filename(psdir, npsdir, psfile, npsf, temp, nt)
	    write(6,'(a,a)')'New output file name is ',temp(1:nt)
        end if

c       Try to open the file.

11	io = check_open(42, psdir, npsdir, psfile, npsf,
     $                  status, 'formatted')

	if (io.ne.0) then

c           If a filename was specified (npsfin > 0), try to open
c           it in /tmp, then move to the default (standard names,
c	    starting in '.').

c           If no particular name was specified (npsfin = 0), this
c	    should *not* occur, as get_std_name should have done
c	    all necessary checks.

	    if (npsfin.gt.0) then

c               A file name was specified and we cannot open it.
c		Try /tmp before switching to "standard" names.

                call make_filename(psdir, npsdir, psfile, npsf,
     $                  	   temp, nt)

                if (psdir(1:npsdir).eq.'./'
     $                  .and. psfile(1:1).ne.'/') then

                    write(6,'(''Unable to open '',a,'//
     &                      ''': trying /tmp.'')') temp(1:nt)
                    npsdir = 5
                    psdir(1:npsdir) = '/tmp/'
                    go to 11

                else

                    write(6,'(''Unable to open '',a'//
     &                      ','': using default.'')') temp(1:nt)
                    npsdir = 2
                    psdir(1:npsdir) = './'

c                   Switch to the default track.

                    npsfin = 0

                end if

		go to 10

            else

		write(6,'(a)')'Can''t open a file in "." or "/tmp"!'
		return

            end if

	end if

	ipsopen = 1
	rewind 42

c	The PostScript output file (unit 42) is now

c			psdir(1:npsdir)//psfile(1:npsf).

c-----------------------------------------------------------------------

	entry ps header

c	UNIX call!
c	---------

	call mygetlog(id)

	font = 'Times-Roman'
	ifsize = 10

c	More UNIX calls!
c	---------------

	call mydate(time(1:10))
	call mytime(time(11:20))

	write(42,100)font,time,id,font,ifsize,
     &			idevpen/256.,.1*idevwt,max(1,min(999,ncopies))
100	format(	'%!PS-Adobe-3.0 EPSF-3.0'/
     &		'%%Title: MCPAK PostScript output'/
     &		'%%DocumentFonts: ',a40/
     &		'%%Creator: Steve McMillan, Drexel University'/
     &		'%%CreationDate: ',a20/
     &		'%%BoundingBox: ',17x,'(atend)'/
     &		'%%Pages: (atend)'//
     &		'%%For: ',a30/
     &		'%%EndComments'/
     &		//'%%%% Begin mcpak setup'/
     &		'/',a40,' findfont ',i3,' scalefont setfont'/
     &		f5.3,' setgray'/f5.1,' setlinewidth'/
     &		'2 setlinejoin'/'2 setlinecap'/
     &		'/m {moveto} def'/'/l {lineto} def'/
     &		'/#copies ',i3,' def')

c    &		//'%%%% Begin mcpak setup'/'gsave'/

       if (isp.eq.0) then

c           Plot a point, on most PostScript devices...

            write(42,101)
101         format('/p {moveto -.125 -.125 rmoveto .25 0 rlineto',
     &             ' 0 .25 rlineto -.25 0 rlineto 0 -.25 rlineto',
     &             ' .125 .125 rmoveto} def')
c            write(42,101)
c101         format('/p {moveto 0 0 rlineto} def')
        else

c           The following definition draws a small box when a point
c           is requested.  This is for the benefit of the SPARCprinter,
c           which refuses to draw a line from a point to itself!

            write(42,102)
102         format('/p {moveto -.5 -.5 rmoveto 1 0 rlineto 0 1 rlineto',
     &		   ' -1 0 rlineto 0 -1 rlineto .5 .5 rmoveto} def')
        end if

c	Another magic number!

	if (ior.eq.2) write(42,110)
110	format( '610 0 translate 90 rotate')

	write(42,120)
120	format( '%%%% End mcpak setup'//
     &		'%%EndProlog'//'%%Page: 1 1')

	nstroke = 0
	nstrpage = 0
	nstroketot = 0
	npage = 1
	ixleft = 10000
	iybot = 10000
	ixright = 0
	iytop = 0

	return

c-----------------------------------------------------------------------

	entry ps check(ich)
	write(6,*)ich,npsf,psdir(1:npsdir)//psfile(1:npsf)
	return	

c-----------------------------------------------------------------------

	entry psquit(iquit)

c	Print PostScript output (UNIX only!) and terminate/reinitialize
c       plotting.

c       iquit = 1 ==> print and reinitialize PS
c       iquit = 2 ==> print and terminate PS
c       iquit = 3 ==> print prior to quitting graphics

	if (ipsopen.eq.0) return

	call ps page(0)

	write(42,'(''%%Trailer'')')

c       These two lines would be OK if (atend) were universally recognized.

c	write(42,'(''%%BoundingBox: '',4i6)')ixleft,iybot,ixright,iytop
c	write(42,'(''%%Pages: '',i7)')npage

c       However, psfig (for example) doesnt understand it, so rewrite
c       the earlier line (note fixed format and space left previously).

c       NOTE 1: The procedure of simply rewriting the BoundingBox line
c       	and scanning to the end of the file DOES NOT WORK on
c		(at least) Solaris 2.x.  Better to do it the clean way,
c		by opening a new file and copying the data...

c	NOTE 2: THE BOUNDING BOX IS WRONG FOR LANDSCAPE ORIENTATION!!

c       Dont rewrite very long files, unless explicitly asked to do so.
c 	Also, it is not clear what the bounding box means for multiple-page
c 	files, so never rewrite those.

        ieps = 0

        if (npage.eq.1.and.nstroketot.gt.0
     &          .and.(iwrite.eq.1.or.nstroketot.le.MAX_LINES)) then

            ieps = 1
            write(6,'(/a)')'Rewriting PostScript file to place'
     &                  //' BoundingBox information at start'
            rewind 42

#ifndef solaris2

c           If it is OK to rewrite a single line, do so...

            i = 0
            do while (i.eq.0)
                read(42,'(a)',iostat=i)line
                if (line(1:14).eq.'%%CreationDate') then
                    write(42,'(''%%BoundingBox: '',4i6)')
     &                      ixleft,iybot,ixright,iytop
                    write(42,'(''%%Pages: '',i7)')npage
                    do while (i.eq.0)
                        read(42,'(a)',iostat=i)line
                    end do
                end if
            end do

#else

c           Otherwise, we have to copy the entire file...

            call get_std_name(npsfiles, MAX_PS_FILES,
     &              	      psdir, npsdir, newfile, newpsf,
     &              	      status, iret)

            if (iret.eq.1) then
                write(6,'(a)')'Unable to rewrite file...'
            else

                io = check_open(43, psdir, npsdir, newfile, newpsf,
     $                  	'new', 'formatted')

                if (io.ne.0) then
                    write(6,'(a)')'Unable to rewrite file...'
                else
                    write(6,'(a)')'Copying output file '//
     &                      psfile(1:npsf)//' to '//newfile(1:newpsf)
                    i = 0
                    do while (i.eq.0)
                        read(42,'(a)',iostat=i)line
                        if (line(1:14).eq.'%%CreationDate') then
                            write(43,'(''%%BoundingBox: '',4i6)')
     &                              ixleft,iybot,ixright,iytop
                            read(42,'(a)',iostat=i)line
                            write(43,'(''%%Pages: '',i7)')npage
                        else
                            j = 132
                            do while (line(j:j).le.' '.and.j.gt.1)
                                j = j - 1
                            end do
                            write(43,'(a)')line(1:j)
                        end if
                    end do
                end if

c               Now the desired file is open as unit 43.
c               Reopen it as unit 42.

                call reopen_file(npsfin, 42, 43, io,
     &                  	 npsfiles, MAX_PS_FILES,
     &              	         psdir, npsdir, psfile, npsf,
     &                           newfile, newpsf, 1)
            end if

#endif

        else

            if (nstroketot.gt.0) write(6,'(a)')
     &              'Note: BoundingBox information is at end'
     &              //' of PostScript file.'
            write(42,'(''%%BoundingBox: '',4i6)')
     &              ixleft,iybot,ixright,iytop
            write(42,'(''%%Pages: '',i7)')npage

        end if

c	write(42,'(''grestore''/''/#copies 1 def'')')

	if (nstroketot.gt.0) then
	    close(42)
	else
	    close(42,status='delete')
	end if

c       Print the PostScript file and optionally delete it.

400	if (nstroketot.gt.0.and.iprintps.eq.1) then

c	    Yet more UNIX!!
c	    --------------

	    call mygetenv('PS_PRINTER',lpr)
	    if (index(lpr,'NO_PRINT').gt.0) then
		write(6,*)'No printer available.'
		iprintps = 0
		ikeepps = 1
		go to 400
	    end if

	    do 500 j = 50,1,-1
		if (lpr(j:j).gt.' ') then
		    l = j
		    go to 501
		end if
500	    continue
	    l = 0

501	    continue

c	    Deal with different printer names on different (UNIX) systems:
c 	    -------------------------------------------------------------

#define ___SYSV___
#ifndef hpux
#ifndef solaris2
#undef ___SYSV___
#endif
#endif

c	    	  Note: "mysystem" will pick up user-defined aliases.

c		  Explicitly suppress aliases for lp and lpr here.

c                 Not all UNIX f77s interpret backslash the same way,
c                 so set it explicitly by ASCII code (92) here:

            print(1:1) = char(92)
#ifdef ___SYSV___
	    np = 3
	    print(2:np) = 'lp'
#else
	    np = 4
	    print(2:np) = 'lpr'
#endif

	    if (l.gt.0) then

c		Use the print queue specified in the environment.

#ifdef ___SYSV___
		print(np+1:np+4) = ' -d '
		np = np + 4
#else
		print(np+1:np+4) = ' -P'
		np = np + 3
#endif
		print(np+1:np+l) = lpr(1:l)
		np = np + l
	    else
		lpr(1:1) = ' '
	    end if

c	    Using pointers is the default on System V, so suppress that
c	    here.  Also suppress the "request id is..." message.

c	    Avoiding banner pages seems to be very  system and printer-
c           dependent!
c	    To prevent a banner page, attach "-h" or "-o nobanner" to the 
c	    environment variable "PS_PRINTER".

#ifdef ___SYSV___
	    print(np+1:np+19) = ' -c -s '
	    np = np + 7
#endif

	    print(np+1:np+1) = ' '
	    np = np + 1

c           Print the file.

            if (ieps.eq.0) then
                write(6,
     &          '(''Sending PostScript output to printer '',a,''...'')')
     &          lpr(1:max(1,l))
            else
                write(6,
     &          '(''Sending Encapsulated PostScript output to printer'
     &                  //' '',a,''...'')')
     &          lpr(1:max(1,l))
            end if

            write(6,'(a)')'Print command is: '//
     &              print(1:np)//psdir(1:npsdir)//psfile(1:npsf)

	    call mysystem(print(1:np)//psdir(1:npsdir)
     &                    //psfile(1:npsf)//'; exit')

	    if (ikeepps.ne.0) then
		write(6,
     &		'(''PostScript output saved in '',a)')
     &                  psdir(1:npsdir)//psfile(1:npsf)
	    else

c		Delete the PostScript file.

                io = check_open(42, psdir, npsdir, newfile, npsf,
     $                  	'unknown', 'formatted')
		close(42, status='delete')
	    end if
	else
	    if (nstroketot.gt.0) then
                call make_filename(psdir, npsdir, psfile,
     $                  	   npsf, temp, nt)
                if (ieps.eq.0) then
                    write(6,'(''PostScript output saved in '',a)')
     &                      temp(1:nt)
                else
                    write(6,'(''Encapsulated PostScript output '
     &                      //'saved in '',a)')
     &                      temp(1:nt)
                end if
            end if
	end if	

c
c       Decide on further action.

        if (iquit.eq.1) then

c           Reinitialize the PostScript file.

c           Use the existing (saved) directory parameters from the last
c           time a PS file was explicitly initialized, but dont use the
c           last file name, if one was specified.

            write(6,'(a)')'Reinitializing PostScript output.'
            status = 'new'
            npsfin = 0
            iflag = 1

            go to 10

        else

            ipsopen = 0

c           It is up to the calling routine (e.g. mcdraw) to clean up
c           the loose ends on return...
c           
        end if

	return

c-----------------------------------------------------------------------

 	entry ps rewrite

c	Improvement over "erasepage" to clear the page...

	if (ipsopen.eq.0) return

c	Try to open a new output file with a standard name.

550	call get_std_name(npsfiles, MAX_PS_FILES,
     &              	  psdir, npsdir, newfile, newpsf,
     &              	  status, iret)
        if (iret.eq.1) then
            io = 1
            go to 560
        end if

        io = check_open(43, psdir, npsdir, newfile, newpsf,
     $          	'new', 'formatted')

560	if (io.ne.0) then
	    if (npsfiles.lt.MAX_PS_FILES) go to 550

c	    Use the less efficient method in case of problems.

	    write(6,*)'Unable to open a new output file...'
	    write(42,*)'% Unable to open a new output file...'
	    write(42,'(''erasepage'')')
	else

c	    Copy the PostScript header from the old file to the new one
c	    (everything up to the line starting with "%%Page: 1 1").
c	    Locate the last "showpage" command in the old file, then copy
c	    the last page of the old file to the new one, with omissions.

	    rewind 42
	    nline = 0
	    nlpage = 0
	    ihead = 1

600	    line = ' '
	    read(42,'(a)',err = 620,end = 620)line
	    if (ihead.eq.1) then
		if (line(1:11).eq.'%%Page: 1 1') ihead = 0
		do i = 132,1,-1
		    if (line(i:i).gt.' ') then
			write(43,'(a)')line(1:i)
			go to 610
		    end if
		end do
	    end if

610	    nline = nline + 1

c	    Note expectation of "pspage" output format.

	    if (line(7:14).eq.'showpage') nlpage = nline

	    go to 600

c	    Now actually erase the page.

620	    if (nlpage.gt.0) then
		rewind 42
		nline = 0

630		read(42,'(a)',err = 650,end = 650)line
		nline = nline + 1

		if (nline.gt.nlpage.and.line(14:14).eq.'%') then
		    if (line(13:13).eq.'p'
	1		    .or.line(13:13).eq.'m'
	2		    .or.line(13:13).eq.'l') go to 630
		end if

		do i = 132,1,-1
		    if (line(i:i).gt.' ') then
			write(43,'(a)')line(1:i)
			go to 630
		    end if
		end do
		go to 630
	    end if

650	    call reopen_file(npsfin, 42, 43, io,
     &              	     npsfiles, MAX_PS_FILES,
     &              	     psdir, npsdir, psfile, npsf,
     &                       newfile, newpsf, 2)

	end if

	end


	logical function ps_open()
	common /psfile open/ ipsopen
	save

	if (ipsopen.eq.0) then
	    ps_open = .false.
	else
	    ps_open = .true.
	end if

	end


        integer function check_open(unit, dir, nd, file, nf,
     $        			    status, form)

        integer unit, nd, nf, io
        character*(*) dir, file, status, form
        character*200 temp
        
        check_open = -1
        if (nf .le. 0) return

        call make_filename(dir, nd, file, nf, temp, nt)

c       Check whether file name is absolute or relative before opening.

        open(unit, file = temp(1:nt),
     &       status = status, form = form, iostat = check_open)

        end


        subroutine get_std_name(nfiles, max_files,
     &        			dir, ndir, file, nf,
     &        			status, iret)

        implicit none

        integer nfiles, max_files, ndir, nf, iret
        character*(*) dir, file, status

        character*200 temp
        integer nt

        integer io

c       Return a directory, name, and open status for the next
c	output file.

        iret = 0

100     nfiles = nfiles + 1
        if (nfiles.ge.max_files) then

            if (dir(1:ndir).eq.'./') then

c               Switch to /tmp and start the sequence again.

                write(6,'(''Trying to open PostScript '
     &                    //'file in /tmp'')')
                ndir = 5
                dir(1:ndir) = '/tmp/'
                nfiles = 0

            else

c               Give up -- too many files!

                write(6,'(''Can''''t open an output file '//
     &                  'in . or /tmp!'')')
                iret = 1
                return

            end if
        end if

        status = 'new'
        nf = 14
        write(file(1:nf),'(''MCD_OUT.'',i3.3,''.ps'')')nfiles

c       Check that the file (1) can be opened in this directory
c		        and (2) doesnt already exist.

        call make_filename(dir, ndir, file, nf, temp, nt)
        open(49, file=temp(1:nt), status=status, iostat=io)

	if (io.ne.0) then

c           If io is nonzero, it means that either the file exists,
c           or we dont have permission to open a file in this
c           directory.  In the latter case, there is no point trying
c           999 other files before switching to /tmp!

            if (dir(1:ndir).eq.'./') then

                call make_filename(dir, ndir, file, nf, temp, nt)
                open(49, file = temp(1:nt),
     &                  status = 'old', form = 'formatted',
     &                  iostat = io)

                if (io.eq.0) then

c                   File already exists.  Keep trying.

                    close(49)
                else

c                   Force an immediate switch to /tmp.

                    write(6,'(a)')'Cannot open files in "."'

                    nfiles = max_files
                end if
            end if

            go to 100

        end if

        close(49, status='delete')

        end


        subroutine reopen_file(npsfin, iunit, junit, io,
     &          	       nfiles, max_files,
     &        		       dir, ndir, file, nf, newfile, newf,
     &			       icall)

        implicit none
        integer npsfin, iunit, junit, io, nfiles, max_files
	integer ndir, nf, newf, icall
        character*(*) dir, file, newfile

        character*200 temp
        integer nt

        character*7 status

        integer i, iret
        character*132 line

c       On entry, the desired file is newfile(1:newf), opened as unit
c       junit.  This must be closed and reopened as unit iunit.  However,
c       if a filename was specified initially, the new file must be copied
c       to the original name!  Otherwise, dont bother with the expense
c	of rewriting.

        close(iunit, status = 'delete')

660     if (npsfin.gt.0) then

c           A name was specified.  Do it the hard way...

            call make_filename(dir, ndir, file, nf, temp, nt)
            open(iunit, file=temp(1:nt),
     &           status='unknown', form='formatted', iostat=io)
            if (io.ne.0) then
                write(6,'(a,a,a,a)')'Can''t reopen ',dir(1:ndir),
     &                  file(1:nf),' (sorry)'
                write(6,'(a,a,a)')'PostScript output file is now ',
     &                    dir(1:ndir), newfile(1:newf)

c               Revert to standard naming scheme.

                npsfin = 0
                ndir = 2
                dir(1:ndir) = './'

                call get_std_name(nfiles, max_files,
     &				  dir, ndir, file, nf,
     &				  status, iret)
                if (iret.eq.1) return

                go to 660
            end if

            rewind iunit
            rewind junit

c           Copy the new file to the old one.

            write(6,'(a,a,a,a)')
     &          'Copying file ',newfile(1:newf),' to ',file(1:nf)

665         read(junit,'(a)',err = 670,end = 670)line
            do 666 i=len(line),1,-1
                if (line(i:i).gt.' ') go to 667
666         continue
            go to 665
667         write(iunit,'(a)')line(1:i)
            go to 665

670         write(iunit,*)'%Continuing...'
            close(junit, status = 'delete')

        else

c           No specific name.  Simply reopen the new file as unit iunit.

            close(junit)
            file = newfile
            nf = newf

            call make_filename(dir, ndir, file, nf, temp, nt)
            open(iunit, file = temp(1:nt),
     &           status = 'old', form = 'formatted', iostat=io)

            if (io.ne.0) then
                write(6,'(a,a,a)')'Can''t reopen ',
     &                    dir(1:ndir),file(1:nf)
                write(6,*)'Terminating in hopeless confusion...'
                stop
            end if

c           Go to the end of the file.

675         read(iunit,'(a)',err = 680,end = 680)line
            go to 675

680         if (icall.eq.2) write(iunit,'(a)')'%Continuing...'

        end if

      end


      subroutine make_filename(dir, nd, file, nf, temp, nt)
      character*(*) dir, file, temp
      integer nd, nf, nt

c     Combine a directory and a file name to make a full pathname.
c     Ignore the directory if file is an absolute file name.

      nt = 0

      if (file(1:1).ne.'/') then
          do i=1,nd
              temp(i:i) = dir(i:i)
          end do
          nt = nd
      end if

      do i=1,nf
          temp(nt+i:nt+i) = file(i:i)
      end do

      nt = nt + nf

      end
