
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

	subroutine printid
        save
c
c	Put identifying information at top of page.
c
	parameter (lid=150)
	character*150 id
	character*11 page
        character*10 tim, dat
        common /plot origin/ ro,so
        common /dev status/ idevon,idevpen,idevwt
	common /ps pages/ npage
c
	character*40 font
	common /ps font/ font,ifsize
	save /ps font/
c
	common /ps head flag/ iheadflag
c
	common /ps orient/ior
c
	save /ps head flag/,/ps orient/
c
c	Draw any outstanding lines:
c
	call ps stroke
c
	if(iheadflag.eq.0)return
c
c	Set up the identifying strings:
c
	id=' '
	nid=6
	id(1:nid+2)='User:   '
c
	call mygetlog(id(nid+3:lid))
	do 10 i=lid,nid+2,-1
	    if(id(i:i).gt.' ')go to 20
10	continue
20	nid=i+3
	id(nid:nid)='('
c
	call myhostname(id(nid+1:lid))
	do 30 i=lid,nid,-1
	    if(id(i:i).gt.' ')go to 40
30	continue
40	nid=i+11
	id(nid-10:nid)=').         '
c
	call mydate(dat)
	write(id(nid+1:nid+27),'(''Date:  '',a10,''.'',9x)')dat
	nid=nid+27
c
	call mytime(tim)
	write(id(nid+1:nid+18),'(''Time:  '',a10,''.'')')tim
	nid=min(lid,nid+18)
c
	write(page,'(''Page'',i7)')npage
c
c	Set header font, pen color and weight (weight probably unnecessary):
c
	write(42,'(/''/Times-Roman findfont 12 scalefont setfont''/)')
	ipsto=idevpen
	call color(0)
	iwsto=idevwt
	call weight(7)
c
c	Draw the header:
c
c	Some more magic numbers!
c
	if (ior.eq.1) then
	    ix = 50
	    jy = 750
	    kx = 525
	else
	    ix = 50
	    jy = 580
	    kx = 725
	end if
c
	write(42,50)ix,jy,(id(i:i),i=1,nid)
50	format(/2i4,' m ('/100a1)
	write(42,51)kx,jy,page
51	format(	') show ',2i4,' m (',a11,') show')
c
c	Reset font, pen color and weight.
c
	write(42,60)font,ifsize
60	format('/',a20,' findfont ',i3,' scalefont setfont'/)
	call color(ipsto)
	call weight(iwsto)
c
	return
c
	entry noheader
c
	iheadflag=0
c
	return
c
	entry header
c
	iheadflag=1
c
	end
