
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
        
        subroutine setpat(nt1,nb1,nt2,nb2)
        save
c       
c       Set pattern of dashes for dplot.
c       
        character*80 device
        common /plot device/ device,aspect,idev
        common/dash/dpatrn(10),dpat,npatrn,ipat,lpen
c       
        data blank,blank0/0.05,.025/
c       
        lpen=2
        ipat=1
        npatrn=4
        dpatrn(1)=nt1*blank
        dpatrn(2)=nb1*blank
        dpatrn(3)=nt2*blank
        dpatrn(4)=nb2*blank
        do 1 i=1,npatrn
            if (dpatrn(i).eq.0.) then
                dpatrn(i)=blank0
c               if (idev.eq.15) dpatrn(i)=2.*dpatrn(i)
            end if
1       continue
        dpat=dpatrn(1)
c       
        return
c       
        entry getpat(i1,i2,i3,i4)
c       
c       Return pattern of dashes.
c       
        i1=dpatrn(1)/blank
        i2=dpatrn(2)/blank
        i3=dpatrn(3)/blank
        i4=dpatrn(4)/blank
c       
        end
