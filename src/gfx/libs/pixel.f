
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

        subroutine pixel(iin,jin,ipenin)
        save
c
c       Move or draw to a specified pixel. Currently implemented
c       for the Tektronix options (or equivalent) only.
c
        character*80 device
        common /plot device/ device,aspect,idev
        common /framesize/ nxpix,nx0,xfac,nypix,ny0,yfac
        common /dev init/ init
        common /dev details/ itek,ivers
        logical on
        character*1 ctrl(0:31),
     &              null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
        common /ctrlch/ ctrl,
     &                  null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
c
        character*1 vec(0:5),up(3),down(3)
        data up/' ','L','F'/down/' ','L','G'/
c
        if(itek.eq.0)return
        up(1)=esc
        down(1)=esc
        vec(0)=gs
c
        if(init.eq.0)then
            init=-1
            call mcinit
            call devon
            call clear
        end if
        if(.not.on())call devon
        i=iin
        j=jin
        ipen=abs(ipenin)
        if(ipen.ne.2)ipen=3
        if(nxpix.gt.1023)then
            i4=i/4
            ii=i-4*i4
            j4=j/4
            jj=j-4*j4
            vec(2)=char(96+ii+4*jj)
            i=i4
            j=j4
        else
            vec(2)='`'
        end if
        j32=j/32
        vec(1)=char(32+j32)
        vec(3)=char(96+(j-32*j32))
        i32=i/32
        vec(4)=char(32+i32)
        vec(5)=char(64+(i-32*i32))
c
c       (these assignments are ok because none of the values are >127)
c
        if(idev.lt.13)then
            if(ipen.eq.3)then
                call type string(vec(0),6)
            else
                call type string(vec(1),5)
            end if
        else
            if(ipen.eq.3)then
                call type string(up,3)
            else
                call type string(down,3)
            end if
            call type string(vec(1),5)
        end if
        end
