
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

c       Contents:  sdecode -- decode a string into one or more integers
c                  decode  -- decode a character into an integer.
c                  rdecode -- decode a string into a range of integers.

        subroutine sdecode(string,narg,iarg,*)
c
c       Decode narg (< =  3) arguments from the given string.
c       Delimiters are " " or ",".
c
c       The string itself is unaltered on return.
c
        character*(*) string
        integer iarg(3)
c
        if (narg.gt.3) return
c
        ibl = 1
        jarg = 0
        do i = 1,len(string)
            if (string(i:i).eq.' '.or.string(i:i).eq.',') then
                ibl = 1
            else
                if (ibl.eq.1) then
                    jarg = jarg+1
c
                    if (i.eq.len(string).or.string(i+1:i+1).eq.' '
     &                      .or.string(i+1:i+1).eq.',') then
                        call decode(string(i:i),iarg(jarg),*50)
                        go to 75
                    end if
c
c                   Attempt to interpret an illegal character or a longer
c                   string as an integer.
c
50                  call readiq(string,i,len(string),1,
     &                          ia,idum,idum,idum,*1001)
                    iarg(jarg) = ia
c
75                  if (jarg.ge.narg) return
                    ibl = 0
                end if
            end if
        end do
c
        do j = jarg+1,narg
            iarg(j) = 0
        end do
c
        return
1001    return 1
c
        end


        subroutine decode(ch,i,*)
c
c       Return 1, 2, or 3 as the code for the input character.
c
        character*1 ch
c
        if (ch.eq.'x'.or.ch.eq.'X'.or.ch.eq.'1') then
            i = 1
        else if (ch.eq.'y'.or.ch.eq.'Y'.or.ch.eq.'2') then
            i = 2
        else if (ch.eq.'z'.or.ch.eq.'Z'.or.ch.eq.'3') then
            i = 3
        else
            return 1
        end if
c
        end
        

        subroutine rdecode(str,nh,i1,i2,*)
c
c       Extract numbers from the input string.  String format may be
c       i1:i2 or i1#i2 or i1^i2.  In the event of an error reading i2,
c       return i2 = i1 (this allows a dual function for this routine).
c
c       The value of str is not altered.  The (unchanged) variable
c       nh is used to set an upper limit on the values returned.
c
        character*(*) str
c       
        l = len(str)
c
        do if = 1,l
            if (str(if:if).gt.' ') go to 2
        end do
        return 1
c
2       do il = l,if,-1
            if (str(il:il).gt.' ') go to 4
        end do
c
c       String runs from if to il.
c
4       do i = if,il
            if (str(i:i).eq.':'.or.str(i:i).eq.'^'
     &              .or.str(i:i).eq.'#') go to 20
        end do
        i = il+1
c
c       Internal delimiter is at location i.
c
20      if (i.gt.if) then
            read(str(if:i-1),*,err = 999,end = 999)i1
        else
            i1 = 1
        end if
c
        if (i.lt.il) then
            read(str(i+1:il),*,iostat = io)i2
            if (io.ne.0) i2 = i1
        else if (i.eq.il) then
            i2 = nh
        else
            i2 = i1
        end if
c
        if (i1.lt.0.or.i2.lt.0)
     &          write(6,*)'Warning: non-relocatable '//
     &                    'historical reference'

        if (i1.lt.0) i1 = nh+1+i1
        if (str(i:i).eq.'^'.or.str(i:i).eq.'#') then
            i2 = abs(i2)
            i2 = i1+i2-1
        end if
        if (i2.lt.0) i2 = nh+1+i2
        i1 = max(0,min(nh,i1))
        i2 = min(max(i1,i2),nh)
c
        return
c
999     return 1
c       
        end
