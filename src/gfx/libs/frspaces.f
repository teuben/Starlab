
c     
c     Copyright (c) 1986,1987,1988,1989,1990,1991,1992,1993,
c     by Steve McMillan, Drexel University, Philadelphia, PA.
c     
c     All rights reserved.
c     
c     Redistribution and use in source and binary forms are permitted
c     provided that the above copyright notice and this paragraph are
c     duplicated in all such forms and that any documentation,
c     advertising materials, and other materials related to such
c     distribution and use acknowledge that the software was developed
c     by the author named above.
c     
c     THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
c     IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
c     WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
c     


      subroutine fr spaces(ztit,nztit)
      save
c     
c     Calculate the true length, nztit, of a character string, ztit.
c     The length of the string is the position of the character before
c     a block of five blanks or a non-printing character is encountered.
c     If the input string is blank, nztit=1 is returned.
c     
      character*(*) ztit
      character*5 space
      space='     '
c
      nztit=index(ztit,space)-1
      if (nztit.lt.0) then
          nztit=len(ztit)
      else if (nztit.eq.0) then
          nztit=1
      end if
c
      do i=1,nztit
          ich=ichar(ztit(i:i))
          if (ich.lt.32.or.ich.gt.126) go to 20
      end do
      return
20    nztit=max(1,min(i-1,nztit))
c
      end
