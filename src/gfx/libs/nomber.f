
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

      subroutine nomber (xcall,ycall,ht,f,an,ic)
      save
c     
c     Similiar to calcomp number (calling sequence identical),
c     but uses simbol to plot the figures.
c     
c     Created  20.jun.80  plh;
c     Re-written  5.mar.82, 21.may.82  smcm.
c     
      integer*4 i,ic,j,m,nchar
      character outstr*80,outstr1*160
      common/fr int/iframe
      common/str posn/i pos set,frx,fry
c     
      character*1 pct,blank
      data pct/'%'/,blank/' '/
c     
      inom=1
      x=xcall
      y=ycall
      go to 1
c     
      entry unomber(xcall,ycall,ht,f,an,ic)
      entry usernomber(xcall,ycall,ht,f,an,ic)
      inom=1
      call fr inches(xcall,ycall,x,y)
      go to 1
c     
      entry nombr(xcall,ycall,ht,f,an,ic)
      inom=0
      x=xcall
      y=ycall
      go to 1
c     
      entry unombr(xcall,ycall,ht,f,an,ic)
      entry usernombr(xcall,ycall,ht,f,an,ic)
      inom=0
      call fr inches(xcall,ycall,x,y)
c     
1     call numsym(f,ic,outstr,m)
      nchar=0
c     
c     Ensure numerical labels are centered in this case:
c     
      if (iframe.eq.1.and.i pos set.eq.0) nchar=-999
c     
c     Add %% for simbol:
c     
      outstr(m+1:m+2)=pct//pct
c     
      if (f.lt.0.) then
          if (nchar.lt.0) then
c             
c             Extra blank for better centering of negative numbers:
c             
              m=m+1
              outstr(m+2:m+2)=outstr(m+1:m+1)
              outstr(m+1:m+1)=outstr(m:m)
              outstr(m:m)=blank
          end if
      end if
c
      if (inom.eq.0) then
          outstr1(1:m+2)=outstr(1:m+2)
          i=0
          do j=1,m+2
              i=i+1
              if (outstr1(j:j).lt.'0'.or.outstr1(j:j).gt.'9') then
                  outstr(i:i)=outstr1(j:j)
              else
                  outstr(i:i)='@'
                  i=i+1
                  outstr(i:i)=outstr1(j:j)
              end if
          end do
      end if
c
      call simbol (x,y,ht,outstr,an,nchar)
c     
      end
