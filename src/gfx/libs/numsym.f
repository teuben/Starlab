
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

      subroutine numsym(f,ic,outstr,m)

      implicit double precision (d)
      save
      real*8 spacc
      real*4 f
      integer*4 i,ic,j,m,n
      character outstr*80

      character*1 dot,neg,blank,digit(0:9)
      data dot/'.'/neg/'-'/,blank/' '/
     &        digit/'0','1','2','3','4','5','6','7','8','9'/
     &        spacc/1.00000001/

      dex(dpf) = 10.**dpf

      dpf = dble(f)
      n = 0
      if (dpf.ne.0.) then
          n = log10(abs(dpf)) + 1.001		! digits to left of decimal.
      end if
      if (n.le.0) n = 1
      m = n + 1 + ic                  		! total characters plotted.
      dpg = abs(dpf)             		! for discounting.

      if (ic.lt.0) then

         do i = 1,m
            dph = dex(dfloat(n-i))
            j = dpg/dph*spacc   		! discounted digit.
            outstr(i:i) = digit(j)
            dpg = dpg - j*dph     		! discount dpg for next digit.
         end do
         
      else

         do i = 1,n
            dph = dex(dfloat(n-i))
            j = dpg/dph*spacc   		! discounted digit.
            outstr(i:i) = digit(j)
            dpg = dpg-j*dph     		! discount dpg for next digit.
         end do

         outstr(n+1:n+1) = dot  		! insert decimal point.
         if (ic.gt.0) then
            do i = 1,ic          		! add decimal digits.
               j = 10*dpg*spacc  	 	! discounted digit.
               outstr(n+1+i:n+1+i) = digit(j)
               dpg = 10.*dpg-j   	 	! continue discounting dpg.
            end do
         end if

      end if

      j = 10*dpg*spacc
      if(j.lt.5)go to 501       		! round off last digit(s).

      if(ic.gt.0)then
          do i = m,n+2,-1   			! start at right edge of string.
              if (ichar(outstr(i:i)).le.56) then
                  outstr(i:i) = char(ichar(outstr(i:i))+1)
                  go to 501
              end if
              outstr(i:i) = digit(0)
          end do
      end if

      do i = n,1,-1         		    ! then adjust left, if necessary.
          if (ichar(outstr(i:i)).le.56) then
              outstr(i:i) = char(ichar(outstr(i:i))+1)
              go to 501
          end if
          outstr(i:i) = digit(0)
      end do
      do  i = m,1,-1
          outstr(i+1:i+1) = outstr(i:i)
      end do
      outstr(1:1) = digit(1)

c     Put "-" up front if necessary.

501   if (f.lt.0.) then
          do i = 1,m
              if (outstr(i:i).ne.'.'.and.outstr(i:i).ne.'0') go to 506
          end do
          go to 508
506       do i = m+2,1,-1
              outstr(i+1:i+1) = outstr(i:i)
          end do
          outstr(1:1) = neg
          m = m+1
      end if
508   continue

      end
