
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


      subroutine enumbr(xi,yi,ht,fpn,theta,ll)
        save
c
c     leaves xi and yi unchanged.
c     plots a number of the form 7.23(-6).
c
      call compoz(fpn,fpnred,lpow)
      call fr numbr(xi,yi,ht,fpnred,theta,ll)
      call simwhe(xn,yn)
      call fr symbl(xn,yn,ht,'(',theta,1)
      call simwhe(xn,yn)
      call fr numbr(xn,yn,ht,float(lpow),theta,-1)
      call simwhe(xn,yn)
      call fr symbl(xn,yn,ht,')',theta,1)
      end
