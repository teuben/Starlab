
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
c
c     Input and output routines to mimic VAX extensions.
c
      subroutine readin(lun,string,nl,*)
      save
c
c     "read(lun,'(q,a200)')nl,string"
c
      character*(*) string
c
      string=' '
      read(lun,'(a200)',end=1000,err=1000)string
c
      do 1 nl=200,1,-1
          if(string(nl:nl).gt.' ')return
1     continue
      return
c
1000  return 1
      end
