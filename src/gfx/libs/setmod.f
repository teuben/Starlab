
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
        
        subroutine setmod(i1,i2,i3,i4,i5)
        save
c
        logical set,mset
c
        common/fr draw/md/fr bnds/mb/fr lbx/mx/fr lby/my/fr rotn/mr
        data md/0/mb/1/mx/0/my/0/mr/0/set/.false./
c
        md=i1
        mb=i2
        mx=i3
        my=i4
        mr=i5
        set = .true.
        return
c       
        entry getmod(k1,k2,k3,k4,k5)
        k1=md
        k2=mb
        k3=mx
        k4=my
        k5=mr
        return
c       
        entry getmset(mset)
        mset = set
c
        end
