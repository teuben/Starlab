
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
        
        program simin 
c       
        integer*2 n,m,num,jl,jr,idic,long
        common/sim fc2/n,m,num(288),jl(288),jr(288),idic(288),
     &                 long(16000)
c       
c       Read the data in "simout" format:
c       
        read(5,*)n,m
c       
        j = 0
        l = 1
100     read(5,*)
        read(5,*)
        read(5,*,err=999,end=999)i1,i2,i3
        j = j + 1
        num(j) = i1
        jl(j) = i2
        jr(j) = i3
c
        if (num(j).gt.0) then
            idic(j) = l
            read(5,*)(long(k),k=l,l+num(j)-1)
            l = l + num(j)
        else
            idic(j) = -1
        end if
c
        go to 100
c              
999     write(6,*)n,m,j,l
        end
