
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

c------------------------------------------------------------------------
c
c       Histogram-plotting routines:  histogram
c                                     histo (used by histogram)
c
c------------------------------------------------------------------------

	subroutine histogram(input,istart,nin,x,y,z,nx,ny,nz,
     &                       iherr,ihmode,ihsave,iprompt,*)
        save
c       
c	Draw a histogram from z or y(x).
c       
        character*(*) input
        dimension x(*),y(*),z(*)
c       
        common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
c       
        common/draw params/roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     &                     idevset,jbox,iorig
c       
        common/frame params 1/xmin,xmax,ymin,ymax,modex,modey
        character*80 xttl,yttl
        common/frame params 2/xttl,yttl
c       
        parameter (NHISTOMAX=200)
        dimension xh(0:NHISTOMAX),yh(0:NHISTOMAX),xp(19)
	logical normal
c
        common/weights/wt(NHISTOMAX),nw
c       
	if (input(3:3).eq.'x') then
c           
	    if (nx.le.1) go to 1001
c           
c	    Already have the data in x, y, and the box is drawn.
c           
            dx2 = .5*(x(2) - x(1))
	    yh(0) = 0.
	    xh(0) = max(xl,min(xr,x(1)-1.5*dx2))
	    do 10180 i=1,nx
	        xh(i) = max(xl,min(xr,x(i)-dx2))
	        yh(i) = y(i)
10180	    continue
c
            if (.not.normal) then
                ier = 1
            else
                ier = nz
            end if
c
            call histo(xh,yh,nx+1,ihmode,ier*iherr)
	    return
c           
	end if
c       
        if (nz.le.1) go to 1001
c       
c       Must make the histogram from scratch.
c       
        call minmax(z,nz,z1,z2)
        xmin0 = z1
        xmax0 = z2
c       
c       Determine binning to use:
c       
        if (input(3:3).ne.'2') then
            if (input(3:3).eq.'n'.or.input(4:4).eq.'n') then
                normal = .true.
            else
                normal = .false.
            end if
c           
            ref = -1.e35
            io = 1
            call readrq(input(istart:nin),2,
     &                  del,ref,dum,dum,*100)
            io = 0
c
100         if (io.ne.0.or.ref.lt.-9.e34) then
                read(input(istart:nin),*,err=1001,end=1001)del
                ref = z1
            end if
            if (del.le.0.) go to 1001
c
            xmin = z1
            xmax = z2
c           
            if (modex.lt.0) then
                if (z1.le.0..or.ref.le.0.) go to 1001
                z1 = log10(z1)
                z2 = log10(z2)
                ref = log10(ref)
            end if
c           
            iref = 0.9999*(ref-z1)/del
            nh = 0
            xh(0) = ref - (iref+1)*del
            iref = iref + 2
c           
c           Notes:  xh(i) is the TOP of bin #i.
c           ref will be at the BOTTOM of bin #iref
c           
180         nh = nh+1
            if (nh.gt.NHISTOMAX) then
                if (iprompt.eq.1) write(6,'(''Too many bins.'')')
                go to 1001
            end if
            xh(nh) = xh(nh-1) + del
            if (xh(nh).le.z2) go to 180
        else
            if (del.le.0..or.nh.le.0..or.nh.gt.nhmax) go to 1001
        end if
c       
c       Bin the data.
c       
        do 182 j=0,nh
            yh(j) = 0.
182     continue
        zbar=0.
        var=0.
        do 184 i=1,nz
            zz = z(i)
            if (modex.lt.0) zz = log10(zz)
            zbar = zbar + z(i)
            var = var + z(i)**2
            j = iref + (zz-ref)/del
            j = max(1,min(nh,j))
            yh(j) = yh(j) + 1.
184     continue
c       
        if (normal) then
            do 185 j = 0,nh
                yh(j) = yh(j)/nz
185         continue
        end if
c       
c       Get statistics (percentiles from the binned data)...
c       
        zbar = zbar/nz
        var = var/nz - zbar**2
c       
        ix=0
        sum=0.
        do 190 ip=1,19
            sp = .05*ip*nz
            if (normal) sp = sp/nz
188         ix = ix + 1
            if (ix.gt.nh) go to 191
            sum = sum + yh(ix)
            if (sum.lt.sp) go to 188
            xp(ip) = xh(ix) - del*(sum-sp)/yh(ix)
            sum = sum - yh(ix)
            ix = ix - 1
190     continue
c       
c       ...and (clip and) draw the histogram.
c       
191     call minmax(yh(1),nh,y1,y2)
c       
        if (modey.gt.0) then
            ymin = 0.
        else
            if (.not.normal) then
                ymin = 1.
            else
                ymin = 1./nz
            end if
            ymin1 = log10(ymin) - 10.
            do 192 j = 0,nh
                if (yh(j).lt.ymin) then
                    yh(j) = ymin1
                else
                    yh(j) = log10(max(ymin,yh(j)))
                end if
192         continue
        end if
c       
        if (input(3:3).ne.'1'.and.input(3:3).ne.'2') then
c           
            ymax = min(10.+y2,1.1*y2)
c           
            delz = xp(19) - xp(1)
c           
            if (xp(1)-z1.gt.delz) then
c               
c               Bottom 5% too spread out:
c               
                xmin = xp(19) - sqrt(delz*(xp(19)-z1))
                ix = abs(xmin - ref)/del
                if (xmin.lt.ref) then
                    xmin = ref - (ix+1)*del
                else
                    xmin = ref + ix * del
                end if
                if (modex.lt.0) xmin = 10.**z1
            end if
c           
            if (z2-xp(19).gt.delz) then
c               
c               Top 5% too spread out:
c               
                xmax = xp(1) + sqrt(delz*(z2-xp(1)))
                ix = abs(xmax - ref)/del
                if (xmax.lt.ref) then
                    xmax = ref - ix*del
                else
                    xmax = ref + (ix+1) * del
                end if
                if (modex.lt.0) xmax = 10.**z2
            end if
c
c           Expand the range slightly (linear case only):
c
            if (modex.gt.0) then
                xav = .5*(xmin+xmax)
                xmin = xmin - .1*(xav-xmin)
                xmax = xmax + .1*(xav-xmin)
            end if
c           
            if (idevset.eq.0) call setup(' ','histogram')
            call devon
            call clear
c           
            if (.not.normal) then
                yttl = '@H'
            else
                yttl = ' fraction'
            end if
c           
            call eframe(xmin,xmax,xlen,modex,xttl,
     &                  ymin,ymax,ylen,modey,yttl)
            ibox=1
c           
c           More clipping...
c           
            do 193 j=0,nh
                if (xh(j).lt.xl) xh(j)=xl
                if (xh(j).gt.xr) xh(j)=xr
193         continue
        end if
c       
        if (.not.normal) then
            ier = 1
        else
            ier = nz
        end if
c
        call histo(xh,yh,nh+1,ihmode,ier*iherr)
        call devoff
c
        if (ihsave.ne.0) then
c
c           Save the histogram data (top center of each bar) in y(x).
c
            do 200 j=2,nh
                x(j-1) = .5*(xh(j-1)+xh(j))
                y(j-1) = yh(j)
                if (ihsave.eq.2.and.nw.gt.0) z(j-1) = wt(j)
200         continue
            nx = nh-1
            ny = nx
        end if
        if (ihsave.eq.2.and.nw.gt.0) nz = nx
c       
        if (iprompt.eq.1) then
            if (modex.lt.0) then
                do 194 ip=1,19
                    xp(ip) = 10.**xp(ip)
194             continue
            end if
            write(6,195)nz,xmin0,xmax0,zbar,
     &              sqrt(var),(xp(ip),ip=2,18,2)
195         format(i6,' points: ',
     &              ' xmin = ',1pe10.3,',  xmax = ',e10.3/
     &              15x,' mean = ',e10.3,',  s.d. = ',e10.3/
     &              ' %tiles:',5e12.3/8x,4e12.3)
        end if
c       
        return
1001    return 1
c       
        end

        
        subroutine histo(x,y,n,ihmode,iherr)
        save
c       
c       Draw y(x) as a histogram or as data points with error bars.
c       
        dimension x(*),y(*)
c
        common/frame params 1/xmin,xmax,ymin,ymax,modex,modey
        common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
        common/draw params/roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     &                     idevset,jbox,iorig
c
        parameter (NHISTOMAX=200)
        common/weights/wt(NHISTOMAX),nw
c       
        nw = 0
c
        if (iherr.le.0) then
            r = (x(1)-xl)*dinchx
            so = 0.
            do 100 i=2,n
                s = max(y(i)-ybot,0.)*dinchy
                if (ihmode.eq.0) then
                    call plot(r,0.,3)
                else
                    call plot(r,min(s,so),3)
                end if
                call plot(r,max(s,so),2)
                call plot(r,s,3)
                r = (x(i)-xl)*dinchx
                call plot(r,s,2)
                so = s
                if (r.ge.rlen) go to 101
100         continue
101         if (r.lt.rlen) call plot(r,0.,2)
        else
            dx = .5*(x(3)-x(2))
            if (hp.le.0.) then
                nn = 4
                hh = .025
            else
                nn = 10
                hh = .5*hp
            end if
c
            do 250 i = 2,n
                if (x(i).ge.xr.or.y(i).lt.ybot) go to 250
                r = (x(i)-dx-xl)*dinchx
                s = (y(i)-ybot)*dinchy
                call ngon(r,s,hh,nn,0.)
c
c               Find the actual number of points in the bin, zbin.
c
                zbin = y(i)
                if (modey.lt.0) zbin = 10.**zbin
                zbin = zbin*iherr
c
c               Determine the standard error.
c
                nw = nw + 1
                wt(i-1) = sqrt(zbin)
c
c               (Note that wt(1) corresponds to x(2)!)
c
                zfac = 1./wt(i-1)
                if (zbin.le.1.) zfac = .999
c
                if (modey.gt.0) then
                    su = (y(i)*(1.+zfac) - ybot)*dinchy
208                 sl = (y(i)*(1.-zfac) - ybot)*dinchy
                else
                    su = (y(i) + log10(1.+zfac) - ybot)*dinchy
                    sl = (y(i) + log10(1.-zfac) - ybot)*dinchy
                end if
c
c               Plot the error bars.
c
                call plotin(r-hh,sl,3)
                call plotin(r+hh,sl,2)
                call plotin(r,sl,3)
                call plotin(r,su,2)
                call plotin(r-hh,su,3)
                call plotin(r+hh,su,2)
c
250         continue
c
        end if
c       
        end
