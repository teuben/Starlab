
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

        subroutine powerspectrum(t,x,n,iout)
c
c       Replace x by its power spectrum, t by frequency.
c       Note that it is assumed that t is uniformly spaced.
c
c       Since we will usually want to take the logarithm on return,
c       SUPPRESS the zero-frequency component.
c
        real t(n),x(n)
c
        if (n.le.1) return
c
c       Discard data to get to a power of 2.
c
        i2 = log(float(n))/log(2.)
        n = nint(2.**float(i2))
        if (iout.ne.0) write(6,*)'Array length truncated to ',n,
     &          ' for FFT'
c
c       Determine the frequency (note storage convention):
c
        dti = 1./(t(n) - t(1))
        do i=1,n/2
            t(i) = i*dti
        end do
c
c       Take out the mean signal.
c
        sum = 0.
        do i=1,n
            sum = sum + x(i)
        end do
c
        sum = sum/n
        do i = 1,n
            x(i) = x(i) - sum
        end do
c
c       Fourier transform (note packing of result).
c
c       Initial data are x(i), i = 1,...,n
c       Results are x(1) = F(0), x(2) = F(n/2), x(3) + i x(4) = F(1), etc.
c
        call realft(x,n/2,1)
c
c       Replace x by the power spectrum.
c
c       Zero-frequency component is x(1), note, but suppressed here:
c
c        x(1) = x(1)**2
c
        x2 = x(2)**2
c
        do 50 k=2,n/2
            x(k-1) = x(2*k-1)**2 + x(2*k)**2
50      continue
c
        x(n/2) = x2
        n = n/2
c       
99999   end
        
        
        subroutine realft(data,n,isign)
c
c       Routine from Numerical Recipes.
c
        real*8 wr,wi,wpr,wpi,wtemp,theta
        dimension data(2*n)
c
        theta=6.28318530717959d0/2.0d0/dfloat(n)
        c1=0.5
        if (isign.eq.1) then
            c2=-0.5
            call four1(data,n,+1)
        else
            c2=0.5
            theta=-theta
        endif
        wpr=-2.0d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.0d0+wpr
        wi=wpi
        n2p3=2*n+3
        do 95000 i=2,n/2+1
            i1=2*i-1
            i2=i1+1
            i3=n2p3-i2
            i4=i3+1
            wrs=sngl(wr)
            wis=sngl(wi)
            h1r=c1*(data(i1)+data(i3))
            h1i=c1*(data(i2)-data(i4))
            h2r=-c2*(data(i2)+data(i4))
            h2i=c2*(data(i1)-data(i3))
            data(i1)=h1r+wrs*h2r-wis*h2i
            data(i2)=h1i+wrs*h2i+wis*h2r
            data(i3)=h1r-wrs*h2r+wis*h2i
            data(i4)=-h1i+wrs*h2i+wis*h2r
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
95000   continue
        if (isign.eq.1) then
            h1r=data(1)
            data(1)=h1r+data(2)
            data(2)=h1r-data(2)
        else
            h1r=data(1)
            data(1)=c1*(h1r+data(2))
            data(2)=c1*(h1r-data(2))
            call four1(data,n,-1)
        endif
c
        end

        
        subroutine four1(data,nn,isign)
c
c       Routine from Numerical Recipes.
c
        real*8 wr,wi,wpr,wpi,wtemp,theta
        dimension data(2*nn)
c
        n=2*nn
        j=1
        do 95000 i=1,n,2
            if(j.gt.i)then
                tempr=data(j)
                tempi=data(j+1)
                data(j)=data(i)
                data(j+1)=data(i+1)
                data(i)=tempr
                data(i+1)=tempi
            endif
            m=n/2
1           if ((m.ge.2).and.(j.gt.m)) then
                j=j-m
                m=m/2
                go to 1
            endif
            j=j+m
95000   continue
        mmax=2
2       if (n.gt.mmax) then
            istep=2*mmax
            theta=6.28318530717959d0/(isign*mmax)
            wpr=-2.d0*dsin(0.5d0*theta)**2
            wpi=dsin(theta)
            wr=1.d0
            wi=0.d0
            do 95001 m=1,mmax,2
                do 95002 i=m,n,istep
                    j=i+mmax
                    tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
                    tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
                    data(j)=data(i)-tempr
                    data(j+1)=data(i+1)-tempi
                    data(i)=data(i)+tempr
                    data(i+1)=data(i+1)+tempi
95002           continue
                wtemp=wr
                wr=wr*wpr-wi*wpi+wr
                wi=wi*wpr+wtemp*wpi+wi
95001       continue
            mmax=istep
            go to 2
        endif
c
        end

