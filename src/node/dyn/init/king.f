
c       Generate a King model.
c
c       Parameters:	-n number
c			-w depth
c			-i random seed
c
c       Note: Uses functions qromb, rkode, and random numbers
c             from "Numerical Recipes"

	program kingmodel

	implicit real*8 (a-h, o-z)
	parameter (NMAX = 50000)
	
	real*8 mass(NMAX), pos(NMAX, 3), vel(NMAX, 3)
        real*8 w0
        real*8 r, random
	integer n, iargc 
        integer iseed
        character*80 arg

c       Settable parameters:	w0 = dimensionless depth
c				n  = number of particles

  	iout = 1

c	Initialize random number generator using iseed.
 	iseed = 42

        n = 0
        do i=1,iargc()
            call getarg(i, arg)
            if (arg(1:2) .eq. '-w') then
                call getarg(i+1, arg)
		read(arg, '(i10)', iostat=io) iw
		if (io.eq.0) then
		   w0 = iw
		else
		   read(arg, '(f8.3)', iostat=io) w0
		   if (io.ne.0) then
		      write(6,'(a)')'Error reading w0.'
		      stop
		   end if
		end if
                call setking(w0)
            else if (arg(1:2) .eq. '-n') then
                call getarg(i+1, arg)
                read(arg, '(i10)', iostat=io) n
                if (io.ne.0) then
		    write(6,'(a)')'Error reading N.'
		    stop
		 end if
            else if (arg(1:2) .eq. '-i') then
                call getarg(i+1, arg)
                read(arg, '(i10)', iostat=io) iseed
                if (io.ne.0) then
		    write(6,'(a)')'Error reading iseed.'
		    stop
		end if
		r = random(iseed)
            end if
	 end do

        if (n .le. 0) then
	    write(6,'(a)')'Number of stars (N) not set: use "-n #".'
	    stop
	end if

	call king(mass, pos, vel, n, NMAX, iout)
        
	do i=1,n
           print *, i, mass(i), (pos(i,k),k=1,3),  (vel(i,k),k=1,3)
	end do
        
	end
        
        
        block data init_king

        implicit real*8 (a-h,o-z)

        common/kpars1/w0,v0,e0,pi,twopi
        common/kpars2/nprof
        common/therm/g(0:1000)

        data w0/0./
        data g(0)/-1./

        end
        
        
        subroutine setking(w)

        implicit real*8 (a-h,o-z)

        common/kpars1/w0,v0,e0,pi,twopi
        common/kpars2/nprof

        w0=w

        end
        
        
        subroutine king(body,x,xdot,n,nmx,iout)

c       Set up a King model, with total mass = 1, core radius = 1..

        implicit real*8 (a-h,o-z)

        dimension body(nmx),xdot(nmx,3)
        real*8 x(nmx,3)

        parameter (nm=2500)
        common/profile/rr(0:nm),d(0:nm),v2(0:nm),psi(0:nm),
     &       zm(0:nm),indx(0:100)
        common/kpars1/w0,v0,e0,pi,twopi
        common/kpars2/nprof

        if (w0.eq.0.) then
           write(6,'(a)')'Dimensionless depth of'//
     &         ' cluster potential (w0) not set: use "-w #".'
           stop
	end if

c       Compute the cluster density/velocity/potential profile.

        call poisson(rr,d,v2,psi,zm,nm,w0,nprof,v20,iout)

c	do i=0,nprof
c	    write(0,*), i, rr(i), d(i)
c	end do

c       Index the mass distribution.

        do 10 i=0,nprof
           zm(i)=zm(i)/zm(nprof)
 10     continue
        indx(0)=0
        indx(100)=nprof

c       Determine core and half-mass radii.

        dz=.01
        z=dz
        iz=1
        do 20 j=1,nprof-1
           if (rr(j).lt.1.)jcore=j
           if (zm(j).lt..5)jhalf=j
           if (zm(j).gt.z) then
              indx(iz)=j-1
              z=z+dz
              iz=iz+1
           end if
 20     continue
        rhalf=rr(jhalf)+(rr(jhalf+1)-rr(jhalf))*(.5-zm(jhalf))
     &       /(zm(jhalf+1)-zm(jhalf))
        zmcore=zm(jcore)+(zm(jcore+1)-zm(jcore))*(1.-rr(jcore))
     &       /(rr(jcore+1)-rr(jcore))

        if (iout.ne.0)write(6,30)w0,nprof,v20,rr(nprof),rhalf,zmcore
30      format(/' King model:  W0 =',f6.2,',  NPROF = ',i4/
     &          ' V20/E0 =',f7.3,', Rt/Rc =',f7.2,
     &          ', Rh/Rc =',f6.2,', Mc/M =',f7.4/)

        zmass=0.
        body(1)=1.

        do 50 i=1,n
            body(i)=body(1)
            zmass=zmass+body(i)
50      continue
        d0=zmass/zm(nprof)
        pi=4.*atan(1.)
        twopi=2.*pi
        four3pi=4.*pi/3.
        e0=four3pi*d0/v20
        v0=.004*sqrt(2.*e0)

c       Assign positions and velocities. Note that it may actually
c       be preferable to do this in layers instead.

c       Note also that rcore is taken to be unity.

        do 100 i=1,n
            call posvel(xx,yy,zz,uu,vv,ww)
            x(i,1)=xx
            x(i,2)=yy
            x(i,3)=zz
            xdot(i,1)=uu
            xdot(i,2)=vv
            xdot(i,3)=ww
100     continue

        end
        
        
        subroutine posvel(x,y,z,vx,vy,vz)

        implicit real*8 (a-h,o-z)

        parameter (nm=2500)
        common/profile/xx(0:nm),d(0:nm),v2(0:nm),psi(0:nm),
     &          zm(0:nm),indx(0:100)

        common/kpars1/w0,v0,e0,pi,twopi
        common/kpars2/nprof
        common/therm/g(0:1000)

        dimension v33(0:1000)
        real*8 random,r
        data v33(1000)/0./

        if (v33(1000).eq.0.) then
            do 5 i=0,1000
	       v33(i)=(.004*i)**3/3.
5           continue
        end if

c       Choose radius randomly from the mass distribution.

        r=random(0)
        i=100*r
        do 10 i1=indx(i),indx(i+1)+1
            if (zm(i1).gt.r) go to 15
10      continue

        write(6,*)' **** ERROR in king.'
        stop

15      rfac=(r-zm(i1-1))/(zm(i1)-zm(i1-1))
        r=xx(i1-1)+rfac*(xx(i1)-xx(i1-1))

c       Angular position random.

        cth=2.*random(0)-1.
        sth=sqrt(1.-cth**2)
        ph=twopi*random(0)
        cph=cos(ph)
        sph=sin(ph)
        x=r*sth*cph
        y=r*sth*sph
        z=r*cth

c       Choose speed randomly from the distribution at this radius.

        p=-(psi(i1-1)+rfac*(psi(i1)-psi(i1-1)))
        pfac=exp(p)
        il=0
        rl=0.
        iu=250.*sqrt(p)
        ru=pfac*g(iu)-v33(iu)
        r=random(0)*ru
100     if (iu-il.gt.1) then
            im=(il+iu)/2
            rm=pfac*g(im)-v33(im)
            if (rm.gt.r) then
                iu=im
                ru=rm
            else
                il=im
                rl=rm
            end if
            go to 100
        end if
        v=v0*(il+(r-rl)/(ru-rl))

c       Direction is random.

        cth=2.*random(0)-1.
        sth=sqrt(1.-cth**2)
        ph=twopi*random(0)
        cph=cos(ph)
        sph=sin(ph)
        vx=v*sth*cph
        vy=v*sth*sph
        vz=v*cth

        end

        
        subroutine poisson(x,d,v2,psi,zm,nmax,w0,nprof,v20,iout)

c       Self-contained 1-D (spherical) Poisson's equation solver.
c       Currently knows only about King models.
c
c       Input:	nmax is the maximum number of points allowed
c       	w0 is the dimensionless central potential
c       	iout allows messages if nonzero
c
c       Output:	x   is scaled radius
c       	d   is scaled density
c       	v2  is scaled velocity dispersion
c       	psi is potential
c       	zm  is cumulative mass (scaling from x, d scalings)
c       	nprof is the actual number of points generated
c       	v20 is the central velocity dispersion

        implicit real*8 (a-h,o-z)

        external rhs
        parameter (RLIN=.25, NLIN=100, RMAX=1.e4, TOL=1.e-6)

        dimension x(0:nmax),d(0:nmax),v2(0:nmax),psi(0:nmax),
     &            zm(0:nmax),y(2)

        common/center/di

        parameter (four3=1.333333333333333, fourpi=12.5663706)

	psi0=-abs(w0)

c       Initialize at center of cluster.

        xn=0.
        y(1)=0.
        y(2)=psi0
        x(0)=0.
        psi(0)=psi0
        v2(0)=1.
        call densvel(psi0,di,v2(0))
        d(0)=di
        di=1./di
        zm(0)=0.
        fac=10.**(log10(RMAX/RLIN)/(NMAX-NLIN))

c       Equation to be solved is
c
c		(x psi)''  = 9 x d
c
c       by defining y(1) = (x psi)
c       	    y(2) = y(1)'
c
c       Cover the first RLIN core radii linearly with NLIN points;
c       the remaining coverage is logarithmic, out to RMAX core radii,
c       if necessary.  Stop when d <= 0.

        do 1000 i=1,nmax
            xo=xn
            if (i.le.NLIN) then
                xn=(RLIN*i)/NLIN
            else
                xn=fac*xo
            end if
            call rkode(y,2,xo,xn,TOL,.1*(xn-xo),.0001*(xn-xo),
     &              nok,nbad,rhs,ier)
c           
c           N.B. Remember that y(1) is x*psi
c           
            if (ier.ne.0) then
                if (iout.ne.0)write(6,100)ier
100             format(1x,'ERROR #',i3,' IN RKODE')
                stop
            end if
            x(i)=xn
            psi(i)=y(1)/xn
            v2(i)=1.
            call densvel(psi(i),d(i),v2(i))
            if (d(i).lt.0.) then
                x(i)=x(i-1)+(x(i)-x(i-1))/(1.-d(i)/d(i-1))
                d(i)=0.
                v2(i)=0.
            end if
            zm(i)=x(i)*y(2)-y(1)
            if (d(i).le.0.) go to 2000
1000    continue
        i=nmax
2000    nprof=i
        v20=v2(0)
        do 2500 i=nprof,0,-1
            d(i)=d(i)/d(0)
            v2(i)=v2(i)/v2(0)
            zm(i)=fourpi/9.*zm(i)
2500    continue

99999   return
        end

c
        subroutine rhs(x,y,ypr)

c       Define RHS of ODE, for use by rkode.

        implicit real*8 (a-h,o-z)

        dimension y(2),ypr(2)
        common/center/di
        data zero/0./

        ypr(1)=y(2)

        if (x.le.0.) then
            d=1.
        else
            call densvel(y(1)/x,d,zero)
            d=d*di
        end if
        ypr(2)=9.*x*d

        end


        subroutine densvel(psi,d,v2)

c       Return scaled density d and velocity dispersion v2,
c       given scaled potential psi.

        implicit real*8 (a-h,o-z)

        d=0.
        if (psi.ge.0.) return

        dw=sqrt(-psi)
        w=dw
        g4=v2
        call gg(w,g2,g4)
        ep=exp(-psi)
        g2e=g2*ep
        d=g2e+dw*psi/3.
        if (v2.gt.0.and.d.gt.0.) v2=2.*(g4*ep-.2*dw*psi**2)/d

        end


        subroutine gg(w,g2,g4)

        implicit real*8 (a-h,o-z)

        common/therm/g(0:1000)

c       Array g contains the integral (0-->x) of y**2 exp(-y**2) dy

        external gaus2
        data eps/1.e-6/

        if (g(0).lt.0.) then
            g(0)=0.
            xo=0.
            dx=.004
            do 100 i=1,1000
                x=xo+dx
                g(i)=g(i-1)+qromb(gaus2,xo,x,eps,ntr)
                xo=x
100         continue
        end if

        ww=250.*w
        iw=ww
        if (iw.ge.1000) then
            g2=g(1000)
        else
            g2=g(iw)+(ww-iw)*(g(iw+1)-g(iw))
        end if
        if (g4.lt.0.) return
        w2=w*w
        ew=exp(-w2)
        wew=w*ew
        g4=1.5*g2-.5*w2*w*ew

        end


        function gaus2(x)

        implicit real*8 (a-h,o-z)

        gaus2=x*x*exp(-x*x)

        end
c
c=======================================================================

c       Routines from "Numerical Recipes".

c       	qromb - Romberg function integrator
c       	rkode - Runge-Kutta ODE integrator

        function qromb(f,a,b,eps,ntrap)

c       Integrate the function f from a to b using Romberg integration
c       with (relative) error tolerance eps.
c       |ntrap| on return is number of trapezoidal subdivisions.
c       ntrap < 0 indicates an error.

c       See "Numerical Recipes," Chapter 4.3


        implicit real*8 (a-h,o-z)

        parameter (one=1.,quarter=.25,zero=0.)
cd      parameter (one=1.d0,quarter=.25d0,zero=0.d0)
        parameter (jmax=20,jmaxp=jmax+1,k=5,km=k-1)
        dimension s(jmaxp),h(jmaxp)
        external f

        ntrap=0
        h(1)=one
        do 10 j=1,jmax
            call trapzd(f,a,b,s(j),j)
            ntrap=j
            if (j.ge.k) then
                call polint(h(j-km),s(j-km),k,zero,qromb,dss,ier)
                if (ier.ne.0) then
                    ntrap=-ntrap
                    qromb=0.
                    return
                end if
                if (abs(dss).lt.eps*abs(qromb)) return
            end if
            s(j+1)=s(j)
            h(j+1)=quarter*h(j)
10      continue
        ntrap=-1-ntrap

        end


        subroutine trapzd(f,a,b,s,n)

        implicit real*8 (a-h,o-z)

        external f
        parameter (half=.5)
        save it

        if (n.eq.1) then
            s=half*(b-a)*(f(a)+f(b))
            it=1
        else
            tnm=it
            del=(b-a)/tnm
            x=a+half*del
            sum=0.
            do 10 j=1,it
                sum=sum+f(x)
                x=x+del
10          continue
            s=half*(s+del*sum)
            it=2*it
        end if

        end


        subroutine polint(xa,ya,n,x,y,dy,ier)

        implicit real*8 (a-h,o-z)

        parameter (nmax=10)
        dimension xa(n),ya(n),c(nmax),d(nmax)

        ns=1
        dif=abs(x-xa(1))
        do 10 i=1,n
            dift=abs(x-xa(i))
            if (dift.lt.dif) then
                ns=i
                dif=dift
            end if
            c(i)=ya(i)
            d(i)=ya(i)
10      continue
        y=ya(ns)
        ns=ns-1
        do 20 m=1,n-1
            do 15 i=1,n-m
                ho=xa(i)-x
                hp=xa(i+m)-x
                w=c(i+1)-d(i)
                den=ho-hp
                if (den.eq.0.) then
                    ier=1
                    return
                end if
                den=w/den
                d(i)=hp*den
                c(i)=ho*den
15          continue
            if (2*ns.lt.n-m) then
                dy=c(ns+1)
            else
                dy=d(ns)
                ns=ns-1
            end if
            y=y+dy
20      continue
        ier=0

        end

c
c-----------------------------------------------------------------------

        subroutine rkode(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,
     &                   derivs,ier)

c       FIFTH-order Runge-Kutta integrator with adaptive step size.
c       See "Numerical Recipes," Chapter 15.2
c       -- Note that their common block /path/ is not used.

        implicit real*8 (a-h,o-z)

        parameter (maxstp=10000,nmax=10,two=1.,zero=0.,tiny=1.e-30)
        dimension ystart(nvar),yscal(nmax),y(nmax),dydx(nmax)
        external derivs

        ier=0
        x=x1
        h=sign(h1,x2-x1)
        nok=0
        nbad=0
        kount=0

        do 10 i=1,nvar
            y(i)=ystart(i)
10      continue

        do 50 nstp=1,maxstp
            call derivs(x,y,dydx)
            do 20 i=1,nvar
                yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny
20          continue
            if ((x+h-x2)*(x+h-x1).gt.zero)h=x2-x
            if (x+h.eq.x) go to 39
            call rk(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,ier)
            if (ier.ne.0) go to 9999
            if (hdid.eq.h) then
                nok=nok+1
            else
                nbad=nbad+1
            end if
            if ((x-x2)*(x2-x1).lt.zero) go to 45
39          do 40 i=1,nvar
                ystart(i)=y(i)
40          continue
            return
45          if (abs(hnext).lt.hmin) go to 99999
            h=hnext
50      continue

9999    ier=ier+1
99999   ier=ier+1

        end


        subroutine rk(y,dydx,n,x,htry,eps,yscal,
     &                hdid,hnext,derivs,ier)

        implicit real*8 (a-h,o-z)

        parameter (nmax=10,pgrow=-.2,pshrink=-.25,fcor=1./15.,
     &          one=1.,safety=.9,errcon=6.e-4)
        external derivs
        dimension y(n),dydx(n),yscal(n),ytemp(nmax),ysav(nmax),
     &            dysav(nmax)

        ier=0
        xsav=x
        do 10 i=1,n
            ysav(i)=y(i)
            dysav(i)=dydx(i)
10      continue
        h=htry
1       hh=.5*h
        call rk4(ysav,dysav,n,xsav,hh,ytemp,derivs)
        x=xsav+hh
        call derivs(x,ytemp,dydx)
        call rk4(ytemp,dydx,n,x,hh,y,derivs)
        x=xsav+h
        if (x.eq.xsav) go to 99999
        call rk4(ysav,dysav,n,xsav,h,ytemp,derivs)
        errmax=0.
        do 20 i=1,n
            ytemp(i)=y(i)-ytemp(i)
            errmax=max(errmax,abs(ytemp(i)/yscal(i)))
20      continue
        errmax=errmax/eps
        if (errmax.gt.one) then
            h=safety*h*(errmax**pshrink)
            go to 1
        else
            hdid=h
            if (errmax.gt.errcon) then
                hnext=safety*h*(errmax**pgrow)
            else
                hnext=4.*h
            end if
        end if
        do 30 i=1,n
            y(i)=y(i)+ytemp(i)*fcor
30      continue
        return

99999   ier=1

        end


        subroutine rk4(y,dydx,n,x,h,yout,derivs)

        implicit real*8 (a-h,o-z)

        parameter (nmax=10)
        dimension y(n),dydx(n),yout(n),yt(nmax),dyt(nmax),dym(nmax)
        external derivs

        hh=.5*h
        h6=h/6.
        xh=x+hh
        do 10 i=1,n
            yt(i)=y(i)+hh*dydx(i)
10      continue
        call derivs(xh,yt,dyt)
        do 20 i=1,n
            yt(i)=y(i)+hh*dyt(i)
20      continue
        call derivs(xh,yt,dym)
        do 30 i=1,n
            yt(i)=y(i)+h*dym(i)
            dym(i)=dyt(i)+dym(i)
30      continue
        call derivs(x+h,yt,dyt)
        do 40 i=1,n
            yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
40      continue

        end

c
c=======================================================================

c     Don't need to define random or anything below this line if a
c     built-in exists...

      function random(iseed)

c     Return a random number in the range [0,1), using the congruential
c     random number generator CONGRNO below.

c     Convention:  iseed = 0 ==> return next number in the sqruence.
c		   iseed > 0 ==> initialize the random generator.

      implicit real*8 (a-h, o-z)
      real*8 random
      real*8 congrno

      save irset
      data irset/0/

      if (irset.eq.0.or.iseed.gt.0) then

c         No system built-in random number generator is known:

          dum = congrno(iseed)

          irset = 1
      end if

c     Return the next random number in the sequence.

      random = congrno(0)

      end


      function congrno(iseed)

c     Random number generator (Press et al. p. 195).
c     ---------------------------------------------

      implicit real*8 (a-h,o-z)
      real*8 congrno

      parameter (M=714025, IA=1366, IC=150889, RM=1.d0/M)

      integer ir(97) 
      save iy,iff,ir,idum
      data iff /0/

      if (iseed.lt.0.or.iff.eq.0) then
          idum = -abs(iseed)
          iff = 1
          idum = mod(IC-idum,M)
          do 11 j = 1,97
              idum = mod(IA*idum+IC,M)
              ir(j) = idum
11        continue
          idum = mod(IA*idum+IC,M)
          iy = idum
      end if

      j = 1 + (97*iy)/m
      if (j.gt.97.or.j.lt.1) write (6,12) j, idum
12    format (/' RAN2:  j, idum = ',2i12)

      iy = ir(j)
      congrno = iy*RM
      idum = mod(IA*idum+IC,M)
      ir(j) = idum

      end 
