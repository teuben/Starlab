      program plot_SeBa

c general plotting program for reduced SeBa output

c      parameter(Nmax=36389,Np=30000)
      parameter(Nmax=100000,Nr=20,Ni=6)
      
      real all_r(Nmax,Nr),arr1(Nmax),arr2(Nmax)
      integer all_i(Nmax,Ni)

      real rmin(Nr),rmax(Nr),xmin,xmax,ymin,ymax
      integer nbin,istat1,istat2
      character*30 label(Nr),xlabel,ylabel,filename,filename2
      character*100 title
      integer pgopen

      data rmin/0.1,1.,1.,0.,0.,0.8,0.5,0.8,0.5,0.1,0.1,
     +     1.e-4,0.,0.,0.1,0.05,0.1,0.05,1.,1./
      data rmax/0.1,1.e4,1.e3,1.,1.,10.,10.,10.,10.,1.e4,1.e2,
     +     1.e2,1.,1.,10.,10.,10.,10.,20.,20./

      data label/'T\di\u (Myr)','a\di\u (R\d\(2281)\u)','P\di\u (d)',
     +     'q\di\u',
     +     'ecc\di\u','M\di\u (M\d\(2281)\u)','R\di\u (R\d\(2281)\u)', 
     +     'm\di\u (M\d\(2281)\u)','r\di\u (R\d\(2281)\u)',
     +     'T\df\u (Myr)','a\df\u (R\d\(2281)\u)','P\df\u (d)',
     +     'q\df\u',
     +     'ecc\df\u','M\df\u (M\d\(2281)\u)','R\df\u (R\d\(2281)\u)', 
     +     'm\df\u (M\d\(2281)\u)','r\df\u (R\d\(2281)\u)',
     +     'id1','id2'/

      title = ' '
      nbin=30


      write(6,*) 'give filename'
      read(5,"(A)") filename

      write(6,*) filename
      n=0
      call read_in(all_r,all_i,filename,n,Nmax,Nr,Ni)  

      write(6,*) n
      do i=1,n
c Pi (in days), qi
         all_r(i,3) = 0.116*sqrt(all_r(i,2)**3/(all_r(i,6)+all_r(i,8)))
         all_r(i,4) = all_r(i,8)/all_r(i,6)
c Pf (in days), qf
         all_r(i,12) = 0.116*sqrt(all_r(i,11)**3
     &               / (all_r(i,15)+all_r(i,17)))
         if (all_r(i,15).eq.0.or.all_r(i,17).eq.0) then
            all_r(i,13) = 0
            all_r(i,17) = max(all_r(i,15),all_r(i,17))
            all_r(i,15) = 0
         else
            all_r(i,13) = all_r(i,17)/all_r(i,15)
         end if 
         all_r(i,19) = all_i(i,5)
         all_r(i,20) = all_i(i,6)
      end do
      
c      istat1=pgopen('?')
      call pgbeg(1,'?',2,2)
      istat1 =1
      print*, istat1
c      call pgsubp(2,2)
      call pgslw(3)
      call pgsch(1.5)

 10   write(6,*) 'Choose which parameters to plot'
      write(6,*) ' 1 Ti        2 ai         3 Pi'
      write(6,*) ' 4 qi        5 ecci       6 Mi'     
      write(6,*) ' 7 Ri        8 mi         9 ri'
      write(6,*) '10 Tf       11 af        12 Pf'
      write(6,*) '13 qf       14 eccf      15 Mf'     
      write(6,*) '16 Rf       17 mf        18 rf'
      write(6,*) '19 id1      20 id2'

      read(5,*) inum1,inum2

      if (inum1.eq.inum2) then         
         do i=1,n
            if (inum1.gt.0) then
               arr1(i)=all_r(i,inum1)
               xmin = rmin(inum1)
               xmax = rmax(inum1)
               xlabel = label(inum1)
            else
               arr1(i)=log10(all_r(i,-inum1))
               xmin = log10(rmin(-inum1))
               xmax = log10(rmax(-inum1))
               xlabel = 'log '//label(-inum1)
            end if
         end do
 20      call pghist(n,arr1,xmin,xmax,nbin,0)
         call pglabel(xlabel,'N',title)
         write(6,*) 'change limits?'
         read(5,"(A)") ch         
         if (ch.eq.'y'.or.ch.eq.'Y') then
            write(6,*) 'give new limits'
            read(5,*) xmin,xmax
            rmin(-inum1) = xmin
            rmax(-inum1) = xmax
            goto 20
         end if
      else
         do i=1,n
            if (inum1.gt.0) then
               arr1(i)=all_r(i,inum1)
               xmin = rmin(inum1)
               xmax = rmax(inum1)
               xlabel = label(inum1)
            else
               arr1(i)=log10(all_r(i,-inum1))
               xmin = log10(rmin(-inum1))
               xmax = log10(rmax(-inum1))
               xlabel = 'log '//label(-inum1)
            end if               
            if (inum2.gt.0) then
               arr2(i)=all_r(i,inum2)
               ymin = rmin(inum2)
               ymax = rmax(inum2)
               ylabel = label(inum2)
            else
               arr2(i)=log10(all_r(i,-inum2))
               ymin = log10(rmin(-inum2))
               ymax = log10(rmax(-inum2))
               ylabel = 'log '//label(-inum2)
            end if
         end do
 30      call spzgray(n,arr1,xmin,xmax,arr2,ymin,
     +        ymax,nbin,nbin,0.,1.,0)
         call pglabel(xlabel,ylabel,title)
         write(6,*) 'change limits?'
         read(5,"(A)") ch         
         if (ch.eq.'y'.or.ch.eq.'Y') then
            write(6,*) 'give new x-limits'
            read(5,*) xmin,xmax
            if (inum1.gt.0) then
               rmin(inum1) = xmin 
               rmax(inum1) = xmax 
            else
               rmin(-inum1) = 10**xmin
               rmax(-inum1) = 10**xmax
            end if               
            write(6,*) 'give new y-limits'
            read(5,*) ymin,ymax
            if (inum2.gt.0) then
               rmin(inum2) = ymin 
               rmax(inum2) = ymax 
            else
               rmin(-inum2) = 10**ymin
               rmax(-inum2) = 10**ymax
            end if               
            goto 30
         end if
      end if
      
      write(6,*) 'Stop? (s) or Print (p)'
      read(5,"(A)") ch

      if (ch.ne.'s'.and.ch.ne.'p') goto 10      

      if (ch.eq.'p') then
         write(6,*) 'give filename'
         read(5,*) filename
c         istat2 = pgopen(filename//'/ps')
         write(6,*) 'give title (optional)'
         read(5,"(A)") title
         
         if (inum1.eq.inum2) then         
            call pghist(n,arr1,xmin,xmax,nbin,0)
            call pglabel(xlabel,'N',title)
         else
            call spzgray(n,arr1,xmin,xmax,arr2,ymin,
     +           ymax,nbin,nbin,0.,1.,0)
            call pglabel(xlabel,ylabel,title)
         end if
c         call pgclos
         print*, istat1, istat2
c         call pgslct(istat1)
         goto 10
      end if

c      call pgclos
c      call end
      end

      SUBROUTINE read_in(all_r, all_i, filename, n, Nmax, Nr, Ni)

      real all_r(Nmax,Nr)
      integer all_i(Nmax,Ni)
      character*30 filename

      open(10,file=filename,status="old")

      n=0
      do i=1,Nmax
         read(10,*,end=99) all_i(i,1),all_r(i,1),all_r(i,2),
     +        all_r(i,5),all_i(i,2),all_r(i,6),all_r(i,7),
     +        all_i(i,3),all_r(i,8),all_r(i,9)
         read(10,*,end=99) all_i(i,4),all_r(i,10),all_r(i,11),
     +        all_r(i,14),all_i(i,5),all_r(i,15),all_r(i,16),
     +        all_i(i,6),all_r(i,17),all_r(i,18)
         n=n+1
      end do

 99   close(1)
      return

      end

      subroutine spzgray(ndata, xdata, xmin, xmax, ydata, ymin, ymax
	1    , nxbin,nybin, iflag)
	integer ndata, nxbin, nybin, iflag
	real xdata, ydata, xy_data
	real xmin, xmax, ymin, ymax 
	dimension xdata(ndata), ydata(ndata), xy_data(nxbin, nybin)
	real dx, dy, tr(6), fg, bg, norm     
	integer nxmin, nxmax, ix, iy, it
	norm=0
	nxmin=1
	nymin=1
	dx = (xmax-xmin)/(1.0*nxbin)
	dy = (ymax-ymin)/(1.0*nybin)
	bg = 0.0
	fg = 1.0
	tr(1) = xmin
	tr(2) = dx
	tr(3) = 0
	tr(4) = ymin
	tr(5) = 0
	tr(6) = dy
	print *, ndata, xmin, xmax, ymin, ymax, nxbin,nybin, iflag
	print *, dx, dy

	do ix=1, nxbin
	   do iy=1, nybin
	      xy_data(ix, iy) = 0
	   enddo
	enddo
	do i=1, ndata
c	   print *, i, xdata(i), ydata(i)
	   ix = 0.5 + (xdata(i)-xmin)/dx
	   iy = 0.5 + (ydata(i)-ymin)/dy
c	   print *, ix, iy, xdata(i), ydata(i)
	   if(ix.ge.1.and.iy.ge.1.and.ix.le.nxbin.and.iy.le.nybin)then
              xy_data(ix, iy) = xy_data(ix, iy) + 1
              norm = max(norm, xy_data(ix, iy))
c	      print *, ix, it, xy_data(ix, it), norm
	   endif
	enddo
	do ix=1, nxbin
	   do iy=1, nybin
	      xy_data(ix, iy) = xy_data(ix, iy)/norm
	   enddo
	enddo

	if(iflag.eq.0)then
	   call pgenv(xmin,xmax, ymin,ymax,0,0)
	endif
	call pggray(xy_data, nxbin, nybin, nxmin, nxbin, nymin,
	1    nybin, fg, bg, tr)
	end


      
