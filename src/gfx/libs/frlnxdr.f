
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
      subroutine fr lnxdr(y,x1,x2,dxs,dxm,dxl,
     $                    iax,ilab,nlab,ndec,lpow)
      save
c     
c     Draw x-axis, tick marks and numbers for linear case (ndec nonzero),
c     and major tick marks and numbers for the logarithmic case (ndec = 0).
c
c     iax	= 1 for bottom axis, 2 for top.
c     ilab	= 0 for no labels, 1 for labels
c     nlab	= estimated number of spaces for labels
c
c     ndec = -1	==> integer format
c     ndec = 0	==> log axes, want to plot 10^n
c     ndec > 0:	number of places to right of decimal point
c     lpow = 0	==> F format, E format (lpow = exponent) otherwise
c     
c     No numbers are drawn for nonzero lmode, regardless of ilab.
c     
      dimension dx(3)
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
      common/dev status/idevon,idevpen,idevwt
      common/fr hts/htl,htn/fr wts/iwts(4)
      common/fr ticks/tiks(3),tikl
      common/fr xnums/xnumbot
c     
      common/fr bare/ibare
      common/fr lbx/lmode
      common/fr tik level/jtik level
c
      character*80 string
      parameter (TOL = 1.e-6, WTOL = 0.75)
c
      data mode/1/jtik level/1/
c     
      data theta/0./
c     
      cinch(xa,xs,di)=(xa-xs)*di
c     
      dx(1)=dxs
      dx(2)=dxm
      dx(3)=dxl
      sax=cinch(y,ybot,dinchy)
      r1=cinch(x1,xl,dinchx)
      r2=cinch(x2,xl,dinchx)
c     
      if (iwts(1).gt.0) then
          jwt=idevwt
          call weight(iwts(1))
      end if
c     
      call plot(r1,sax,3)
      call plot(r2,sax,2)
c     
      ds=1.
      if (iax.eq.2) ds=-1.
c     
c     Draw the tick marks: j = 1, 2, 3 correspond to small, medium and
c     large tick marks, respectively.  Tick sizes are set in eframe.
c     
      do j=jtik level,3
          if (dx(j).ne.0) then
              stik=sax+ds*tiks(j)
c
              call fr lnfnc(x1,x2,dx(j),mode,firstx,nx)
c
              dr=dx(j)*dinchx
              r=cinch(firstx,xl,dinchx)
              do i=1,nx
                  call plot(r,sax,3)
                  call plot(r,stik,2)
                  r=r+dr
              end do
          end if
      end do
c     
      if (iwts(1).gt.0) call weight(jwt)
      if (ilab.eq.0.or.lmode.ne.0) return
c     
c     Now add the numbers.
c     -------------------
c
c     Numbers will be placed at intervals of dx1, starting at firstx.
c     Offset is above or below the axis, depending on iax.
c     Amount of offset scales with htn.
c
c     Determine placement above or below the axis...
c     
      htnsave = htn
c
      x=firstx
      r=cinch(firstx,xl,dinchx)
c
c     First determine the width of a "typical" number string...
c
      if (ibare.eq.1.or.(ndec.ne.0.and.lpow.eq.0)) then
          if (ibare.eq.1) then
c
c             Just use the nominal width of the number string.
c
              wid=.85*htn*(nlab+.5)
          else
              call numsym(x1,ndec,string,nsym)
              call sim size(htn,string,nsym,wid,dum)
              call numsym(x2,ndec,string,nsym)
              call sim size(htn,string,nsym,wid2,dum)
              if (wid2.gt.wid) wid = wid2
          end if
c
      else
          if (ndec.eq.0) then
c
c             Nominal dimensions:
c
c             wid = .85*2.5*htn
c
c             Better:
c
              call sim size(htn,'10^+n',5,wid,dum)
c
          else
              call exp_string(x1,ndec,string,nsym,1)
              call sim size(htn,string,nsym,wid,dum)
              call exp_string(x2,ndec,string,nsym,1)
              call sim size(htn,string,nsym,wid2,dum)
              if (wid2.gt.wid) wid = wid2
          end if
      end if
c     
c     Reduce the scale if the estimated size is too great.
c     
      if (wid.gt.WTOL*dr) htn=htn*WTOL*dr/wid
c
c     ...then determine the vertical offset...
c
      if (ibare.eq.1.or.(ndec.ne.0.and.lpow.eq.0)) then
          if (iax.eq.1) then
              s=sax-1.5*htn
          else
              s=sax+.5*htn
          end if
      else
          if (iax.eq.1) then
              s=sax-1.75*htn
          else
              s=sax+.5*htn
          end if
      end if
c     
c     ...and update the data on the bottom of the number field.
c     
      if (iax.eq.1)then
          xnumbot=s
          if (ndec.ne.0) xnumbot=xnumbot-.5*htn
      end if
c     
      if (iwts(2).gt.0)then
          jwt=idevwt
          call weight(iwts(2))
      end if
c     
c     Plot the numbers.
c     
      call pushstr
      call strpos(.5,0.)
c
      rleft=cinch(x1,xl,dinchx)
      rright=cinch(x2,xl,dinchx)
      if (rleft.gt.rright) then
          temp = rleft
          rleft = rright
          rright = temp
      end if
c
      xscale = max(abs(x1),abs(x2))
      do i=1,nx
          if (lpow.eq.0.or.ibare.eq.1
     $            .or.(ndec.ne.0.and.abs(x)/xscale .lt. TOL)) then
              call fr numbr(r,s,htn,x,theta,ndec)
          else
              if (ndec.eq.0) then
c                 
c                 Logarithmic plot.
c                 
c                 First see if we need any intermediate labels.
c                 Note that there will still be cases where no labels
c                 appear (narrow ranges in logarithmic plots...).
c
                  if (nx.le.3.and.abs(dxl-1.).lt..01) then
                      r3 = r-.5229*dr
                      if (r3.ge.rleft) then
                          call format_string(10.**(x-.5*dxl),-1,
     $                                       string,nsym,1)
                          call simbol(r3,s,htn,string,theta,nsym)
                      end if
                      if (i.eq.nx) then
                          r3 = r+.4771*dr
                          if (r3.le.rright) then
                              call format_string(10.**(x+.5*dxl),-1,
     $                                           string,nsym,1)
                              call simbol(r3,s,htn,string,theta,nsym)
                          end if
                      end if
                  end if
c
                  call format_string(10.**x,ndec,string,nsym,0)
                  call simbol(r,s,htn,string,theta,nsym)
c
              else
c
c                 Exponential format on a linear plot.
c
                  call exp_string(x,ndec,string,nsym,1)
              end if
c
              call simbol(r,s,htn,string,theta,nsym)
c
          end if
          r=r+dr
          x=x+dxl
      end do
c
      call popstr
c                 
      if (iwts(2).gt.0) call weight(jwt)
c
c     DON'T restore the height here if we want changes to propogate
c     to the y-axis labels.
c
      call getyfollowsx(iy)
      if (iy.eq.0) htn = htnsave
c     
      end


      subroutine setyfollowsx(iy)
c
c     If ixy = 0, then the y number sizes will be INDEPENDENT of
c     changes made by frlnxdr.  If ixy = 1, changes in the x-axis
c     number sizes will propogate to the y axis.
c
      save
      common /yfollowsx/ixy
      data ixy/1/
c
      ixy = iy
      if (ixy.ne.0) ixy = 1
      return
c
      entry getyfollowsx(iy)
      iy = ixy
      end


      subroutine format_string(x,ndec,string,nsym,icoef)
      character*(*) string
      parameter (TOL = 0.1)
c
      common /fr plain/ iplain
c
c     Call exp_string, but check for special cases first.
c
      isp = 1
      if (abs(x/.01-1.).lt.TOL) then
          string = '0.01'
          nsym = 4
      else if (abs(x/.03-1.).lt.TOL) then
          string = '0.03'
          nsym = 4
      else if (abs(x/.1-1.).lt.TOL) then
          string = '0.1'
          nsym = 3
      else if (abs(x/.3-1.).lt.TOL) then
          string = '0.3'
          nsym = 3
      else if (abs(x/1.-1.).lt.TOL) then
          string = '1'
          nsym = 1
      else if (abs(x/3.-1.).lt.TOL) then
          string = '3'
          nsym = 1
      else if (abs(x/10.-1.).lt.TOL) then
          string = '10'
          nsym = 2
      else if (abs(x/30.-1.).lt.TOL) then
          string = '30'
          nsym = 2
      else if (abs(x/100.-1.).lt.TOL) then
          string = '100'
          nsym = 3
      else if (abs(x/300-1.).lt.TOL) then
          string = '300'
          nsym = 3
      else
          isp = 0
          call exp_string(x,ndec,string,nsym,icoef)
      end if
c
      if (isp.eq.1.and.iplain.eq.1)
     $        call convert_to_plain(string, nsym)
c
      end


      subroutine exp_string(x,ndec,string,nsym,icoef)
      character*(*)string
      character*10 is
c
      common /fr plain/ iplain
c
c     Convert x into a "simbol" string in exponential format.
c     Skip the leading exponent if icoef = 0.
c
      call compoz(x,f,n)
c
      if (icoef.eq.0) then
          string(1:2) = '10'
          nsym = 2
      else
c
c         We probably don't want exponential format if n = -1, 0, or 1.
c
          if (abs(n).le.1) then
c
c             May need to add or remove decimal digits to maintain
c             precision in this case.
c
              call numsym(x,ndec-n,string,nsym)
          else
              call numsym(f,ndec,string,nsym)
              nt = 11
              string(nsym+1:nsym+nt) = '%@%@ %@@*10'
              nsym = nsym + nt
          end if
      end if
c
      if (icoef.eq.0.or.abs(n).gt.1) then
          write(is,'(i10)')n
          do i=1,10
              if (is(i:i).gt.' ') then
                  string(nsym+1:nsym+1) = '^'
                  string(nsym+2:nsym+2) = is(i:i)
                  nsym = nsym + 2
              end if
          end do
      end if
c
      if (iplain.eq.1) call convert_to_plain(string, nsym)
c
c     write(6,*)'x, str = ',x,'  ',string(1:nsym)
c
      end


      subroutine convert_to_plain(string, nsym)
      character*(*) string
      character*80 strtmp
c
c     Make digits plain font.
c
      ns = nsym
      strtmp = string
      nsym = 0
      do i=1,ns
          if (strtmp(i:i).lt.'0'.or.strtmp(i:i).gt.'9') then
              nsym = nsym + 1
              string(nsym:nsym) = strtmp(i:i)
          else
              nsym = nsym + 1
              string(nsym:nsym) = '@'
              nsym = nsym + 1
              string(nsym:nsym) = strtmp(i:i)
          end if
      end do
c
      end
