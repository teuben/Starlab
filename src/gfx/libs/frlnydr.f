
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
      
      subroutine fr lnydr(x,y1,y2,dys,dym,dyl,
     $                    iax,ilab,nlab,ndec,lpow)
      save
c     
c     As for frlnxdr, but for the y-axis
c     
c     Draw y-axis, tick marks and numbers for linear case (ndec nonzero),
c     and major tick marks and numbers for the logarithmic case (ndec = 0).
c     
c     iax	 =  1 for bottom axis, 2 for top.
c     ilab	 =  0 for no labels, 1 for labels
c     nlab	 =  estimated number of spaces for labels
c     
c     ndec = -1	==> integer format
c     ndec = 0	==> log axes, want to plot 10^n
c     ndec > 0:	number of places to right of decimal point
c     lpow = 0	==> F format, E format (lpow = exponent) otherwise
c     
c     No numbers are drawn for nonzero lmode, regardless of ilab.
c     N.B. Variable "ilab" is redundant, now.
c     
      character*80 device
      common/plot device/device,aspect,idev
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
      common/dev status/idevon,idevpen,idevwt
      common/fr hts/htl,htn/fr wts/iwts(4)/fr rotn/irot
      common/fr ticks/tiks(3),tikl/fr int/iframe
      common/fr conf/scent,rnuml,rnumr,snumt,snumb,dsnums,jrot,stopnum
      common/fr tik level/jtik level
      common/fr bare/ibare
      common/fontc1/dum2(5),rrmax,rrmin,ssmax,ssmin
      common/fr setax/kax,lax
c
      dimension dy(3)
      character*80 string
c     
      common/fr lby/lmode
      data mode/1/
c     
      cinch(xa,xs,di) = (xa-xs)*di
c     
      if (lax.lt.0) return
c     
c     Draw axes.
c     ---------
c     
      dy(1) = dys
      dy(2) = dym
      dy(3) = dyl
      rax = cinch(x,xl,dinchx)
      s1 = cinch(y1,ybot,dinchy)
      s2 = cinch(y2,ybot,dinchy)
      if (iwts(1).gt.0) then
          jwt = idevwt
          call weight(iwts(1))
      end if
      call plot(rax,s1,3)
      call plot(rax,s2,2)
c     
      if (lax.eq.1.and.iax.eq.2) return
      if (lax.eq.2.and.iax.eq.1) return
c     
c     Draw tick marks.
c     ---------------
c     
      dr = 1.
      if (iax.eq.2) dr = -1.
      do j = jtik level,3
          if (dy(j).ne.0.) then
              rtik = rax+dr*tiks(j)
              call fr lnfnc(y1,y2,dy(j),mode,firsty,ny)
              ds = dy(j)*dinchy
              s = cinch(firsty,ybot,dinchy)
              do i = 1,ny
                  call plot(rax,s,3)
                  call plot(rtik,s,2)
                  s = s+ds
              end do
          end if
      end do
      if (iwts(1).gt.0) call weight(jwt)
      jrot = irot
c     
      if (lmode.ne.0) return
      if (lax.eq.0.and.iax.eq.2) return
      if (lax.eq.3.and.iax.eq.1) return
c     
c     Add numeric labels.
c     ------------------
c
      yscale = max(abs(y1),abs(y2))
c
c     Labels are horizontal (left to right) if jrot = 0, horizontal otherwise.
c
      dsnums = ds
c     
c     Get the actual length of the numbers for use below.
c     
c     Definitions:	the maximum width  of the numeric labels is xlnum
c     the maximum height of the numeric labels is dsnum
c     
      if (ibare.eq.1) then
c         
c         A device-specific number-drawing routine will be used, so
c         pretty output is not expected.  Just estimate the space
c         requirements.
c         
          if (ndec.ne.0) then
              nnum = 0
              if (firsty.lt.0.) nnum = 1
              if (ndec.gt.0) nnum = nnum+ndec
              mm = 0
              do k = 1,ny
                  y = firsty+(k-1.)*dyl
                  if (abs(y).lt.1.) then
                      nn = 1
                  else
                      nn = log10(abs(y))+1
                  end if
                  if (nn.gt.mm) mm = nn
              end do
              nnum = nnum+mm
              xlnum = htn*(nnum+.5)
          else
              xlnum = 2.*htn
              mm = 0
              do k = 1,ny
                  y = firsty+(k-1.)*dyl
                  nn = 0
                  if (y.lt.0.) nn = 1
                  if (abs(y).ge.1.) nn = nn+1+log10(abs(y))
                  mm = max(mm,nn)
              end do
              xlnum = xlnum+.5333*mm*htn
          end if
      else
c         
c         Simbol will be used, so we can determine the sizes exactly.
c         
          xlnum = 0.
          do k=1,ny
              y = firsty + (k-1.)*dyl
              call make_number_string(y,ndec,lpow,string,nsym,yscale)
              call sim size(htn,string,nsym,dx,dum)
              xlnum = max(xlnum, dx)
          end do
          if (idev.eq.7.or.idev.eq.8) xlnum = 1.1*xlnum
c
c         Modify xlnum in case of extra labels in the logarithmic case:
c
          if (ndec.eq.0.and.nx.le.3) xlnum = 2.*xlnum
      end if
c
c     Nominal symbol height:
c
      dsnum = 1.5*htn
      if (jrot.ne.0) dsnum = 1.3*xlnum
c
c     Check that the number size is OK.
c
100   if (dsnum.gt.ds) then
c         
c         Label is too tall.
c         
          fac = ds/dsnum
          if (jrot.ne.0.and.fac.lt..4) then
              jrot = 0
              dsnum = 1.5*htn
              go to 100
          end if
          htn = htn*fac
          xlnum = xlnum*fac
          dsnum = ds
      end if
c     
      if (iax.eq.1) then
          call getlhe(space)
          space = -space - htl
      else
          call getrhe(space)
          space = space - x - htl
      end if
c     
      if (space.gt.0.and.jrot.eq.0.and.xlnum.gt..9*space) then
c         
c         Too wide...
c         
          htn = htn*.9*space/xlnum
          xlnum = .9*space
      end if
c
c     See if the first label should be moved up.
c     -----------------------------------------
c
c     (There is a lot of repeated code here--there is probably a much
c      better way of writing this.)
c
      y = firsty
      s1 = cinch(firsty,ybot,dinchy)
      half = .5*htn
      strlen = xlnum
      strlen2 = .5*strlen
c
      call pushstr
      if (jrot.eq.0) then
          if (iax.eq.1) then
              call strpos(1.,.5)
          else
              call strpos(0.,.5)
          end if
      else
          if (iax.eq.1) then
              call strpos(.5,0.)
          else
              call strpos(.5,1.)
          end if
      end if
c
      if (jrot.eq.0) then
          th = 0.
          if (iax.eq.1) then
              r = rax-half
          else
              r = rax+half
          end if
      else
          th = 90.
          if (iax.eq.1) then
              r = rax-htn
          else
              r = rax+htn
          end if
      end if
      s = s1
c
      smin = half			! Only offset the first y-axis
      if (th.eq.90.) smin = strlen2	! label if it really would
      if (s1.lt.smin) s = s+smin	! extend below the x-axis.
c
c     Reference point for first number is (r, s).
c
      if (iwts(2).gt.0) then
          jwt = idevwt
          call weight(iwts(2))
      end if
c     
      if (iax.eq.2) th = -th
c     
c     Draw the number and initialize bookkeeping.
c
      if (ibare.eq.1) then
          if (ndec.eq.0) then
              call lognum(r,s,htn,y,th,-1)
              rnummin = r
          else
              call fr numbr(r,s,htn,y,th,ndec)
              rnummin = rrmin
          end if
      else
          call make_number_string(y,ndec,lpow,string,nsym,yscale)
          call simbol(r,s,htn,string,th,nsym)
      end if
c
      rnummax = rnummin + xlnum
c     
      snumb = 0.
      snumt = slen
      rnuml = 0.
      rnumr = 0.
      inumset = 0
c
      if (jrot.eq.0) then
          sref = s - half
      else
          sref = s - strlen2
      end if
c     
c     sref is the level of the bottom of the current numerical label.
c
      dstrue = ds
      if (ndec.eq.0.and.ny.le.3.and.abs(dyl-1.).lt..01) then
c
c         In this case, we must take into account the fact that
c         intermediate numbers are being inserted, and that the bottom
c         label could be above the centerline...
c
          dstrue = .5*ds
          isign = 1
          if (sref.gt.scent) isign = -1
      end if
c
      if ( (sref-scent) * (sref+isign*dstrue-scent) .le. 0. ) then
c
c         This number and the next straddle the centerline.
c         The vertical whitespace for the label lies between snumb
c         and snumt.  The y-label will be centered there.
c
          inumset = 1
c
          snumt = sref+dstrue
c
          if (jrot.ne.0) then
              rnuml = r-htn
              rnumr = r
              snumb = s+strlen2
          else
              rnuml = r-strlen
              rnumr = r
              if (ndec.eq.0) then
                  snumb = sref+1.2*htn
              else
                  snumb = sref+htn
              end if
          end if
      end if
c
c     Draw the remaining numbers.
c     --------------------------
c
c     Undo any offset:
c
      if (s1.lt.smin) s = s-smin
c
      do i = 2,ny
          s = s+ds
          y = y+dyl
c         
          if (ibare.eq.1) then
              if (ndec.eq.0) then
                  call lognum(r,s,htn,y,th,-1)
                  rnummin = min(rnummin,r)
                  rnummax = max(rnummax,r+strlen)
              else
                  call fr numbr(r,s,htn,y,th,ndec)
                  rnummin = min(rnummin,rrmin)
                  rnummax = max(rnummax,rrmax)
              end if
          else
              call make_number_string(y,ndec,lpow,string,nsym,yscale)
              call simbol(r,s,htn,string,th,nsym)
c
c             Still have to deal with intermediates...
c
              if (ndec.eq.0.and.ny.le.3.and.abs(dyl-1.).lt..01) then
                  s3 = s-.5229*ds
                  if (s3.ge.s1) then
                      call format_string(10.**(y-.5*dyl),-1,
     $                                   string,nsym,1)
                      call simbol(r,s3,htn,string,th,nsym)
                  end if
c
                  if (i.eq.2) then
                      s3 = s-1.5229*ds
                      if (s3.ge.s1) then
                          call format_string(10.**(y-1.5*dyl),-1,
     $                                       string,nsym,1)
                          call simbol(r,s3,htn,string,th,nsym)
                      end if
                  end if
c                 
                  if (i.eq.ny) then
                      s3 = s+.4771*ds
                      if (s3.le.s2) then
                          call format_string(10.**(y+.5*dyl),-1,
     $                            string,nsym,1)
                          call simbol(r,s3,htn,string,th,nsym)
                      end if
                  end if
              end if
c
              rnummin = min(rnummin,rrmin)
              rnummax = max(rnummax,rrmax)
          end if
c         
c         Bookeeping (same as before):
c
          if (inumset.eq.0) then
              if (jrot.eq.0) then
                  sref = s - half
              else
                  sref = s - strlen2
              end if
c
              if (sref.le.scent.and.sref+dstrue.gt.scent) then
                  inumset = 1
c
                  snumt = sref+dstrue
c
                  if (jrot.ne.0) then
                      rnuml = r-htn
                      rnumr = r
                      snumb = s+strlen2
                  else
                      rnuml = r-strlen
                      rnumr = r
                      snumb = sref+htn
c
c                     Old version:
c
c                     if (ndec.eq.0) then
c                         snumb = sref+1.2*htn
c                     else
c                         snumb = sref+htn
c                     end if
c
                  end if
              end if
          end if
      end do
c
c     Restore settings.
c
      if (iwts(2).gt.0) call weight(jwt)
      iframe = 1
      call popstr
c     
      rnuml = rnummin
      rnumr = max(rnummax,rnuml+htn*jrot+strlen*(1-jrot))
c     
      snumb = max(0.,snumb)
      snumt = min(slen,snumt)
      stopnum = s
c
c     Force the label to be above the centerline.
c
      if (.5*(snumb+snumt).lt.scent) then
          snumb = snumb + dstrue
          snumt = snumt + dstrue
      end if
c     
      end


      subroutine make_number_string(y,ndec,lpow,string,nsym,yscale)
      character*(*) string
      common /fr plain/ iplain
c
      parameter (TOL = 1.e-6)
c     
      if (lpow.eq.0.or.(ndec.ne.0.and.abs(y)/yscale.lt.TOL)) then
          call numsym(y,ndec,string,nsym)
          if (iplain.eq.1) call convert_to_plain(string, nsym)
      else
          if (ndec.eq.0) then
              call format_string(10.**y,ndec,string,nsym,0)
          else
              call exp_string(y,ndec,string,nsym,1)
          end if
      end if
c
      end
