
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

      subroutine eframe(xmin,xmax,xlen,modx,xctit,
     *                  ymin,ymax,ylen,mody,yctit)
        save
c
      character*(*) xctit,yctit
      character*100 outbuf
c
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
      common/dev status/idevon,idevpen,idevwt
      common/frdraw/mode/frhts/htl,htn/frwts/iwts(4)
      common/frpens/icolors(3)/frrotn/irot
      common/frticks/tiks(3),tikl/frint/iframe
      common/frconf/scent,rnuml,rnumr,snumt,snumb,dsnums,jrot,stopnum
      common/frbare/ibare
      common/frsetax/kax,lax
      common/dev init/init dev
      common/debug trace/itrace
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c      Draw frame for plot with tick marks, numerical labels, and
c      titles, using the extended (SIMBOL) font set.
c
c
c      input: ( y is similar to x)
c      ------
c
c       xmin plot value at left hand side
c       xmax plot value at right hand side
c       xlen length of x axis in inches
c       modx = 1 linear plot limits correspond to xmin xmax
c       modx = 2 linear plot limits adjusted to contain to xmin xmax
c       modx =-1 log plot limits correspond to xmin xmax
c       modx =-2 log plot limits adjusted to contain to xmin xmax
c       xctit (character) contains x title
c      *** for log plot enter actual variable, i.e. .01 not -2
c
c
c      output:
c      -------
c
c       xl actual value of left limit
c       xr actual value of right limit
c       dinchx inches per plot unit, i.e. xlen/(xr-xl)
c       ybot actual value of bottom limit
c       ytop actual value of top limit
c       dinchy inches per plot unit, i.e. ylen/(ytop-ybot)
c      *** for log plots limit is log10 of variable
c
c
c	switches:
c       ---------
c
c	the following common blocks each may contain one integer*4
c	variable (imode, say), whose effect is as described.
c
c	(i)  /frdraw/  if imode is nonzero, only scaling information is
c 		        returned -- nothing is drawn,
c	(ii) /frbnds/  for nonzero imode, only that part of the
c		        graph (produced by m(d)line) lying within the
c		        "frame"-defined box is actually plotted,
c	(iii)/frlbx/  the x-axis is numbered only for imode=0,
c	(iv) /frlby/  the y-axis is numbered only for imode=0.
c	(v)  /frrotn/	an attempt will be made to keep all y-axis
c			labels horizontal if imode is zero. numerical labels
c			longer than six characters and text labels with
c			length greater than  max( 1.2, 7.5*htl )
c			will still be plotted vertically.
c	(vi) /frplain/ if imode is nonzero, eframe will use "nombr" for the
c			numbers, to save on time.
c       (vii) /frbare/ if imode is nonzero, no labels will be drawn and
c			"numbr" will be used for the numbers.
c			
c
c	The switches in (i) to (v) above may be set with
c	"call setmod(im1,im2,im3,im4,im5)".
c	/frplain/ is set using subroutine setsym.
c
c	Thus, if no switch is set, entire graphs will be drawn and both
c	axes will be numbered (with horizontal labels, if possible).
c
c
c       Other variable parameters:
c	--------------------------
c
c	(vi) htl, htn, in common block /frhts/, give the sizes of the
c	titles and numerical labels, respectively. defaults are .15, .15.
c	The sizes of all tick marks along the axes scale with htn.
c
c	Set heights with "call sethts(ht1,ht2)".
c
c	(vii) the weights of various components of the frame may be set
c	individually via the integer*4 array iwts in common /frwts/:
c
c		iwts(1):   box and tick marks.
c		iwts(2):   numerical labels (excluding exponents, if any).
c		iwts(3):   exponents (default = iwts(2)).
c		iwts(4):   text labels.
c
c	s/r weight is called with argument iwts(.), when necessary.
c	defaults are 0,0,0,0.
c
c	Set weights with "call setwts(iw1,iw2,iw3,iw4)".
c
c	(viii) the pen types (colors) of various frame components may be set
c	individually via the integer*4 array icolors in common /frpens/:
c
c		icolors(1):   box and tick marks.
c		icolors(2):   numerical labels (including exponents, if any).
c		icolors(3):   text labels.
c
c	s/r color is called with argument icolors(.), when necessary.
c	defaults are 0,0,0.
c
c	Set pens with "call setpens(ip1,ip2,ip3)".
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      data kountmax/5/
c
      if(itrace.eq.1)open(2,file='UNIT2')
c
      if(xmax.eq.xmin.or.ymax.eq.ymin.or.xlen.eq.0..or.ylen.eq.0.
     +   .or.(modx.lt.0.and.(xmin.le.0..or.xmax.le.0.))
     +   .or.(mody.lt.0.and.(ymin.le.0..or.ymax.le.0.)))then
          call display text('error in eframe arguments:',26)
          write(outbuf,11111)'x',xmin,xmax,xlen,modx
11111     format(a1,': ',1p3e15.6,i10)
          call display text(outbuf,58)
          write(outbuf,11111)'y',ymin,ymax,ylen,mody
          call display text(outbuf,58)
          xl=0.
          xr=0.
          dinchx=0.
          ybot=0.
          ytop=0.
          dinchy=0.
          rlen=0.
          slen=0.
          return
      end if
c
      if(init dev.eq.0)then
          init dev=-1
          call mcinit
          call devon
          call clear
      end if
      if(idevon.eq.0)call devon
c
      call routine id('eframe')
      iframe=1
c
c     Experimental precautionary measure:
c
      call clrstr
c
c-----------------------------------------------------------------------------
c
c     Adjust label and tick sizes, if the plot is small.
c
      dmin=min(xlen,ylen)
      if(dmin.gt.2.)then
          if(htl.eq.0.)htl=.15
          if(htn.eq.0.)htn=.15
      else
          if(htl.eq.0.)htl=.075*dmin
          if(htn.eq.0.)htn=htl
      end if
c
      tiks(1)=.2*htn
      tiks(2)=.30*htn
      tiks(3)=.45*htn
      tikl=.2*htn
c
c     Save current heights, pen and weight settings, for restoration at the end.
c
      htlsto=htl
      htnsto=htn
      call getstatus(idum1,idum2,icolorinit,iwtinit)
      irotsto=irot
c
c     These will inform a new X-window of the curent colors...
c
      call color(icolorinit)
      call weight(iwtinit)
c
c     Common variables:
c     ----------------
c
      modex=modx
      modey=mody
      rlen=xlen
      slen=ylen
c
c     Initialize parameters set by the fr*dr routines.
c
      iconf=0
      scent=.5*slen
c
      if (icolors(1).gt.0) call color(icolors(1))
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     Set axis, tick, and label info, then draw the axes and labels.
c
c     The x-axis action is controlled by x1, x2, y1, y2, modex, mode, kax.
c     The y-axis action is controlled by x1, x2, y1, y2, modey, mode, lax.
c
c     NOTE that both axes must be "set" first, and the fr*set routines
c     may modify x1, x2, y1, and y2.
c
c     NOTE also that the fr*dr routines return global variables used
c     by the "labels" routine.
c
c     Determine appropriate x tick spacings...
c
      x1=xmin
      x2=xmax
      if (modex.gt.0) then
          call frlnset(x1,x2,xlen,modex,dxs,dxm,dxl,
     &			labxsp,labxdp,lpowx)
c
c         On return from frlnset, 
c
c         labxsp      = total number of spaces for label
c         labxdp=-1   : integer format
c         labxdp.gt.0 : number of places to right of decimal point
c         lpowx.ne.0  : E format, otherwise F format
c
c         (Same for y below...)
c
      else
          call frlgset(x1,x2,xlen,modex,dxs,dxm,dxl)
      end if
c
c     ...and set the x scaling in common /scales/.
c
      xl=x1
      xr=x2
      dinchx=xlen/(xr-xl)
c
c     Do the same for y.
c
      y1=ymin
      y2=ymax
      if (modey.gt.0) then
          call frlnset(y1,y2,ylen,modey,dys,dym,dyl,
     &			labysp,labydp,lpowy)
      else
          call frlgset(y1,y2,ylen,modey,dys,dym,dyl)
      end if
c
      ybot=y1
      ytop=y2
      dinchy=ylen/(ytop-ybot)
c
c     Draw axes, tick marks, labels.
c     -----------------------------
c

c        write(6,*)'iwts = ',iwts
c        write(6,*)'icolors = ',icolors
c        write(6,*)'mode, kax = ',mode,kax

      if (mode.eq.0.and.kax.gt.0) then
          if (modex.gt.0) then
              call frlnxdr(y1,x1,x2,dxs,dxm,dxl,
     $                      1,1,labxsp,labxdp,lpowx)
              call frlnxdr(y2,x1,x2,dxs,dxm,dxl,
     $                      2,0,labxsp,labxdp,lpowx)
          else
              call frlgxdr(y1,x1,x2,dxs,dxm,dxl,1,1)
              call frlgxdr(y2,x1,x2,dxs,dxm,dxl,2,0)
	  end if
          call labl(1,xctit)
      end if
c
      if (mode.eq.0.and.lax.ge.0) then
c
c         Numbers are always drawn horizontally.
c
          irot=0
c
          if (modey.gt.0) then
     	      call frlnydr(x1,y1,y2,dys,dym,dyl,
     $                     1,1,labysp,labydp,lpowy)
     	      call frlnydr(x2,y1,y2,dys,dym,dyl,
     $                     2,0,labysp,labydp,lpowy)
          else
              call frlgydr(x1,y1,y2,dys,dym,dyl,1,1)
              call frlgydr(x2,y1,y2,dys,dym,dyl,2,0)
	  end if
	  irot = irotsto
          call labl(2,yctit)
      end if
c
c-----------------------------------------------------------------------------
c
c     Restore "true" settings (just in case).
c
      htl = htlsto
      htn = htnsto
      call color(icolorinit)
      call weight(iwtinit)
      irot = irotsto
c
      iframe=0
c
      end


      subroutine labl(which,string)
        save
      integer which
      character*(*) string
c
      common /frdraw/mode
      common /frbare/ibare
      common /frwts/iwts(4)
      common /frpens/icolors(3)
c
c     Draw the x- or the y-label.
c
      if (mode.eq.0.and.ibare.ne.1.and.iwts(4).ge.0) then
c
c	  Save current settings.
c
	  call getstatus(idum1,idum2,ipsto,iwsto)
c
          if (iwts(4).gt.0) call weight(iwts(4))
          if (icolors(3).gt.0) call color(icolors(3))
c
	  if (which.eq.1) then
	      call labels(string,' ')
	  else
	      call labels(' ',string)
	  end if
c
c	  Restore settings.
c
	  call color(ipsto)
	  call weight(iwsto)
c
      end if
c
      end
