
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

c     
c     *****************************************************************
c     *								      *
c     *     GET DEVICE:   Get/set graphics device characteristics.    *
c     *								      *
c     *****************************************************************
c     
      subroutine get device
      save
c     
      character*80 device,opt
      common /plot sizes/ xsize,ysize
c     
      common /plot device/ device,aspect,idev
      common /forced/ opt
c     
      common /framesize/ nxpix,nx0,xfac,nypix,ny0,yfac
      common /ncar/ nxpix1,nypix1,nx01,ny01,xfac1,yfac1
      common /plain font/ wid
      common /dev status/ idevon,idevpen,idevwt
      common /dev details/ itek,ivers
c     
      common /hp plot/ ivdef
      character*1 hpinit1(34),hpinit2(34)
c     
      character*1 ctrl(0:31),
     &        null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
      common /ctrlch/ ctrl,
     &        null,ctrla,tab,lf,ff,cr,ctrlx,ctrlz,esc,gs,del
c     
      common /graphproc/ igp
      common /mcpak_colormap/ ncolor,red(0:255),green(0:255),
     &                        blue(0:255)
c     
      common /ps head flag/ iheadflag
      save /ps head flag/
      common /ps keep/ ikeepps,iprintps
      common /ps enforced/ ibounds,ps rmax,ps smax
      common /ps copies/ ncopies
      character*80 psfile,temp
      logical ps_open
c     
      common /sunscreen/ isun
c     
      external mcdxnopen        !$pragma C (mcdxnopen)
      external mcdxinit         !$pragma C (mcdxinit)
      external mcdxquit         !$pragma C (mcdxquit)
c
      common /x input/ interact
      common /xwin init/ no_new_xwin
c     
      data iheadflag/1/no_new_xwin/-1/
c     
      data hpinit1/' ','.','Y',' ','.','P','1',':',
     &        ' ','.','T','1','0','0','0','0',';','1','0','0',
     &        ';',';','1','0','0',';','2','0','0','0',':',
     &        ' ','.','L'/
      data hpinit2/' ','.','@','9','0','0','0',';',':',
     &        'S','P','1',';','I','N',';','R','O','0','0',';',
     &        'I','P',';','C','S','0',';','C','V','1',
     &        ' ','.','Z'/
c     
      do 10001 i=1,34
          if (hpinit1(i).eq.' ')hpinit1(i)=esc
          if (hpinit2(i).eq.' ')hpinit2(i)=esc
10001 continue
c     
      ireaderror=0
      ienv = 0
c     
10101 do 10201 nd=80,1,-1
          if (device(nd:nd).gt.' ')go to 10003
10201 continue
c     
c     No device specified.  See if the environment contains an
c     "MCD_DEVICE" setting.
c     
      if (ienv.eq.0) then
          call mygetenv('MCD_DEVICE', device)
          ienv = 1
c         write(6,'(a,a,a)')'env dev = "',device,'"'
          go to 10101
      end if
      ienv = 2
c     
c     Prompt for a device ID.
c     
1     device=' '
      write(6,'(''Device: ''$)')
      read(5,'(a80)',end=2,err=2)device
      ienv = 0
      do 10002 nd=80,1,-1
          if (device(nd:nd).gt.' ')go to 10003
10002 continue
c     
2     if (ireaderror.gt.3) then
          write(6,*)'Too many errors!'
          stop
      end if
c     
      write(6,3)
3     format(/' Options:'/
     &        ' d: store vectors on disk (PLOT.DAT)'/
     &        ' n: NCAR plotting package (not implemented)'/
     &        ' t01: Tektronix 4010, 10x10'/
     &        ' t02: Tektronix 4010, 10x7.6'/
     &        ' t41: Tektronix 4014, 10x10 - XTERM'/ 
     &        ' t, t42: Tektronix 4014, 10x7.6 - XTERM'/
     &        ' t5:  Tektronix 4105'/ 
     &        ' h: HP plotter 7550a'/
     &        '    ht ==> transparency'/
     &        '    h1 ==> 10x10 display (10x7.7 otherwise)'/
     &        '    hr ==> 90 degree rotation',
     &        ' (--> 10x10 or 10x13.0)'/
     &        ' v0: as t0, but for mac/versaterm'/
     &        ' v, v4: as t4, but for versaterm-pro'/
     &        ' v5: as t5, but for versaterm-pro'/
     &        ' s: SUN window output'/
     &        ' p: PostScript output'/
     &        ' x: X window output (test version)'/
     &        )
c     
c     Devices and associated IDs:
c     --------------------------
c
c		d:		1
c		n:		2
c		t01:		3
c		t02:		4
c		ht:		5
c		hr:		5
c		h1:		6
c		v01:		7
c		v02:		8
c		t41:		9
c		t = t42:	10
c		v41:		11
c		v = v42:	12
c		t5:		13
c		v5:		14
c		s:		15
c		p:		16
c		x:		17
c
      ireaderror=ireaderror+1
      go to 1
c     
10003 if (ienv.eq.1) write(6,'(a)')'Adopted device = "'//
     &        device(1:nd)//'" from the environment'
      nd=nd+1
      device(nd:nd)=' '
c     
c     ONLY CONVERT THE FIRST NONBLANK CHARACTER TO LOWERCASE!!!
c     
4     do 5 i=1,nd
          if (device(i:i).ge.'A'.and.device(i:i).le.'Z')
     &            device(i:i)=char(ichar(device(i:i))+32)
          if (device(i:i).gt.' ') go to 6
5     continue
c     
6     do 10 i=1,nd
          if (device(i:i).gt.' ')go to 15
10    continue
      go to 2
15    do j=i,nd
          device(j-i+1:j-i+1)=device(j:j)
      end do
      nd=nd-i+1
c     
c------------------------------------------------------------------------
c
c     The rest of this routine checks the device and sets the following
c     generally-used variables (as well as any device-specific ones):
c
c     itek	= 1 if we are a Tektronix
c     ivers	= 1 if we are Versaterm
c     wid	= nominal character width
c     aspect	= device aspect ratio
c     ncolor	= number of colors available
c     nxpix	= number of pixels in the x direction
c     nypix	= number of pixels in the y direction
c     nx0	= x offset of the plotting area
c     ny0	= y offset of the plotting area
c     idev	= device ID
c
      itek=0
      ivers=0
      wid=1.
      aspect=1.
      ncolor=2
c     
      if (device(1:1).eq.'d') then
c         
c         Generic data file (not used).
c         ----------------------------
c         
          idev=1
          open(60,file='PLOT.DAT',status='unknown',form='formatted')
          rewind 60
          nx0=0
          nxpix=262143
c         = 64**3 - 1
          ny0=nx0
          nypix=nxpix
c         
      else if (device(1:1).eq.'n') then
c         
c         NCAR graphics (not implemented now).
c         -----------------------------------
c         
          idev=2
c         
c         for centered output (e.g. to lca0):
c         
          nx0=384
          nxpix=32000
c         
c         for output going to lpa0, smaller by 0.877 and
c         shifted left:
c         
c         nxpix=28064
c         nx0=10
c         
          nypix=nxpix
          ny0=nx0
          wid=.86
          nxpix1=877
          nypix1=nxpix1
          nx01=12
          ny01=nx01
c         
      else if (device(1:1).eq.'t') then
c         
c         Tektronix options (slightly buggy?).
c         -----------------------------------
c         
          if (device(2:2).eq.'0') then
c             
c             4010: 1024x781 pixels.
c             
              if (device(3:3).eq.'1') then
c                 
c                 Square output area.
c                 
                  idev=3
                  nxpix=781
              else
c                 
c                 Rectangular output area.
c                 
                  idev=4
                  nxpix=1023
                  aspect=.7625
              end if
c             
          else if (device(2:2).eq.'4'.or.device(2:2).eq.' ') then
c             
c             4010: 4096x3132 pixels.
c             
              if (device(3:3).eq.'1') then
c                 
c                 Square output area.
c                 
                  idev=9
                  nxpix=3132
              else
c                 
c                 Rectangular output area.
c                 
                  idev=10
                  nxpix=4095
                  aspect=.764
              end if
c             
          else if (device(2:2).eq.'5') then
c             
c             Tektronix 4105.
c             --------------
c             
              idev=13
              nxpix=4095
          else
              go to 2
          end if
c         
          itek=1
          nx0=0
          ny0=0
          nypix=nxpix*aspect
c         
      else if (device(1:1).eq.'h') then
c         
c         HP plotter.
c         ----------
c         
          i1=index(device(1:nd),'1')
          if (i1.gt.0)i1=1
          i2=1-i1
          it=index(device(1:nd),'t')
          if (it.gt.0)it=1
          ir=index(device(1:nd),'r')
          if (ir.gt.0)ir=1
          nx0=0
          ny0=0
          idev=5+i1
          if (ir.eq.0) then
              if (i1.eq.1) then
                  nxpix=7840
              else
                  nxpix=10170
                  aspect=.7709
              end if
              nypix=7840
              hpinit2(19)='0'
          else
              if (i1.eq.1) then
                  nypix=7840
              else
                  nypix=10170
                  aspect=1./.7709
              end if
              nxpix=7840
              hpinit2(19)='9'
          end if
          wid=.9
          write(6,'(1x,40a1)')hpinit1
          read(5,*)idummy
          idevpen=1
          hpinit2(12)='1'
          write(6,'(1x,40a1)')hpinit2
          if (it.eq.1)ivdef=1
c         
      else if (device(1:1).eq.'v') then
c         
c         Versaterm-PRO Tektronix emulation (buggy with latest VT release!).
c         -----------------------------------------------------------------
c         
          if (device(2:2).eq.'0') then
c             
c             Versaterm/Mac Tek 4010 emulation apparently allows pixel
c             addresses from 0 to 1023 in both x and y, but maps any
c             y-values above 781 onto 781! The actual output region is
c             a 17.3 cm by 10.1 cm screen.
c             (Similarly for versaterm-pro in 4014 mode with y > 3132.)
c             
c             For versaterm, only bother with 4010 emulation.
c             
              if (device(3:3).eq.'1') then
c                 
c                 Square output area.
c                 
                  idev=7
                  nxpix=593
                  nypix=781
              else
c                 
c                 Rectangular output area.
c                 
                  idev=8
                  nxpix=1023
                  nypix=781
                  aspect=.58
              end if
c             
          else if (device(2:2).eq.'4'.or.device(2:2).eq.' ') then
c             
c             Tektronix 4014 emulation.
c             ------------------------
c             
              if (device(3:3).eq.'1') then
c                 
c                 Square output area.
c                 
                  idev=11
                  nxpix=3132
              else
c                 
c                 Rectangular output area.
c                 
                  idev=12
                  nxpix=4095
                  aspect=.7648
              end if
              nypix=nxpix*aspect
c             
          else if (device(2:2).eq.'5') then
c             
c             Tektronix 4105 emulation.
c             ------------------------
c             
              idev=14
              nxpix=4095
              nypix=4095
          else
              go to 2
          end if
c         
          nx0=0
          ny0=0
          itek=1
          ivers=1
c         
      else if (device(1:1).eq.'s') then
c         
c         SunCore (under Sunview) -- now obsolete and no longer supported.
c         ---------------------------------------------------------------
c         
c         SunCore will not coexist peacefully with X or PostScript!
c         
          if (idev.eq.16) call psquit(2)
          if (idev.eq.17) call mcdxquit
c         
          if (isun.ne.0) then
              write(6,'(''Enter <CR> to delete window'',
     &                '' and reinitialize graphics.'')')
              call plstop
          end if
c         
          if (nd.eq.1) then
              nd=2
              device(2:2)=' '
          end if
c         
c         Most Suncore initialization is done in plinit, which determines
c         the type of frame buffer we have, and the size of the color
c         map, opens a window and initializes the internals of the Core
c         sraphics package.  Note that the interpretation of the input
c         command string is done by plinit, too.
c         
          call plinit('s -b 255 -a 1. -s .5 '//opt
     &            //' '//device(2:nd)//' ',
     &            aspect,icolor,igp,ncolor,ierr)
          if (ierr.ne.0) go to 2
c         
          idev=15
          isun=1
c         
          wid=.85
          idevwt = 1
          idevpen = 1
c         
      else if (device(1:1).eq.'p') then
c         
c         PostScript output.
c         -----------------
c         
c         Only allow one open PostScript file.  Note that opening
c         PostScript will not terminate any X windoes.
c
c         SunCore will not coexist peacefully with X or PostScript!
c         
          if (idev.eq.15) then
              call plstop
              isun = 0
          end if
c
c         If we already have an open PostScript file, append to it
c         (and set graphics parameters from stored values).  If we
c         specify a file name that differs from the previous name,
c         close the current file and open a new one.
c
          if (ps_open()) then
              call ps_filename(device(1:nd)//opt//' ',nd,temp,nt)
              if (temp(1:1).gt.' '
     $                .and.temp(1:nt).ne.psfile(1:npsf)) then
c
c                 A new file name was explicitly specified, and it is
c                 not the same as the current one.  Close out the current
c                 file and get ready to open the new one.
c
                  call psquit(2)
c
              end if
          end if
c
c         Decipher options from command line if necessary
c         -----------------------------------------------
c
          if (.not.ps_open())
     $            call ps_parse(device(1:nd)//opt//' ',nd,
     $                          ikeepps,iprintps,iheadflag,ncopies,
     $                          psfile,npsf,iorient,isparc,psaspect)
c         
c         Offset of origin:
c
c         ***** Beware of magic numbers! *****
c         
          nx0=22
          ny0=nx0
          aspect = psaspect
c         
          if (iorient.eq.1) then
              nxpix = 570
              if (iaspect.gt.0) then
                  nypix = aspect*nxpix
              else
                  nypix = nxpix
                  aspect = 1.
              end if
          else
              nxpix = 750
              if (iaspect.gt.0) then
                  nypix = aspect*nxpix
              else
                  nypix = 510
              end if
              aspect = nypix/float(nxpix)
          end if
c         
          ncolor=256
          psrmax=nx0+nxpix
          pssmax=ny0+nypix
          ibounds=1
c         
c         Attempt to guess a "standard" character width:
c         
          wid=.75
c
c         NOTE change to defaults (overwrite previous PS settings):
c         Note also modification to ps color so color 1 is black in greyscale.
c
          idevpen=1
          idevwt=5
c
c         Note that we do NOT append to an existing file if we open
c         it with psinit!
c
          idev = 16
          if (.not.ps_open())
     $            call psinit(psfile(1:npsf),iorient,isparc)
c         
      else if (device(1:1).eq.'x') then
c         
c         X-windows.
c         ---------
c         
c         SunCore will not coexist peacefully with X or PostScript!
c         
          if (idev.eq.15) then
              call plstop
              isun = 0
          end if
c
          if (mcdxnopen().le.0.or.no_new_xwin.lt.0) then
c
c             Open a new X window.
c         
c             Non-interactive flag for use with scripts.
c         
              interact = 1
              if (index(device(1:nd)//' '//opt,'-i').gt.0
     &                .or. index(device(1:nd)//' '//opt,'-I').gt.0)
     &                interact = 0
c
              xaspect = aspect
              call getaspect(device(1:nd)//' '//opt,iaspect,xaspect)
              call mcdxinit(xaspect,nxpix,nypix,ncolor,ierr)
              if (ierr.gt.1) go to 2
          end if
c
c         (Don't bother with pixels...)
c
          aspect = xaspect
          idev = 17
c         
c         Defaults:
c         
c         background color = black (1)
c         foreground color = white (0)
c         pen width = 0
c
c         NOTE: These may overwrite current X settings with defaults.
c
          call background(0)
          idevpen = 1
          idevwt = 1
          wid=1.
c         
      else
          go to 2
      end if
c     
      end


      subroutine getaspect(string,iaspect,aspect)
      save
      character*(*) string
c     
c     Return the location and value of the aspect ratio in the input
c     string.  Only the first "-a" with a legal number following counts.
c     
      iaspect = index(string,'-a')
      if (iaspect.le.0) iaspect = index(string,'-A')
c     
      if (iaspect.gt.0) then
c         
c         Locate and read the aspect ratio from the input string.
c         
          ibl = 1
          do 10 i=iaspect+2,len(string)
              if (string(i:i).le.' ') then
                  if (ibl.eq.0) then
                      if (index(string(iaspect+2:i),'.').gt.0) then
                          read(string(iaspect+2:i),'(f20.10)',
     &                            iostat=ii)aa
                      else
                          read(string(iaspect+2:i),'(i20)',
     &                            iostat=ii)ia
                          aa = ia
                      end if
                      if (ii.eq.0) then
                          aspect = aa
                          return
                      end if
                  end if
              else
                  ibl = 0
              end if
10        continue
      end if
c     
      end


      subroutine ps_parse(dev,nd,ikeepps,iprintps,iheadflag,ncopies,
     $                    psfile,npsf,iorient,isparc,aspect)
c     
c     Read PostScript parameters from the "device" line.
c     
      character*(*) dev,psfile
c     
c     (1) Print/store:
c     
      ikeepps=0
      iprintps=1
c
      if (index(dev,'-k').gt.0.or.index(dev,'-K').gt.0) ikeepps=1
      if (index(dev,'-n').gt.0.or.index(dev,'-N').gt.0) iprintps=0
c
      if (iprintps.eq.0) ikeepps=1
      if (ikeepps.eq.0) iprintps=1
c     
c     (2) Suppress standard header line:
c     
      if (index(dev,'-h').gt.0.or.index(dev,'-H').gt.0) iheadflag=0
c     
c     (3) Specify # of copies:
c     
      ic=index(dev,'-c')
      if (ic.le.0) ic=index(dev,'-C')
      ncopies=1
      if (ic.gt.0.and.nd.gt.ic+2) then
          inum=0
          do 1600 ii=ic+3,nd
              if (dev(ii:ii).lt.'0'.or.dev(ii:ii).gt.'9') then
                  if (inum.ne.0) then
c                     
c                     NOTE assumption of a trailing blank here!
c                     
                      j=ii-1
                      if (ii.eq.nd) j=ii
                      read(dev(inum:j),*,err=1610,end=1610)
     &                        ncopies
                      go to 1610
                  end if
              else
                  if (inum.eq.0) inum=ii
              end if
1600      continue
      end if
c     
c     (4) PostScript file name:
c     
1610  if=index(dev,'-f')
      if (if.eq.0) if=index(dev,'-F')
      psfile=' '
      npsf=1
      if (if.gt.0) then
          inbl=0
          do 1650 ii=if+3,nd
              if (inbl.eq.0.and.dev(ii:ii).gt.' ') inbl=ii
              if (inbl.gt.0.and.dev(ii:ii).eq.' ') go to 1660
1650      continue
          ii=nd+1
1660      if (inbl.gt.0) then
              npsf=ii-inbl
              psfile(1:npsf)=dev(inbl:ii-1)
              ikeepps=1
              iprintps=0
          end if
      end if
c     
c     (5) Landscape mode?
c     
      iorient = 1
      if (index(dev,'-l').gt.0 .or. index(dev,'-L').gt.0) iorient = 2
c     
c     (6) SPARCprinter point-plotting bug fix:
c     
      isparc = 0
      if (index(dev,'-S').gt.0) isparc = 1
c     
c     (7) Force the output aspect ratio:
c     
      call getaspect(dev,iaspect,aspect)
c     
      end


      subroutine ps_filename(dev,nd,psfile,npsf)
c     
c     Attempt to read a PostScript file name from the "device" line.
c     
      character*(*) dev,psfile
c     
      psfile=' '
      npsf=1
c
      if=index(dev,'-f')
      if (if.eq.0) if=index(dev,'-F')
c
      if (if.gt.0) then
          inbl=0
          do ii=if+3,nd
              if (inbl.eq.0.and.dev(ii:ii).gt.' ') inbl=ii
              if (inbl.gt.0.and.dev(ii:ii).eq.' ') go to 100
          end do
          ii=nd+1
100       if (inbl.gt.0) then
              npsf=ii-inbl
              psfile(1:npsf)=dev(inbl:ii-1)
              ikeepps=1
              iprintps=0
          end if
      end if
c
      end
