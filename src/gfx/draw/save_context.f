c
c     The following is an alphabetized list of all common blocks
c     necessary for the proper functioning of mcdraw, when graphics
c     are being switched (i.e. these blocks must be saved and restored
c     to preserve continuity).
c
c     Only non-mcpak arrays are listed here.
c
c     NOTE: We need not store device-specific material if only a
c     single instance of a given device is permitted--the info is
c     saved anyway, and cannot be overwritten.
c
c     NOTE: Most history information will be lost when the graphics
c     context changes.
c
      subroutine save_mcdraw_context(id)
      save
c
c     Save all relevant commons in a character string.
c
      parameter (NCMAX = 100)
      character*1000 save_string(NCMAX)
c
      common/colorlimits/icmin,icmax
      common/compressframes/icompress
      common/dataoffset/delx,dely,delz,facx,facy,facz
      common/drawparams/roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     $                  idevset,jbox,iorig
      common/fileinput/inpmode
      common/frameparams1/xmin,xmax,ymin,ymax,modex,modey
      character*80 xttl,yttl
      common/frameparams2/xttl,yttl
      common/framesize/nxpix,nx0,xfac,nypix,ny0,yfac
      common/initial/roff0,soff0,hn0,hs0,hp0
      common/inputposn/iposn,jposn
      common/localoffset/xoff,yoff
      common/mcdcolor/icolor
      character*80 colormapfile
      common /mcdraw_colormap_file/ colormapfile
      character*1 plot_symbol
      common/mcd_local/ offx,offy,offxsave,offysave,offlabel,
     $                  angle,anglesave,rloc,sloc,
     $                  iweight,iwtsto,jth,jsym,itype,
     $                  ibox,ierbox(0:4),plot_symbol
c
      integer ixttl(80),iyttl(80),icmap(80)
c
c     Note special treatment of character strings xtty and yttl.
c
      nxttl = 0
      nyttl = 0
      ncmap = 0
      do i=80,1,-1
          if (nxttl.eq.0.and.xttl(i:i).gt.' ') nxttl = i
          if (nyttl.eq.0.and.yttl(i:i).gt.' ') nyttl = i
          if (ncmap.eq.0.and.colormapfile(i:i).gt.' ') ncmap = i
      end do
c
      write(save_string(id),*)
     $        offx,offy,offxsave,offysave,offlabel,
     $        angle,anglesave,rloc,sloc,
     $        iweight,iwtsto,jth,jsym,ichar(plot_symbol),itype,
     $        ibox,ierbox,
     $        icmin,icmax,icompress,delx,dely,delz,
     $        facx,facy,facz,roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     $                  idevset,jbox,iorig,
     $        inpmode,xmin,xmax,ymin,ymax,modex,modey,
     $        nxpix,nx0,xfac,nypix,ny0,yfac,roff0,soff0,hn0,hs0,hp0,
     $        iposn,jposn,xoff,yoff,icolor,nxttl,nyttl,ncmap,
     $        (ichar(xttl(i:i)),i=1,nxttl),(ichar(yttl(i:i)),i=1,nyttl),
     $        (ichar(colormapfile(i:i)),i=1,ncmap)
      return
c
      entry restore_mcdraw_context(id)
c
c     Restore a saved graphics context.
c
      read(save_string(id),*)
     $        offx,offy,offxsave,offysave,offlabel,
     $        angle,anglesave,rloc,sloc,
     $        iweight,iwtsto,jth,jsym,jsymbol,itype,
     $        ibox,ierbox,
     $        icmin,icmax,icompress,delx,dely,delz,
     $        facx,facy,facz,roff,soff,aspect1,xlen,ylen,hs,hn,hp,
     $                  idevset,jbox,iorig,
     $        inpmode,xmin,xmax,ymin,ymax,modex,modey,
     $        nxpix,nx0,xfac,nypix,ny0,yfac,roff0,soff0,hn0,hs0,hp0,
     $        iposn,jposn,xoff,yoff,icolor,nxttl,nyttl,ncmap,
     $        (ixttl(i),i=1,nxttl),(iyttl(i),i=1,nyttl),
     $        (icmap(i),i=1,ncmap)
c
c     Reconstruct the strings.
c
      plot_symbol = char(jsymbol)
      do i=1,80
          if (i.le.nxttl) then
              xttl(i:i) = char(ixttl(i))
          else
              xttl(i:i) = ' '
          end if
          if (i.le.nyttl) then
              yttl(i:i) = char(iyttl(i))
          else
              yttl(i:i) = ' '
          end if
          if (i.le.ncmap) then
              colormapfile(i:i) = char(icmap(i))
          else
              colormapfile(i:i) = ' '
          end if
      end do
c
      end
