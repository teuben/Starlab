c
c     The following is an alphabetized list of all common blocks
c     necessary for the proper functioning of the mcpak graphics
c     library, when graphics contexts are being switched (i.e. these
c     blocks must be saved and restored to preserve continuity).
c
c     NOTE: We need not store device-specific material if only a
c     single instance of a given device is permitted--the info is
c     saved anyway, and cannot be overwritten.
c
      subroutine save_mcpak_context(id)
      save
c
c     Save all relevant commons in a character string.
c
      parameter (NCMAX = 100)
      character*1000 save_string(NCMAX)
c
      common/backgcolor/ibackgcolor
      common/dash/dpatrn(10),dpat,npatrn,ipat,lpen
      common/devdetails/itek,ivers
      common/devinit/initdev
      common/devstatus/idevon,idevpen,idevwt
      common/findex/index
      common/framesize/nxpix,nx0,xfac,nypix,ny0,yfac
      common/frbare/ibare
      common/frbnds/modeb
      common/frconf/scent,rnuml,rnumr,snumt,snumb,dsnums,jrot,stopnum
      common/frdraw/moded
      common/frhts/htl,htn
      common/frint/iframe
      common/frlbx/lmodex
      common/frlby/lmodey
      common/frpens/icolors(3)
      common/frplain/iplain
      common/frrotn/irot
      common/frsetax/kax,lax
      common/frsord/idash
      common/frticks/tiks(3),tikl
      common/frtiklevel/jtiklevel
      common/frwts/iwts(4)
      common/frxnums/xnumbot
      common/frylabpos/slab
      common/lastpoint/rl,sl
      common/lhlimit/rlhe
      common/lowlimit/sbot
      common/mcpak_colormap/ncolor
      common/mlineon/imline
      common/ngonstars/istar
      common/numbron/inumbr
      common/numsymint/rs,ss,dx,dy,sint,cost
      common/penposn/npo
      common/plainfont/wid
      common/plotinvert/inv,ipensto
      common/plotoffset/iin
      common/plotorigin/ro,so
      common/plotsizes/xsize,ysize
      common/sboxdata/iborder,ierase,fraction
      common/scales/xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
      common/strlimits/offxmin,offxmax,offymin,offymax
      common/strposn/iposset,frx,fry
      common/yfollowsx/ixy
c
      save_string(id) = ' '
      write(save_string(id),*)
     $        ibackgcolor,dpatrn,dpat,npatrn,ipat,lpen,
     $        itek,ivers,initdev,idevon,idevpen,idevwt,
     $        index,nxpix,nx0,xfac,nypix,ny0,yfac,ibare,
     $        modeb,scent,rnuml,rnumr,snumt,snumb,dsnums,
     $        jrot,stopnum,moded,htl,htn,iframe,lmodex,
     $        lmodey,icolors,iplain,irot,kax,lax,idash,ixy,
     $        ncolor,
     $        tiks,tikl,jtiklevel,iwts,xnumbot,slab,rl,sl,
     $        rlhe,sbot,imline,istar,inumbr,rs,ss,dx,dy,
     $        sint,cost,npo,wid,inv,ipensto,iin,ro,so,
     $        xsize,ysize,iborder,ierase,fraction,xl,xr,
     $        dinchx,ybot,ytop,dinchy,rlen,slen,offxmin,
     $        offxmax,offymin,offymax,iposset,frx,fry
      return
c
      entry restore_mcpak_context(id)
c
c     Restore a saved graphics context.
c
      read(save_string(id),*)
     $        ibackgcolor,dpatrn,dpat,npatrn,ipat,lpen,
     $        itek,ivers,initdev,idevon,idevpen,idevwt,
     $        index,nxpix,nx0,xfac,nypix,ny0,yfac,ibare,
     $        modeb,scent,rnuml,rnumr,snumt,snumb,dsnums,
     $        jrot,stopnum,moded,htl,htn,iframe,lmodex,
     $        lmodey,icolors,iplain,irot,kax,lax,idash,ixy,
     $        ncolor,
     $        tiks,tikl,jtiklevel,iwts,xnumbot,slab,rl,sl,
     $        rlhe,sbot,imline,istar,inumbr,rs,ss,dx,dy,
     $        sint,cost,npo,wid,inv,ipensto,iin,ro,so,
     $        xsize,ysize,iborder,ierase,fraction,xl,xr,
     $        dinchx,ybot,ytop,dinchy,rlen,slen,offxmin,
     $        offxmax,offymin,offymax,iposset,frx,fry
c
      end
