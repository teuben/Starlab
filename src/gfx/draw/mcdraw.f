
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

        program mcdraw main
c
c       Call the mcdraw routines via a C routine which allocates
c       memory for all internal arrays.
c
        character*200 arg,string
        character c1,c2
        equivalence (c1,arg(1:1)),(c2,arg(2:2))
c
        integer which,s,size
c
c       Defaults:
c
        size = 10000
        string = ' '
        which = 0
c
c       Note use of the generic functions nparams and getparam.
c
        i = 0
        do while (i.lt.nparams())
            i = i + 1
            call getparam(i,arg)
            io = 0
c
            if (c1.eq.'-') then
c
c               Make the second character lowercase.
c
                if (c2.le.'Z') c2 = char(ichar(c2)+32)
c
                if (c2.eq.'n'.or.c2.eq.'s') then
                    i = i + 1
                    call getparam(i,arg)
                    read(arg,'(i10)',iostat=io)s
                    if (io.eq.0.and.s.gt.1.and.s.le.10000000)
     &                      size = s
                else if (c2.eq.'f') then
                    which = 1
                    i = i + 1
                    call getparam(i,string)
                else if (c2.eq.'c'.or.c2.eq.'i'.or.c2.eq.'l') then
                    which = 2
                    i = i + 1
                    call getparam(i,string)
                else if (c2.eq.'h') then
                    call input help
                    go to 99999
                else
                    io = 1
                end if
c
            else
c
                read(arg,'(i10)',iostat=io)s
                if (io.eq.0.and.s.gt.1.and.s.le.10000000)
     &                  size = s
c
            end if
c
            if (io.ne.0) then
                call input help
                go to 99999
            end if
c
        end do
c
        call mcdrawc(size,string,which)
c
99999   end


        subroutine input help
c
c       Print out basic command-line help.
c
        write(6,'(a)')
     &          'Options:   -h(H)	help'
        write(6,'(a)')
     &          '           [-n(s)] #	'//
     &          'specify maximum array length'
        write(6,'(a)')
     &          '           -f(F) file	'//
     &          'take initial input from a file'
        write(6,'(a)')
     &          '           -c(l) line	'//
     &          'specify initial command line'
c       
        end


        block data init mcd data
        save
c
c       Version number:
c
        character*20 version
        common /mcd version/ nversion,version
c
c       Basic mcdraw plotting information:
c       ---------------------------------
c
        common /mcd color/ icolor
        common /initial/ roff0,soff0,hn0,hs0,hp0
        common /draw params/ roff,soff,aspect1,xlen,ylen,
     &                       hs,hn,hp,idevset,jbox,iorig
        common /fr plain/iplain
c
c       Mcdraw frame parameters:
c       -----------------------
c
        common /frame params 1/ xmin,xmax,ymin,ymax,modex,modey
        character*80 xttl,yttl
        common /frame params 2/ xttl,yttl
c
c       History:
c       -------
c
        parameter (NHMAX = 500)
        common /histnums/ lhist(NHMAX),nhist,ishist(NHMAX),
     &                    izhist(NHMAX),ibhist(NHMAX)
c
c       Internal modes:
c       --------------
c
        common /replay/ ireplay
        common /prompt/ iprompt
c
c       File input:
c       ----------
c
        common /file input/ inpmode
c
	common /data offset/ delx,dely,delz,facx,facy,facz
c
c       "Local" data:
c       ------------
c
        character*1 plot_symbol
        common /mcd_local/ offx,offy,offxsave,offysave,offlabel,
     $                     angle,anglesave,rloc,sloc,
     $                     iweight,iwtsto,jth,jsym,itype,
     $                     ibox,ierbox(0:4),plot_symbol
c
        data nversion/3/version/'2.1'/
c
        data icolor/1/
        data roff0/.175/soff0/.125/hn0/.2/hs0/.25/hp0/.1/
        data aspect1/-1./idevset/0/jbox/0/iorig/0/
        data iplain/0/
c       
        data xmin/0./xmax/1./modex/2/xttl/' '/
     &          ymin/0./ymax/1./modey/2/yttl/' '/
c
        data nhist/0/
c
        data ireplay/0/
        data iprompt/1/
        data inpmode/0/
	data delx/0./dely/0./delz/0./facx/1./facy/1./facz/1./
c
        data offx/0./offy/0./offxsave/0./offysave/0./
        data offlabel/0.5/
        data angle/0./anglesave/0./
        data rloc/0./sloc/0./
        data iweight/1/iwtsto/1/
        data jth/0/
        data jsym/0/
        data plot_symbol/' '/
        data itype/0/
        data ibox/0/
        data (ierbox(k),k=0,4)/0,1,1,1,1/
c
        end


        subroutine mcdraw(x,y,z,arr,w,nmax,instring,iwhich)
        save
c       
c       Interactive interface to the MCPAK plotting routines.
c       
        integer nmax,iwhich
        dimension x(1),y(1),z(1),arr(nmax,3),w(1)
        character*(*) instring
c       
c       >>>>>  For now, w really IS of dimension 1!  <<<<<
c              
c       NOTE:  arr and (x, y, z) actually share the same space,
c       so they are effectively equivalenced to one another.
c
c......................................................................
c       
c       General mcpak common information (used "readonly"):
c       --------------------------------------------------
c       
c       Plotting information:
c       --------------------
c       
        character*80 device
        common /plot device/ device,aspect,idev
        common /plot sizes/ xsize,ysize
        common /plot origin/ ro,so
        common /framesize/ nxpix,nx0,xfac,nypix,ny0,yfac
        common /input posn/ iposn,jposn
        common /mcpak_colormap/ ncolors
        common /dev details/ itek,ivers
        common /xwin init/ no_new_xwin
        character*80 colormapfile
        common /mcdraw_colormap_file/ colormapfile
c
c       Eframe options and parameters:
c       -----------------------------
c
        common /scales/ xl,xr,dinchx,ybot,ytop,dinchy,rlen,slen
        common /fr ticks/ tiks(3)
c       
c......................................................................
c
c       Basic mcdraw plotting information:
c       ---------------------------------
c
        common /mcd color/ icolor
        common /initial/ roff0,soff0,hn0,hs0,hp0
        common /draw params/ roff,soff,aspect1,xlen,ylen,
     &                       hs,hn,hp,idevset,jbox,iorig
        common /compress frames/ icompress
        common /fr plain/iplain
c
c       Mcdraw frame parameters:
c       -----------------------
c
        common /frame params 1/ xmin,xmax,ymin,ymax,modex,modey
        character*80 xttl,yttl
        common /frame params 2/ xttl,yttl
c
c       History:
c       -------
c
        parameter (NHMAX = 500)
c       
        character*200 history(NHMAX),lsto
        common /histchars/ history
        common /histnums/ lhist(NHMAX),nhist,ishist(NHMAX),
     &                    izhist(NHMAX),ibhist(NHMAX)
c
c       Internal modes:
c       --------------
c
        common /replay/ ireplay
        common /prompt/ iprompt
c
c       File input:
c       ----------
c
        common /file input/ inpmode
        character*15 modename(0:1)
        common /token delim/ nsep
c
        character*1 comment
        common /file comment/ icomment,comment
c
c       Optional offsets and scales for data (rescale to single precision).
c
	common /data offset/ delx,dely,delz,facx,facy,facz
c
c......................................................................
c       
c       Uninitialized mcdraw common blocks:
c       ----------------------------------
c
c       Help:
c       ----
c
        common /hlevel/ ihelp
        character*2 hindx
        common /hindex/ hindx
c       
        common /resid/ resi
c
c......................................................................
c       
c       Local variables and parameters:
c       ------------------------------
c
        character*1 plot_symbol
        common /mcd_local/ offx,offy,offxsave,offysave,offlabel,
     $                     angle,anglesave,rloc,sloc,
     $                     iweight,iwtsto,jth,jsym,itype,
     $                     ibox,ierbox(0:4),plot_symbol
c
        dimension narr(3)
c
c       Zooming and rescaling:
c       ---------------------
c
        parameter (NZMAX = 50)
        dimension xlzoom(NZMAX),xrzoom(NZMAX),
     &            ybzoom(NZMAX),ytzoom(NZMAX)
c
        parameter (NBMAX = 50)
        dimension xbox(NBMAX),ybox(NBMAX),bscale(NBMAX)
c
        character line*2000,loopstring*2000,
     &            datafile*300,
     &            home*200,stringsto*200,label*200,labelsto*200,
     &            cwd*150,oldcwd*150,tempcwd*150,
     &            input*80,temp*80,keyword*80,tmp*400,
     &            c1*1,c2*1,c3*1,c4*1
c
        character*40 token(20)
        character*120 pwd
c
        dimension xpoly(4),ypoly(4)
        dimension iarg(3),ipat(4)
c       
        integer mysystem
        logical audio,gfxin
c       
        equivalence (n,nx,narr(1)),(ny,narr(2)),(nz,narr(3))
        equivalence (c1,input(1:1)),(c2,input(2:2)),
     $              (c3,input(3:3)),(c4,input(4:4))
c       
        parameter (REDUCE = 0.9)
        parameter (ISEED = 42)
c
        integer num_win, curr_win
        real*8 random
c       
c......................................................................
c
c       Local initialization:
c       --------------------
c
        data cwd(1:1)/'.'/ncwd/1/oldcwd(1:1)/'.'/noldcwd/1/
        data iopen/0/
        data inunit/5/
        data nrange/10000000/
c
        data modename /'fast but stupid','clever but slow'/
c
        data nhfull/0/
        data nstring/1/
        data nzoom/0/
        data nbox/0/
        data nhead0/0/nhead/0/
        data ihmode/0/
        data icsave/1/
        data jreplay/0/
        data npl/1/
        data audio/.true./
        data iherr/0/ihsave/0/
c
        data loop/0/loop1/1/loop2/1/linc/1/
        data nlabel/1/nlabelsto/0/
c
        data rgin/0./sgin/0./
        data narr/3*0/
        data ix/0/iy/0/iz/0/
c
c......................................................................
c       
        hslast = hs
        hnlast = hn
        hplast = hp
c       
        dummy = random(ISEED)
c       
c       do 2 i=1,nmax
c       w(i)=1.
c2      continue
c       
        w(1)=1.
c       *****  For now, w really is of dimension 1.  *****
c       
        nl = -1
c       
c       Initialize some box parameters.
c       
        hn = hn0
        hs = hs0
        hp = hp0
c       
c       Allow the possibility of starting off with file or command-line
c       input instead of using any ~/.mcdrc file.
c       
7       if (iwhich.eq.1) then
            inunit = 7
            open(inunit,file=instring,
     &              status='old',form='formatted',iostat=io)
            if (io.ne.0) then
                write(6,*)'Error reading initial input file.'
                go to 99999
            else
                write(6,*)'Initial input taken from ',instring
            end if
        else if (iwhich.eq.2) then
            line = instring
            nl = len(instring)
        else
c           
c           Is there a standard startup file?
c           
            call mygetenv('HOME',home)
c           
            io1 = 0
            do 5 i=80,1,-1
                if (io1.eq.0.and.inunit.eq.5
     &                  .and.home(i:i).gt.' ') then
                    open(7,file='.mcdrc',status='old',iostat=io)
                    if (io.ne.0) open(7,file=home(1:i)//'/.mcdrc',
     &                      status='old',iostat=io1)
                    if (io1.eq.0) inunit = 7
                end if
5           continue
c           
        end if
c       
c***************************************************************************
c***************************************************************************
c       
9       iscolon=0
c       
c       In operation, iscolon points to the end of the next (sub)command.
c       
c       Note:
c       
c       "go to 10"	==>  get next subcommand
c       "go to 12"	==>  get next command line (use discouraged)
c       "go to 1001"	==>  error handling (such as it is)
c       "go to 99999"	==>  quit
c       
10      if (iscolon.gt.nl) then
c           
            if (iprompt.eq.1) call devoff
c           
c           Get the next command line
c           -------------------------
c           
12          if (loop.ge.0) loop = loop + linc
c
            xloop = float(loop-loop1)*float(loop2-loop)
            if (xloop.lt.0.) then
c
                if (ireplay.gt.0.and.ireplay.le.jreplay) then
c
                    jh = ireplay
                    do while (jh.gt.NHMAX)
                        jh = jh - 1
                    end do
c
                    line = history(jh)
                    nl = lhist(jh)
                    call setisave(ishist(jh))
                    izoom = izhist(jh)
c               
                    write(6,10013)ireplay,line(1:nl)
10013               format('Replay (',i5,'): ',a)
                    ireplay = ireplay + 1
                else
                    ireplay = 0
                    jreplay = 0
                    call getline(inunit,line,nl,iprompt,nhist)
                end if
c
                ihset = 1
c
            else
c
c               We are in a loop, so just use the loop string.
c
                line = loopstring
                nl = nlstring
c
c               Suppress addition to the history list.
c
                ihset = 0
c
            end if
c           
c           Start processing the ENTIRE line
c           --------------------------------
c           
c           Note that historical references are NOT expanded at this stage, 
c           unless it is explicitly requested, or needed for substitution.
c           This is in order to save space, in case of multiple or nested
c           references.
c
c           1a) Expand repeat (!!) references immediately (!! --> !nhist).
c           
            call expandr(line,nl)
c
c           1b) Expand references of the form "!string" (--> !ihist).
c
            call expands(line,nl)
c           
c           2) Strip off trailing blanks and non-significant semicolons.
c           
            call stripbl(line,nl,*10,*12)
c           
c           3) Check for verify (echo-only) mode (trailing \).
c           
            call chktail(line,nl,'\\',noexe,1)
c           
c           4) Check for explicit suppression of recursive expansion
c              (leading ~)
c           
            call chkhead(line,nl,'~',isuppr,1)
c
c           5) Check for forced explicit expansion (leading %).
c           
            call chkhead(line,nl,'%',iexpand,1)
c
c           Note that explicit expansion and substitution are not
c           completely compatible, as the substitution should only be
c           performed after a single level of history expansion.
c           Thus, if recursive expansion is forced, DO IT NOW, and don't
c           worry about its effects on substitutions.
c
            if (iexpand.eq.1) call expandl(line,nl,isuppr,*9)
c               
c           6) Check for substitutions (.....^xxx"yyy" substitutes yyy
c              for xxx).
c
c              Substitutions are performed on a command-by-command basis.
c              Continue making substitutions until none remain.
c
            call applysubs(line,nl,*9)
c           
c           7) Finally, update the history list
c           
            if (ireplay.eq.0.and.ihset.eq.1)
     &              call updhist(line,nl,nzoom,nbox)
c           
            go to 9
c           
        end if
c       
c***************************************************************************
c***************************************************************************
c       
c       Split the line up into individual commands.
c       ------------------------------------------
c       
15      call getnext(line,nl,iscolon,ic0,ic1,input,nin)
c
c       Note that input does NOT have a trailing semicolon.
c       
c       Recursively expand any historical references.
c       --------------------------------------------
c
        do i=1,nin
            if (input(i:i).gt.' ') then
                if (input(i:i).ne.'!') go to 50
                if (i.ge.nin) go to 1001
c               
c               Decode the range reference in the line:
c               
                call rdecode(input(i+1:nin),nhist-1,ihmin,ihmax,*1001)
                go to 35
            end if
        end do
        go to 10
c       
c       Save the rest of the line...
c       
35      nlsto = 0
        if (ic1.lt.nl) then
            nlsto = nl - ic1 + 1
            lsto(1:nlsto) = line(ic1:nl)
        end if
c       
c       (N.B. the leading ";" at position ic1 IS saved here.)
c       
c       Expand the history into the line after ic0.  Note that expandh 
c       does not place a trailing semicolon at the end of the inserted 
c       string, and does not alter nl or ic0.  The process of saving,
c       expanding, and restoring the rest of the line overwrites the
c       original reference.
c
        call expandh(line,nl,ic0,ihmin,ihmax,lextra,iret)
c
        if (iret.ne.0) then
            nhist = nhist - 1
            iscolon = nl + 1
            go to 10
        end if
c       
c       Restore the rest of the line.  (Note: this would be a good
c       place for a consistency check, as the string lengths are not
c       forced to agree...)
c       
        nl = nl - nin + lextra
        if (nlsto.gt.0) then
            line(ic0+1+lextra:nl) = lsto(1:nlsto)
        end if
c       
c       Go back and continue decomposing the line.  This process
c       continues until the next command found contains no historical
c       references.
c       
        iscolon = ic0
        go to 15
c       
c***************************************************************************
c       
c       Strip leading blanks, convert the current command to lowercase,
c       and locate the start of the argument list.
c
50      nino = nin
        call cleanup(input,nin,istart,*10)
c
c       Also shift line left, to reduce the (small) chance of overflow.
c
        if (ic0.gt.0) then
            call shiftstr(line,nl,ic0)
            iscolon = iscolon - ic0
            ic1 = ic1 - ic0
            ic0 = 0
        end if
c       
c***************************************************************************
c***************************************************************************
c       
c       Now decode the command line.
c       ---------------------------
c       
c       Notes:	The entire command is input(1:nin)
c       The argument list is input(istart:nin)
c       The characters c1, c2, and c3 are, respectively,
c       input(1:1), input(2:2), and input(3:3)
c       
c       >>>> (The following long  "if...then...else" construction should be
c       >>>> replaced by a  "case"  statement, if FORTRAN ever gets one...)
c
c       Clean up any pending X graphics commands (force synchronism).
c
        call win_flush
c
        if (noexe.eq.1) then
c           
c           Echo only.
c           
            write(6,*)input(1:nin)
c           
        else if (c1.eq.'=') then
c           
c           =N: Subdivision of the plotting area into quarters.
c           First flush the box stack, if any.
c           
            if (nbox.gt.0) then
                write(temp,'(i4)')nbox
                nt = 4
                if (c2.eq.'0') then
c                   
c                   (Necessary to make sure the screen is cleared.)
c                   
                    temp(5:6) = ' 1'
                    nt = 6
                end if
                call popbox(temp,nt,ierbox,xbox,ybox,bscale,nbox,*1001)
            end if
c           
c           Choose and initialize a new output area.
c           
            call newbox(input,nin,ierbox,*1001)
c           
        else if (c1.eq.'*') then
c           
c           Toggle small/large amount of output.
c           
            iprompt = 1 - iprompt
c           
        else if (c1.eq.'2') then
c           
c           2D: Two-dimensional plot of the specified array.
c               Details in subroutine "ez2dplot."
c           
            call sdecode(input(istart:nin),1,iarg,*1001)
            if (iarg(1).le.0) iarg(1) = 3
            if (narr(iarg(1)).le.0) then
                if (iprompt.eq.1) write(6,*)'No points !'
                go to 1001
            end if
c           
            if (idevset.eq.0.or.xlen.le.0.) call setup(' ','2d')
            if (ibox.eq.0) call box(0.,xlen,0.,ylen)
c           
            call ez2dplot(arr(1,iarg(1)),narr(iarg(1)),
     &              xlen,ylen,iprompt,itype,*1001)
c           
        else if (c1.eq.'3') then
c           
c           3D: Three-dimensional plot of the specified array.
c               Details in subroutine "ez3dplot."
c           
            call sdecode(input(istart:nin),1,iarg,*1001)
            if (iarg(1).le.0) iarg(1) = 3
            if (narr(iarg(1)).le.0) then
                if (iprompt.eq.1) write(6,*)'No points !'
                go to 1001
            end if
c           
            call getorigin(rsave,ssave)
            call setorigin(0.,0.)
            call ez3dplot(arr(1,iarg(1)),narr(iarg(1)),iprompt,*1001)
            call setorigin(rsave,ssave)
c           
        else if (c1.eq.'a') then
c           
            if (c2.eq.' ') then
                go to 1001
c               
            else if (c2.eq.'n') then
c
                if (c3.ne.'-') then
c               
c                   AN: Alter drawing angle (for strings).
c
                    call readrq(input(istart:nin),1,
     &                          aa,dum,dum,dum,*1001)
c
                    anglesave = angle
                    angle = aa
                else
c               
c                   AN-: Restore previous drawing angle (for strings).
c               
                    aa = angle
                    angle = anglesave
                    anglesave = aa
                end if
c               
            else if (c2.eq.'s') then
c               
c               AS: Alter default aspect ratio, relative to the physical
c                   shape of the display.  The value of aspect1 set here
c                   actually goes into effect as soon as the display is
c                   (re)initialized ("de") or resized ("=1," "bo," etc.).
c               
                if (ireplay.eq.0) then
                    call readrq(input(istart:nin),1,
     &                          aspect1,dum,dum,dum,*1001)
                    call newbox('=-0',3,ierbox,*1001)
                end if
c               
            else if (c2.eq.'r') then
c               
c               AR: Draw an arrow, with positions measured in "inches"
c                   or user units, depending on the box status.
c
c               ARH: Specify the size of the arrowhead, in "inches".
c
                if (c3.eq.'h') then
                    call readrq(input(istart:nin),1,
     &                          head_size,dum, dum, dum,*1001)
                    call set_arrow_head(head_size)
                else
                    call readrq(input(istart:nin),4,
     &                          xtail,ytail,xhead,yhead,*1001)
c
                    call arrow(xtail,ytail,xhead,yhead,1-ibox)
                end if
c               
            else if (c2.eq.'u') then
c               
c               AU: Take autocorrelation of specified array (y or z).
c                   Assume that x holds equally-spaced data points.
c               
                call sdecode(input(istart:nin),1,ia,*1001)
c
                if (ia.eq.1) go to 1001
c
                ny0 = narr(ia)
                call cautocorrel(arr(1,ia),narr(ia))
c
c               On return, the specified array has been replaced by its
c               (normalized) autocorrelation function, and reduced in length
c               by a factor of two.  If (as will often be the case), the
c               array is y and x contained time (in equally-spaced steps),
c               simply dividing nx by two will leave an x-array suitable
c               for plotting (i.e. time --> tau).
c
                if (nx.eq.ny0) nx = narr(ia)
                if (iprompt.ne.0) write(6,*)'Array size now = ',narr(ia)
c               
            else
c               
c               A* (etc.): Perform vector arithmetic on two arrays.
c               
                call varith(input,istart,nin,arr,nmax,narr,iprompt,
     &                      *1001)
c               
            end if
c           
        else if (c1.eq.'b') then
c           
            if (c2.eq.'i'.or.c2.eq.'+') then
c               
c               BI/B+: Increase the box size.
c               
                xlen = xlen/REDUCE
                ylen = ylen/REDUCE
c               
            else if ((c2.eq.' '.or.c2.eq.'n'.or.c2.eq.'o'.or.c2.eq.'b')
     &                  .and.(nin.le.3.or.istart.ge.nin)) then
c               
c               B/BB/BN: Draw a box and set up scalings for plotting.
c               
                if (idevset.eq.0) call setup(' ','box')
                call devon
                call clrstr
c               
                call getmod(i1,i2,i3,i4,i5)
                if (c2.eq.'n'.or.c2.eq.'b') then
c                   
c                   Scaling only.
c                   
                    call setmod(1,i2,i3,i4,i5)
                else
                    call setmod(0,i2,i3,i4,i5)
                end if
c               
                call eframe(xmin,xmax,xlen,modex,xttl,
     &                  ymin,ymax,ylen,modey,yttl)
c               
                if (c2.eq.'b') call box(0.,xlen,0.,ylen)
c               
                ibox = 1
                call setmod(i1,i2,i3,i4,i5)
c               
            else if (c2.eq.'c') then
c
c               BC:  Toggle compression of frames in "4-box" format.
c
                icompress = 1 - icompress
c               
            else if ((c2.eq.'o'.or.c2.eq.' ').and.istart.le.nin) then
c               
c               BO:  Offset and rescale the box parameters.
c                    -- Generalization of "=1," etc.
c               
                if (nbox.ge.NBMAX) then
                    if (iprompt.eq.1) write(6,*)' Box stack overflow.'
                    go to 1001
                end if
c               
                if (idevset.eq.0) call setup(' ','bo')
                call offbox(input(istart:nin),nin-istart+1,
     &                      ierbox,xbox,ybox,bscale,nbox,*1001)
c               
            else if (c2.eq.'-') then
c               
c               B-: Pop previous box parameters from the stack.
c               
                call popbox(input(istart:nin),nin-istart+1,
     &                  ierbox,xbox,ybox,bscale,nbox,*1001)
c               
            else if (c2.eq.'s') then
c               
c               BS: Display the box stack.
c               
                if (ireplay.eq.0.and.iprompt.eq.1) then
                    do 80 ib=1,nbox
                        write(6,*)ib,xbox(ib),ybox(ib),bscale(ib)
80                  continue
                end if
c               
            else
                go to 1001
            end if
c           
        else if (c1.eq.'c') then
c           
            if (c2.eq.'o') then
                if (c3.eq.'-') then
c
c                   CO-: Restore the previous plotting color.
c
                    ic = icsave
                else
c               
c                   CO: Set the current plotting color.
c               
                    call readiq(input(istart:nin),1,
     &                          ic,idum,idum,idum,*1001)
                end if
c
                icsave = icolor
                icolor = ic
                call color(icolor)
c               
            else if (c2.eq.'m') then
c               
c               CM: Set or display the color map, if appropriate.
c               
                if (idevset.eq.0)call setup(' ','cm')
                if (idev.lt.15.or.idev.gt.17) go to 1001
c               
                if (ireplay.eq.0) then
c                   
                    idisp = 1
                    ier = 0
                    if (istart.le.nin) then
                        if (idev.eq.15) then
                            call readmap(cwd(1:ncwd)//
     $                                   '/'//input(istart:nin))
                            ncolors = 256
                            idisp = 0
                        else if (idev.eq.17) then
c
c                           Palette names are relative to cwd, unless
c                           name begins with a "/".  In addition, we will
c                           search Starlab if the environment is set.
c
                            if (input(istart:istart).ne.'/') then
c
c                               See if the file exists:
c
                                open(99,file=cwd(1:ncwd)//'/'
     $                                  //input(istart:nin),
     $                                  status='old',iostat=io)
                                if (io.eq.0) then
                                    close(99)
                                    call set_colormap(cwd(1:ncwd)//
     $                                      '/'//input(istart:nin),ier)
                                else
c
c                                   Try Starlab:
c                                   (This mess should be in a separate routine!)
c
                                    temp = ' '
                                    call mygetenv('STARLAB_PATH',temp)
                                    i = len(temp)
                                    do while (i.gt.0)
                                        if (temp(i:i).gt.' ') then
                                            temp(i+1:i+17)
     $                                          = '/src/gfx/palettes'
                                            i = -i - 17
                                        end if
                                        i = i - 1
                                    end do
                                    if (i.ge.0) then
                                        if (iprompt.eq.1) write(6,'(a)')
     $                                      'Color map file not found.'
                                        ier = 1
                                    else
                                        open(99,file=temp(1:-i-1)//'/'
     $                                          //input(istart:nin),
     $                                          status='old',iostat=io)
                                        if (io.eq.0) then
                                            close(99)
                                            call set_colormap(
     $                                          temp(1:-i-1)//'/'//
     $                                          input(istart:nin),ier)
                                        else
                                            if (iprompt.eq.1)
     $                                          write(6,'(a)')
     $                                      'Color map file not found.'
                                            ier = 1
                                        end if
                                    end if
                                end if
                            else
                                call set_colormap(input(istart:nin),ier)
                            end if
                            idisp = 0
                        end if
                        if (ier.eq.0) colormapfile = cwd(1:ncwd)//
     $                                       '/'//input(istart:nin)
                    end if
c
                    if (idisp.eq.1) then
                        call getfill(ifsave)
c
                        dr = .2
                        xpoly(1) = -ro
                        xpoly(2) = xpoly(1) + dr
                        xpoly(3) = xpoly(2)
                        xpoly(4) = xpoly(1)
c
                        ds = ysize/ncolors
                        do 90 j=0,ncolors-1
                            s = j*ds
                            ypoly(1) = s - so
                            ypoly(2) = ypoly(1)
                            ypoly(3) = ypoly(2) + ds
                            ypoly(4) = ypoly(3)
                            call setfill(j)
                            call polyfill(xpoly,ypoly,4)
                            call color(icolor)
                            if (ncolors.lt.64)
     $                              call polydraw(xpoly,ypoly,4)
90                      continue
c
                        if (ifsave.ge.0) then
                            call setfill(ifsave)
                        else
                            call unsetfill
                        end if
                    end if
c
                end if
c               
            else if (c2.eq.'b') then
c               
c               CB: Set the background color.
c               
                if (idev.ne.15.and.idev.ne.17) go to 1001
c
                call readiq(input(istart:nin),1,
     &                      jcolor,idum,idum,idum,*1001)
                call background(jcolor)
c               
            else if (c2.eq.'d') then
c               
c               CD: Change working directory.
c
                if (c3.eq.'-') then
                    tempcwd = cwd
                    ntempcwd = ncwd
                    cwd = oldcwd
                    ncwd = noldcwd
                    oldcwd = tempcwd
                    noldcwd = ntempcwd
                else
                    if (istart.gt.nin) then
                        oldcwd = cwd
                        noldcwd = ncwd
                        cwd = '.'
                        ncwd = 1
                    else
                        if (input(istart:istart).eq.'/') then
                            oldcwd = cwd
                            noldcwd = ncwd
                            cwd(1:nin-istart+1) = input(istart:nin)
                            ncwd = nin-istart+1
                        else if (input(istart:istart).eq.'~') then
c
c                           Must explicitly expand this reference.
c
                            call mygetenv('HOME',tempcwd)
                            do i=len(cwd),1,-1
                                if (tempcwd(i:i).gt.' ') go to 100
                            end do
                            if (iprompt.eq.1)
     $                              write(6,*)'Can''t expand ~'
                            go to 1001
c
100                         oldcwd = cwd
                            noldcwd = ncwd
                            ncwd = i
                            cwd(1:ncwd) = tempcwd(1:ncwd)
                            if (istart.lt.nin) then
                                cwd(ncwd+1:ncwd+nin-istart+1)
     $                                  = input(istart+1:nin)
                                ncwd = ncwd + nin-istart
                            end if
                        else
                            oldcwd = cwd
                            noldcwd = ncwd
                            cwd(ncwd+1:ncwd+nin-istart+2)
     &                              = '/'//input(istart:nin)
                            ncwd = ncwd+nin-istart+2
                        end if
                    end if
                end if
c               
            else if (c2.eq.'u') then
c               
c               CU: Replace an array by its cumulative sum.
c               
                call sdecode(input(istart:nin),1,iarg,*1001)
                if (iarg(1).le.0) go to 1001
                if (narr(iarg(1)).le.0) go to 1001
c               
                do i=2,narr(iarg(1))
                    arr(i,iarg(1)) = arr(i-1,iarg(1)) + arr(i,iarg(1))
                end do
c               
            else if (c2.eq.'x'.or.c2.eq.'y'.or.c2.eq.'z') then
c               
c               CX/Y/Z: Read a standard image into the specified array.
c               
                call decode(c2,ia,*1001)
c               
                if (input(istart:istart).eq.'/') then
c                   
c                   Absolute path name.
c                   
                    datafile = input(istart:nin)
                    ndata = nin-istart+1
                else
c                   
c                   Relative path name.
c                   
                    datafile = cwd(1:ncwd)//'/'//input(istart:nin)
                    ndata = nin-istart+1+ncwd+1
                end if
c               
                call ctox(datafile(1:ndata),arr(1,ia),narr(ia),nmax)
c               
                if (iprompt.eq.1) write(6,*)narr(ia),' points read'
c               
            else if (c2.eq.' ') then
c               
c               C: Read x and y columns from the input file.
c               
                if (iopen.eq.0) then
                    call devoff
                    write(6,*)'c: No file open'
                    go to 10
                end if
c
                call readiq(input(istart:nin),2,
     &                            ix,iy,idum,idum,*1001)
                if (ix.le.0.or.iy.le.0.or.ix.eq.iy) go to 1001
c               
                call readcols(10,nhead,nrange,ix,iy,x,y,z,nx,nmax,
     &                  c1,1,iprompt,*1001)
                ny = nx
c               
            else
                go to 1001
            end if
c           
        else if (c1.eq.'d') then
c           
            if (c2.eq.' '.or.c2.eq.'r'.or.c2.eq.'a') then
c               
c               D: Draw a line to the specified point.
c               
                if (nin.lt.istart) then
                    xloc = rgin
                    yloc = sgin
                else
                    call readrq(input(istart:nin),2,
     &                          xloc,yloc,dum,dum,*1001)
c
                    if (c2.ne.'a') then
c                       
c                       User coordinates.
c                       
                        if (ibox.ne.0)
     &                          call fr inches(xloc,yloc,xloc,yloc)
                    else
c                       
c                       DA: Absolute coordinates.
c                       
                        xloc = xloc - ro
                        yloc = yloc - so
                    end if
                end if
c
                call segment(rloc,sloc,xloc,yloc)
                rloc = xloc
                sloc = yloc
c               
            else if (c2.eq.'d') then
c               
c               DD: Draw a dashed line to the specified point.
c               DDA: Draw a dashed line to the specified absolute point.
c               
                call dplot(rloc,sloc,3)
                call readrq(input(istart:nin),2,
     &                      xloc,yloc,dum,dum,*1001)
c
                if (c3.ne.'a') then
c                       
c                   DD: User coordinates.
c
                    if (ibox.ne.0)
     &                      call fr inches(xloc,yloc,xloc,yloc)
                else
c                       
c                   DDA: Absolute coordinates.
c                       
                    xloc = xloc - ro
                    yloc = yloc - so
                end if
c
                call dplot(xloc,yloc,2)
                rloc = xloc
                sloc = yloc
c               
            else if (c2.eq.'e') then
c                   
c               DE: Specify (new) device and initialize.
c                   
                if (ireplay.eq.0) then
                    if (nin.lt.istart) then
                        nin=istart
                        input(istart:nin)=' '
                    end if
                    call setup(input(istart:nin),' ')
                    call nobounds
                end if
c               
            else if (c2.eq.'i') then
c
c               DI(FF): Differentiate second array with respect to the first,
c                       placing the result in the third (default: second).
c
                call sdecode(input(istart:nin),3,iarg,*1001)
c
                if (iarg(1).eq.0) then
c
c                   Default is y(x) --> y.
c
                    iarg(1) = 1
                    iarg(2) = 2
                else if (iarg(2).le.0) then
                    if (iarg(1).eq.1) then
                        iarg(2) = 2
                    else
                        iarg(2) = 1
                    end if
                end if
c
                if (iarg(3).le.0) iarg(3) = iarg(2)
c
                if (narr(iarg(1)).le.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No points !'
                    end if
                    go to 1001
                end if
c
c               Note that array iarg(1) is the independent variable,
c               whose length determines the number of elements used.
c
                call differentiate(arr(1,iarg(1)),arr(1,iarg(2)),
     &                             arr(1,iarg(3)),narr(iarg(1)),
     &                             iprompt,*1001)
                narr(iarg(3)) = narr(iarg(1))
c
            else if(c2.eq.'x'.or.c2.eq.'y'.or.c2.eq.'z') then
c               
c               DX(Y,Z): Take 10**(x, y, or z).
c               
                call decode(c2,ia,*1001)
                do 150 i=1,narr(ia)
                    arr(i,ia)=10.**arr(i,ia)
150             continue
c               
            end if
c           
        else if (c1.eq.'e') then
c           
            if (c2.eq.'a'.or.(c2.eq.' '.and.jbox.eq.0)) then
c               
c               EA: Erase entire screen.
c               
159             if (idevset.eq.0)call setup(' ','erase')
                call devon
                call clear
                ibox = 0
c
c               Explicit "go to" here because of the slightly illegal
c               jump into this loop from below.
c
                go to 10
c               
            else if (c2.eq.' ') then
c               
c               E: Erase current plotting area.
c               
160             if (idevset.eq.0.or.ibox.eq.0)call setup(' ','erase')
                call devon
                if (jbox.eq.0) then
                    call clear
                else
                    call erase(-roff,.5*xsize-roff,
     &                         -soff,.5*ysize-soff)
                end if
c
c               Explicit "go to" here because of the slightly illegal
c               jump into this loop from below.
c
                go to 10
c               
            else if (c2.eq.'c') then
c               
c               EC: Echo input.
c               
                if (ireplay.eq.0) write(6,*)input(istart:nin)
c               
            else if (c2.eq.'b') then
c               
c               EB: Erase interior of box.
c               
                if (ibox.eq.0)go to 1001
                call devon
                call erase(tiks(3),xlen-tiks(3),tiks(3),ylen-tiks(3))
c               
            else if (c2.eq.'l') then
c               
c               EL: Erase label area.
c               
                call erase(0.,xlen,ylen+.01,ylen+(offlabel+2.)*hs)
c               
            else if (c2.eq.'r') then
                if (c3.eq.'r'.or.c3.eq.'p') then
                    if (n.le.0) go to 1001
c                   
c                   ERR: Error bars (from z):
c                   ERP: Error bars (from z), limited by current point size:
c
                    if (c3.eq.'p'.and.jth.ne.0) then
                        hh = hp
                    else
                        hh = 0.
                    end if
                    call errors(input,istart,nin,x,y,z,n,hh,*1001)
c                   
                else if (c3.eq.'c') then
                    if (n.le.0) go to 1001
c                   
c                   ERC: Specify error-bar cap size.
c                   
                    call readrq(input(istart:nin),1,
     &                          capsize,dum,dum,dum,*1001)
c
                    call setercap(capsize)
c                   
                else
c                   
c                   ER: Erase all or part of the display.
c                   
                    do 165 i=nin,istart,-1
                        if (input(i:i).gt.' ') go to 166
165                 continue
c                   
                    if (jbox.eq.0) then
                        go to 159
                    else
                        go to 160
                    end if
c                   
166                 call readrq(input(istart:nin),4,
     &                          x1,x2,y1,y2,*1001)
c
                    if (idevset.eq.0) call setup(' ','erase')
                    if (ibox.eq.1) call fr inches(x1,y1,x1,y1)
                    if (ibox.eq.1) call fr inches(x2,y2,x2,y2)
                    call devon
                    call erase(x1,x2,y1,y2)
                end if
c               
            else
c               
c               EX(Y,Z): Exponentiation of x, y, or z.
c               
                call decode(c2,ia,*1001)
                do 170 i=1,narr(ia)
                    arr(i,ia)=exp(arr(i,ia))
170             continue
c               
            end if
c           
        else if (c1.eq.'f') then
c           
            if (c2.eq.' '.or.c2.eq.'i') then
c               
c               F: Open an input file (as unit 10).
c               
                if (iopen.eq.1) close(10)
c               
                inb=0
                do 171 i=nin,istart,-1
                    if (input(i:i).eq.char(92)) then
c
c                       (This is a backquote...).
c
                        if (inb.eq.0) then
                            nhead0 = 0
                        else
                            call readiq(input(i+1:nin),1,
     &                                  nhead0,idum,idum,idum,*1001)
                        end if
                        nhead = nhead0
                        nin = i - 1
                        go to 172
                    end if
                    if (input(i:i).gt.' ') inb=inb+1
171             continue
172             if (input(istart:istart).eq.'/') then
c                   
c                   Absolute path name.
c                   
                    datafile = input(istart:nin)
                    ndata = nin-istart+1
                else
c                   
c                   Relative path name.
c                   
                    datafile = cwd(1:ncwd)//'/'//input(istart:nin)
                    ndata = nin-istart+1+ncwd+1
                end if
c               
                open(10,file=datafile(1:ndata),status='old',
     &                  form='formatted',err=1001)
                iopen=1
c               
            else if (c2.eq.'c') then
c               
c               FC: Close the currently open data file.
c               
                if (iopen.eq.1) then
                    close(10)
                    iopen = 0
                end if
c               
            else if (c2.eq.'o') then
c               
c               FO: Toggle use of fancy font.
c               
                iplain = 1 - iplain
c               
            else if (c2.eq.'1'.or.c2.eq.'2') then
c               
c               F1(2): Perform a linear least-squares fit to y(x).
c               
                if (n.le.1) go to 1001
                if (c3.eq.'p'.and.ibox.eq.0) go to 1001
c               
                if (c2.eq.'1') then
                    call fitxy(input,istart,nin,x,y,z,w,n,iprompt,ier)
                else
                    call fitxyz(input,istart,nin,x,y,z,w,n,iprompt,ier)
                end if
                if (ier.ne.0) go to 1001
c               
            else if (c2.eq.'p'.or.c2.eq.'t') then
c               
c               FP(T): Perform a general (trig) polynomial fit to y(x)
c               (equal weighting).
c               
                if (ibox.eq.0.or.n.le.1) go to 1001
c               
                call genfit(input,istart,nin,x,y,z,n,1,iprompt,*1001)
c               
            else if (c2.eq.'z') then
c               
c               FZ: Perform a general polynomial fit to y(x) (z-weighting).
c               
                if (ibox.eq.0.or.n.le.1) go to 1001
c               
                call genfit(input,istart,nin,x,y,z,n,n,iprompt,*1001)
c               
            else
                go to 1001
            end if
c           
        else if (c1.eq.'g') then
c           
            if (c2.eq.'s') then
c               
c               GS: Read a character string from stdin or the save stack.
c               
                call getstring(input,istart,nin,stringsto)
c               
c               Get the length of the string.
c               
                do 175 i=200,1,-1
                    if (stringsto(i:i).gt.' ') go to 176
175             continue
                i = 1
c               
176             nstring = i
c               
            else
c               
c               G: Get graphic input.
c               
                if (idevset.eq.0) call setup(' ','gin')
c               
                if (c2.eq.'c'.and.(idev.lt.15.or.idev.gt.17))
     &                  go to 1001
c               
                if (itek.eq.1.and.ivers.eq.0.and.iprompt.eq.1)
     &                  write(6,*)'Click mouse and hit <CR> for ',
     &                            'graphic input in Tek mode.'
c               
                call getgfx(rgin,sgin)
                call devoff
c               
                if (c2.ne.'c') then
c                   
c                   Location:
c                   
                    if (idev.lt.15.or.idev.gt.17) then
                        if (ibox.eq.0) then
                            write(6,179)iposn,jposn,rgin,sgin
179                         format(' Pixels ',2i5,
     &				    ',  coordinates ',2f7.3)
                        else
                            call fr users(rgin,sgin,xx,yy)
                            write(6,10179)iposn,jposn,rgin,sgin,xx,yy
10179                       format(' Pixels ',2i5,',  coords ',
     &                              2f7.3,',  values ',1p2e13.5)
                        end if
                    else
                        if (ibox.eq.0) then
                            write(6,20179)ro+rgin,so+sgin,rgin,sgin
20179                       format(' Position ',2f7.3,
     &                              ',  coordinates ',2f7.3)
                        else
                            call fr users(rgin,sgin,xx,yy)
                            write(6,21179)ro+rgin,so+sgin,rgin,sgin,
     &                              xx,yy
21179                       format(' Position',2f7.3,', coords',
     &                              2f7.3,', values',1p2e13.5)
                        end if
                    end if
c                   
                else
c                   
c                   GC: Color input from mouse:
c                   
                    icolor=ncolors*(sgin+so)/ysize
                    write(6,'('' color = '',i3)')icolor
                    call color(icolor)
c                   
                end if
c               
c               GM(D): Move/draw options:
c               
                if (c2.eq.'d') call segment(rloc,sloc,rgin,sgin)
c               
                if (c2.eq.'d'.or.c2.eq.'m') then
                    rloc=rgin
                    sloc=sgin
                end if
c               
            end if
c           
        else if (c1.eq.'h'.or.c1.eq.'?') then
c           
            if (c2.eq.' '.or.c2.eq.'e') then
c               
                if (ireplay.eq.0) then
c                   
c                   H: Print out some helpful information.
c                   
                    call devoff
                    ihelp = 1
                    if (c1.eq.'h') ihelp = 2
                    if (istart.gt.nin) then
                        hindx = '  '
                    else
                        if (istart.lt.nin) then
                            hindx = input(istart:istart+1)
                        else
                            hindx(1:1) = input(istart:istart)
                            hindx(2:2) = ' '
                        end if
                    end if
                    call help(' ')
                end if
c               
            else if (c1.eq.'?') then
c               
                go to 1001
c               
            else if (c2.eq.'h') then
c               
c               HH: Full help information!
c               
                call fullhelp
c               
            else if (c2.eq.'k') then
c               
c               HK: Keyword  help.
c               
                if (ireplay.eq.0) then
c                   
                    call devoff
                    ihelp = 2
c
                    if (istart.gt.nin) then
                        hindx = '  '
                        keyword = ' '
                    else
                        keyword = input(istart:nin)
                    end if
c
                    call help(keyword)
c
                end if
c
            else if (c2.eq.'+') then
c               
c               H+: Audio (HAL) help on.
c               
                audio = .true.
c               
            else if (c2.eq.'-') then
c               
c               H-: Audio (HAL) help off.
c               
                audio = .false.
c               
            else if (c2.eq.'?') then
c               
c               H?: Height information.
c               
                call devoff
                write(6,*)'hp = ',hp
                write(6,*)'hn = ',hn
                write(6,*)'hs = ',hs
c               
            else if (c2.eq.'g') then
c               
                if (c3.eq.'e') then
c                   
c                   HGE: Toggle/set histogram error-bar mode.
c                   
                    if (istart.gt.nin) then
                        iherr = 1 - iherr
                    else
                        call readiq(input(istart:nin),1,
     &                              iherr,idum,idum,idum,*1001)
                        if (iherr.ne.0) iherr = 1
                    end if
c                   
                else if (c3.eq.'m') then
c                   
c                   HGM: Toggle/set histogram display mode.
c                   
                    if (istart.gt.nin) then
                        ihmode = 1 - ihmode
                    else
                        call readiq(input(istart:nin),1,
     &                              ihmode,idum,idum,idum,*1001)
                        if (ihmode.ne.0) ihmode = 1
                    end if
c                   
                else if (c3.eq.'s') then
c                   
c                   HGS: Toggle/set histogram-save mode.
c                   
                    if (istart.gt.nin) then
                        if (ihsave.eq.0) then
                            ihsave = 1
                        else
                            ihsave = 0
                        end if
                    else
                        call readiq(input(istart:nin),1,
     &                              ihsave,idum,idum,idum,*1001)
                        if (ihsave.ne.0.and.ihsave.ne.2) ihsave = 1
                    end if
c                   
                else
c                   
c                   HG: Draw a histogram.
c                   
                    call histogram(input,istart,nin,x,y,z,nx,ny,nz,
     &                      iherr,ihmode,ihsave,iprompt,*1001)
                    ibox = 1
                end if
c               
            else if (c2.eq.'i') then
c               
                if (ireplay.eq.0) then
c                   
c                   HI: List history.
c                   
                    if (nin.lt.istart) then
                        i1=1
                        i2=nhist
                    else
                        call rdecode(input(istart:nin),nhist-1,
     &                          i1,i2,*1001)
                        i1=max(1,min(nhist,i1))
                        i2=max(i1,min(nhist,i2))
                    end if
                    call devoff
                    do 210 ih=i1,i2
                        jh=ih
209                     if (jh.gt.NHMAX) then
                            jh=jh-1
                            go to 209
                        end if
                        write(lsto(1:6),'(i4,'': '')')ih
                        write(6,*)lsto(1:6),history(jh)(1:lhist(jh))
210                 continue
                end if
c               
            else if (c2.eq.'l') then
c               
c               HL: Specify horizontal y-axis label.
c               
                call getmod(i1,i2,i3,i4,i5)
                call setmod(i1,i2,i3,i4,0)
c               
            else if (c2.eq.'n') then
c               
c               HN: Set number height.
c               
                call modifyh(hnlast,hn,c3,input(istart:nin),*1001)
                call sethts(hs,hn)
c               
            else if (c2.eq.'p') then
c               
c               HP: Set point size.
c               
                call modifyh(hplast,hp,c3,input(istart:nin),*1001)
c               
            else if (c2.eq.'s') then
c               
c               HS: Set character height.
c               
                call modifyh(hslast,hs,c3,input(istart:nin),*1001)
                call sethts(hs,hn)
c               
            else
                go to 1001
            end if
c           
        else if (c1.eq.'i') then
c
            if (c2.eq.'n'.and.c3.eq.'t') then
c
c               INT: Integrate second array with respect to the first.
c
                call sdecode(input(istart:nin),3,iarg,*1001)
c
                if (iarg(1).eq.0) then
c
c                   Default is y(x) --> y.
c
                    iarg(1) = 1
                    iarg(2) = 2
                else if (iarg(2).le.0) then
                    if (iarg(1).eq.1) then
                        iarg(2) = 2
                    else
                        iarg(2) = 1
                    end if
                end if
c
                if (iarg(3).le.0) iarg(3) = iarg(2)
c
                if (narr(iarg(1)).le.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No points !'
                    end if
                    go to 1001
                end if
c
                call integrate(arr(1,iarg(1)),arr(1,iarg(2)),
     &                         arr(1,iarg(3)),narr(iarg(1)),
     &                         iprompt,*1001)
c
                write(6,*)' Integral = ',arr(narr(iarg(1)),iarg(3))
c
            else if (c2.eq.'c') then
c
c               IC: Allow/disallow comments in input line (smart input only).
c
                icomment = 1 - icomment
c
                if (icomment.ne.0) then
c
                    if (istart.gt.nin) then
c                       
c                       Default is '%'
c                       
                        comment = '%'
                    else
c
c                       Use first nonblank character as comment indicator.
c
                        comment = ' '
                        do i = istart,nin
                            if (comment.eq.' '.and.input(i:i).gt.' ')
     &                              comment = input(i:i)
                        end do
                        if (comment .eq. ' ') comment = '%'
                    end if
c
                    if (iprompt.ne.0)
     &                   write(6,*)'Comment character is "',comment,'"'
c
                end if
c
            else if (c2.eq.'d') then
c
c               ID: Wait for graphic input before continuing.
c                   (Only wait for devices with a visible display).
c
                call devoff
                if (gfxin()) then
                    write(6,'(a)')'Click right-hand mouse button in'//
     &                            ' graphics window to continue'
                    call graphin(r,s)
                else if (idev.gt.2.and.idev.ne.5.and.idev.ne.6
     &                      .and.idev.ne.16) then
                    write(6,'(a)')'Enter <CR> to continue'
                    read(5,'(a)')c1
                end if
c
            else if (c2.eq.'m') then
c
c               IM: Toggle input mode to accept non-numeric columns.
c
                inpmode = 1 - inpmode
                if (iprompt.ne.0) write(6,*)'Input mode = ',inpmode,
     &                                      ' ('//modename(inpmode)//')'
c
            else if (c2.eq.'n'.and.ireplay.eq.0) then
c               
c               IN: Take input from file.
c               
                if (istart.le.nin) then
                    if (inunit.ne.5) close(inunit)
                    inunit = 7
                    if (input(istart:istart).eq.'/') then
c                       
c                       Absolute path name.
c                       
                        datafile = input(istart:nin)
                        ndata = nin-istart+1
                    else
c                       
c                       Relative path name.
c                       
                        datafile = cwd(1:ncwd)//'/'//
     &                             input(istart:nin)
                        ndata = nin-istart+1+ncwd+1
                    end if
c                   
                    open(inunit,file=datafile(1:ndata),
     &                      status='old',form='formatted',err=1001)
                    if (iprompt.eq.1)
     &                      write(6,*)' Reading input from ',
     &                      datafile(1:ndata)
                else if (inunit.ne.5) then
                    close(inunit)
                    inunit = 5
                end if
c               
            else if (c2.eq.'w') then
c
c               IW: Iconify an X window.
c
                if (num_win(17).le.0) then
                    if (iprompt.ne.0) write(6,'(a)')'No X windows!'
                else
c
c                   Iconify an X window.
c
                    if (istart.le.nin) then
                        call readiq(input(istart:nin),1,
     &                              iwin,idum,idum,idum,*1001)
c
                        call iconify_win(iwin)
                    else
                        call iconify_all
                    end if
                end if
c
            else
c               
                if (ireplay.eq.0) then
                    call devoff
                    if (istart.gt.nin) then
c                       
c                       I: Print information on x, y and/or z.
c                       
                        call iminmax(x,narr(1),ix1,ix2)
                        call iminmax(y,narr(2),iy1,iy2)
                        call iminmax(z,narr(3),iz1,iz2)
                        call devoff
                        write(6,196)'x',nx,x(max(1,ix1)),ix1,
     &                          x(max(1,ix2)),ix2,
     &                          xmin,xmax
                        write(6,196)'y',ny,y(max(1,iy1)),iy1,
     &                          y(max(1,iy2)),iy2,
     &                          ymin,ymax
                        write(6,196)'z',nz,z(max(1,iz1)),iz1,
     &                          z(max(1,iz2)),iz2
196                     format(1x,a1,' (',i6,'):  min, max = ',1p,
     &                          2(e12.4,' (',i6,') '):/
     &                          14x,'plot:',6x,e12.4,10x,e12.4)
                    else
c                       
c                       Information on specified array.
c                       
                        call decode(input(istart:istart),ia,*1001)
                        call iminmax(arr(1,ia),narr(ia),i1,i2)
                        write(6,196)char(ia+119),narr(ia),
     &                          arr(max(1,i1),ia),i1,
     &                          arr(max(1,i2),ia),i2
                    end if
                end if
c               
            end if
c           
        else if (c1.eq.'j') then
c           
c           J: Specify point type/connection.
c           
            call readiq(input(istart:nin),1,
     &                  jth,idum,idum,idum,*1001)
c           
        else if (c1.eq.'k') then
c           
c           K: Delete an X-window.
c
            if (idev.ne.17) then
                if (iprompt.ne.0) write(6,'(a)')'Not using X!'
            else
                call readiq(input(istart:nin),1,
     &                      iwin,idum,idum,idum,*1001)
                call kill_win(iwin)
                if (num_win(idev).le.0) idevset = 0
            end if
c           
        else if (c1.eq.'l') then
c           
c           L: Get/set limits.
c           
            if (c2.eq.' '.or.c2.eq.'i') then
c               
                if (istart.gt.nin) then
                    if (n.le.0) go to 1001
                    call minmax(x,n,xmin,xmax)
                    call minmax(y,n,ymin,ymax)
                else
                    call gettokens(input(istart:nin),token,nt)
                    if (n.le.0.) then
                        xlo = 0.
                        xhi = 0.
                        ylo = 0.
                        yhi = 0.
                    else
                        call minmax(x,n,xlo,xhi)
                        call minmax(y,n,ylo,yhi)
                    end if
c
                    if (nt.eq.1.and.token(1)(1:1).eq.'+') then
c
c                       Same limits on both axes, +/- L.
c
                        call readrtoken(token(1),xmax,xhi)
                        xmin = -xmax
                        ymin = xmin
                        ymax = xmax
                    else
                        if (nt.gt.0) call readrtoken(token(1),xmin,xlo)
                        if (nt.gt.1) call readrtoken(token(2),xmax,xhi)
                        if (nt.gt.2) call readrtoken(token(3),ymin,ylo)
                        if (nt.gt.3) call readrtoken(token(4),ymax,yhi)
                    end if
                end if
                nzoom = 0
c               
            else if (c2.eq.'a') then
c               
c               LA: Add a label.
c               
                if (istart.gt.nin.and.nstring.le.0) go to 1001
c
                call strpos(.5,0.)
                if (idevset.eq.0) call setup(' ','la')
                call devon
                if (istart.le.nin) then
c
c                   Use the specified label.
c
                    nlabel = nin-istart+1
                    label(1:nlabel) = input(istart:nin)
c
                else if (nlabelsto.gt.0) then
c
c                   Use the last label.
c
                    nlabel = nlabelsto
                    label(1:nlabel) = labelsto(1:nlabelsto)
c
                else
c
c                   Just use the contents of the current string.
c
                    nlabel = nstring
                    label(1:nlabel) = stringsto(1:nstring)
c
                end if
c
c               Save the label, for future use.
c
                nlabelsto = nlabel
                labelsto(1:nlabelsto) = label(1:nlabel)
c
c               Check for font modifications:
c
                ntmp = nlabel
                tmp(1:ntmp)=label(1:nlabel)
                if (iplain.eq.1) call plain_sim(tmp, ntmp)
c
                call simbol(.5*xlen,ylen+offlabel*hs,hs,
     &                      tmp(1:ntmp)//'%%',0.,999)
                call clrstr
c               
            else if (c2.eq.'o') then
c
                if (c3.eq.'o') then
c
c                   LOOP: Loop over the rest of the line (after the next ";").
c
                    loop1 = 0
                    loop2 = 0
                    linc = 1
c
                    if (istart.le.nin) then
                        call readiq(input(istart:nin),3,
     &                              loop1,loop2,linc,idum,*1001)
c
                        if (loop1.gt.0.and.loop2.le.0) then
c
c                           Special syntax:
c
                            loop2 = loop1
                            loop1 = 1
                        end if
c
                        if (loop1.le.0) loop1 = 1
                        if (loop2.le.0) loop2 = 1000000
                    else
                        loop1 = 1
                        loop2 = 1000000
                    end if
c
                    if (loop1.gt.0) then
c
                        nlstring = nl - iscolon
                        loopstring(1:nlstring) = line(iscolon+1:nl)
                        iscolon = nl + 1
c
                        if (loop1.le.loop2) then
                            linc = abs(linc)
                        else
                            linc =  -abs(linc)
                        end if
c
                        loop = loop1 - linc
c
                    end if
c
                else
c               
c                   LO: Label vertical offset.
c                   
                    call readrq(input(istart:nin),1,
     &                          offlabel,dum,dum,dum,*1001)
c
                end if
c
            else if (c2.eq.'=') then
c           
c               L=: Subdivision of the plotting area into quarters,
c                   according to the loop counter.
c
c               First flush the box stack, if any.
c           
                if (nbox.gt.0) then
                    write(temp,'(i4)')nbox
                    nt = 4
                    if (c2.eq.'0') then
c                   
c                       (Necessary to make sure the screen is cleared.)
c                   
                        temp(5:6) = ' 1'
                        nt = 6
                    end if
                    call popbox(temp,nt,ierbox,xbox,ybox,bscale,
     &                          nbox,*1001)
                end if
c
                ib = (loop - loop1) / linc + 1
                do while (ib.gt.4)
                    ib = ib - 4
                end do
c           
c               Choose and initialize a new output area.
c
                write(temp(1:2),'(''='',i1)')ib
                call newbox(temp,2,ierbox,*1001)
c           
            else if (c2.eq.'c') then
c               
c               LC: Set color from loop counter.
c               
                ic = loop
                if (ncolors*(ic/ncolors).eq.ic) ic = ic + 1
c
                icsave = icolor
                icolor = ic
                call color(icolor)
c               
            else if (c2.eq.'n'.and.c3.ne.'x'.and.c3.ne.'y'
     &                        .and.c3.ne.'z') then
c
c               LN: Set number of n-gon sides from loop counter.
c
                jsym = loop
                plot_symbol = ' '
                do while (jsym.gt.20)
                    jsym = jsym - 20
                end do
c
            else if (c2.eq.'p') then
c               
c               LP: Print loop counter.
c               
                write(6,*)'Loop = ',loop
c               
            else if (c2.eq.'w') then
c               
c               LW: Set weight from loop counter.
c               
                iwt = 10*(loop-1)
                do while (iwt.gt.50)
                    iwt = iwt - 50
                end do
c               
                if (iwt.le.0) iwt = 1
c               
                iwtsto = iweight
                iweight = iwt
                call weight(iweight)
c               
            else if (c2.eq.'y') then
c               
c               LY: Read the y array from a column specified by
c                   the loop counter.  Exit the loop on error.
c               
                if (iopen.eq.0) then
                    if (loop.gt.0) loop = 0
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No file open'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c
                iy = loop
c
                call readcols(10,nhead,nrange,iy,iy,x,y,z,ny,nmax,
     &                        c2,0,iprompt,*1001)
c               
            else if (c2.eq.'g'.or.c2.eq.'n') then
c
c               LG(N): Logarithms.
c               
                call decode(c3,ia,*1001)
                nerr=0
                do 200 i=1,narr(ia)
                    s=arr(i,ia)
                    if (s.le.0.) then
                        nerr=nerr+1
                        arr(i,ia)= -50.
                    else
                        if (c2.eq.'n') then
                            arr(i,ia)=log(arr(i,ia))
                        else
                            arr(i,ia)=log10(arr(i,ia))
                        end if
                    end if
200             continue
                if (nerr.gt.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,'(i5,'' error(s)'')')nerr
                    end if
                end if
c
            else
                go to 1001
            end if
c           
        else if (c1.eq.'m') then
c           
            if (c2.eq.'o') then
c               
c               MO: Specify plot modes.
c               
                call readiq(input(istart:nin),2,
     &                      mx,my,idum,idum,*1001)
                modex = mx
                modey = my
                if (my.eq.0) modey = modex
c
            else if (c2.eq.'v') then
c               
c               MV: Move an array.
c               
                ibl=1
                narg=0
                do 250 i=istart,nin
                    if (input(i:i).eq.' '.or.input(i:i).eq.',') then
                        ibl=1
                    else
                        if (ibl.eq.1) then
                            narg=narg+1
                            call decode(input(i:i),iarg(narg),*1001)
                            if (narg.eq.2) then
                                if (narr(iarg(1)).le.0) go to 1001
                                do 240 ii=1,narr(iarg(1))
                                    arr(ii,iarg(2))=arr(ii,iarg(1))
240                             continue
                                narr(iarg(2)) = narr(iarg(1))
                                go to 251
                            end if
                            ibl=0
                        end if
                    end if
250             continue
251             continue
c               
            else if (c2.eq.' '.or.c2.eq.'a') then
c               
c               M: Move to the specified point.
c               
                if (nin.lt.istart) then
                    rloc=rgin
                    sloc=sgin
                else
                    call readrq(input(istart:nin),2,
     &                          rloc,sloc,dum,dum,*1001)
c
                    if (c2.ne.'a') then
                        if (ibox.ne.0)
     &                          call fr inches(rloc,sloc,rloc,sloc)
                    else
                        rloc = rloc - ro
                        sloc = sloc - so
                    end if
                end if
            end if
c           
        else if (c1.eq.'n') then
c           
            if (c2.eq.'c') then
c
c               NC: Toggle use of comma as delimiter.
c                   Default value of nsep is 3, make it 2
c                   to allow commas in get_token.
c
                nsep = 5 - nsep
c
                if (iprompt.ne.0) then
                    if (nsep.eq.2) then
                        write(6,*)'Commas ignored in input data'
                    else
                        write(6,*)'Commas accepted as delimiters'
                    end if
                end if
            else if (c2.eq.'p') then
c               
c               NP: Plot a polygon at the specified point.  Size, number
c                   of sides, and fill color are the same as for mline.
c               
                if (nin.lt.istart) then
                    rloc=rgin
                    sloc=sgin
                else
                    call readrq(input(istart:nin),2,
     &                          rloc,sloc,dum,dum,*1001)
                    if (ibox.ne.0) call fr inches(rloc,sloc,rloc,sloc)
                end if
c               
                call ngon(rloc,sloc,.5*hp,jsym,0.)
c               
            else if (c2.eq.'x'.or.c2.eq.'y'.or.c2.eq.'z') then
c               
c               NX(Y,Z): Explicitly set narr.
c               NXY: set nx and ny to the same value
c               NXYZ: set nx, ny, and nz to the same value
c
                call readiq(input(istart:nin),1,
     &                      newn,idum,idum,idum,*1001)
c
                if (c2.eq.'x'.and.c3.eq.'y') then
                    nn = 2
                    if (c4.eq.'z') nn = 3
                    do ia=1,nn
                        narr(ia) = newn
                    end do
                else
                    call decode(c2,ia,*1001)
                    narr(ia) = newn
                end if
c               
            else
c               
c               N: Set polygon style.
c
                if (istart.le.nin) then
c
c                   n s ==> plot symbols at each point.
c                   n # ==> set jsym (polygons, spokes, etc.)
c
                    if ((input(istart:istart).lt.'0'
     $                      .or.input(istart:istart).gt.'9')
     $                    .and.input(istart:istart).ne.'-') then
                        plot_symbol = input(istart:istart)
                    else
                        plot_symbol = ' '
                        call readiq(input(istart:nin),1,
     &                              jsym,idum,idum,idum,*1001)
                        if (jsym.ge.0) then
                            call unsetstars
                        else
                            call setstars
                        end if
                    end if
                end if
            end if
c           
        else if (c1.eq.'o') then
c
            if (c2.ne.'-') then
c
                if (c2.eq.'f') then
c
c                   OF: set data delsets.
c
                    delx = 0.
                    dely = 0.
                    delz = 0.
                    call readrq(input(istart:nin),3,
     &                          delx,dely,delz,dum,*1001)
c
                else
c           
c                   O: Set string offset.
c           
                    call readrq(input(istart:nin),2,
     &                          fx,fy,dum,dum,*1001)
c
                    if (ierr.ne.0) then
                        call clrstr
                    else
                        offxsave = offx
                        offysave = offy
                        offx = fx
                        offy = fy
                        call strpos(offx,offy)
                    end if
                end if
c
            else
                dum = offxsave
                offx = offxsave
                offxsave = dum
                dum = offysave
                offy = offysave
                offysave = dum
                call strpos(offx,offy)
            end if
c           
        else if (c1.eq.'p') then
c           
            if (c2.eq.'a') then
c               
c               PA: New page.
c               
                call devon
                call newpage
c               
            else if (c2.eq.'f') then
c               
c               PF: Set polygon fill index.
c               
                call readiq(input(istart:nin),1,
     &                      i1,idum,idum,idum,*1001)
c
                if (i1.ge.0) then
                    call setfill(i1)
                else
                    call unsetfill
                end if
c
            else if (c2.eq.'o') then
c
c               PO: Plot points only (same as "j -1", but more obvious)
c
                jth = -1
c
            else if (c2.eq.'r') then
c               
c               PR(C): Print current Postscript and (optionally) close
c		       PS output.
c
                if (idev.eq.16) then
                    if (c3.eq.'c') then
                        call psquit(2)
c
c                       In this case, we have to reinitialize the graphics
c                       device -- either X or from scratch.
c
                        if (num_win(17).gt.0) then
c
c                           Restore the last-used X-window.
c
                            if (iprompt.gt.0) write(6,'(a)')
     $                              'Reverting to X-window output.'
c
                            idev = 17
                            iwin = curr_win()
                            no_new_xwin = iwin
                            call setup('x','psquit')
                            no_new_xwin = -1
                            call set_win(iwin, ierr)
                        else
                            idevset = 0
                            call setup(' ', 'psquit')
                        end if
                    else
                        call psquit(1)
                    end if
                end if
c
            else if (c2.eq.'s') then
c
                if (c3.eq.'c') then
c
c                   PSC: Toggle use of color in PostScript:
c
                    call set ps color
c
                    if (iprompt.eq.1) then
                        call get ps color(ipsc)
                        if (ipsc.eq.0) then
                            write(6,*)'Using ',ncolors,
     &                                '-greyscale PostScript'
                        else
                            write(6,*)'Using ',ncolors,
     &                                '-color PostScript'
                        end if
                    end if
c
                else if (c3.eq.'q') then

c                   Very like 'prc', but don't attempt to repoen
c                   a graphics device of none is left.

                    if (idev.eq.16) then
                        call psquit(2)
c
c                       In this case, we have to reinitialize the graphics
c                       device -- either X or from scratch.
c
                        if (num_win(17).gt.0) then
c
c                           Restore the last-used X-window.
c
                            if (iprompt.gt.0) write(6,'(a)')
     $                              'Reverting to X-window output.'
c
                            idev = 17
                            iwin = curr_win()
                            no_new_xwin = iwin
                            call setup('x','psquit')
                            no_new_xwin = -1
                            call set_win(iwin, ierr)
                        else
                            idevset = 0
                        end if
                    end if

                else if (c3.eq.'r') then

c                   PSR: Force rewriting of long Postscript files
c                        to place BoundingBox at start.

                    call set ps write

                else
c
c                   PS: Replace y(x) by its power spectrum, x by frequency.
c                       Assume that the data are evenly spaced.
c
                    call powerspectrum(x,y,nx,iprompt)
                    ny = nx
c
                end if
c
            else if (c2.eq.'w'.and.c3.eq.'d') then
c               
c               PWD: UNIX "pwd" lookalike.
c               
                call mygetenv('PWD',pwd)
                do i=len(pwd),1,-1
                    if (pwd(i:i).gt.' ') go to 300
                end do
                i = 1
                pwd(1:1) = '?'
300             write(6,*)cwd(1:ncwd),' (. = ',pwd(1:i),')'
c               
            else if (c2.eq.'z') then
c               
c               PZ: Plot y(x), clipped according to z.
c               
                if (n.le.0) then
                    if (iprompt.eq.1) write(6,*)'No points !'
                    go to 1001
                end if
c               
                if (ibox.le.0) then
                    if (iprompt.eq.1) write(6,*)'No box !'
                    go to 1001
                end if
c               
                call xyzplot(input,istart,nin,
     &                       x,y,z,n,jth,jsym,plot_symbol,hp,*1001)
c               
            else
c               
c               P: Plot y(x).
c               
                if (n.le.0) then
                    if (iprompt.eq.1) write(6,*)'No points !'
                    go to 1001
                end if
c               
                if (ibox.le.0) then
                    if (iprompt.eq.1) write(6,*)'No box !'
                    go to 1001
                end if
c               
                if (c2.eq.'c') then
                    call xyplotc(input,istart,nin,x,y,z,n,
     &                      itype,jth,jsym,plot_symbol,hp,*1001)
                else
                    call xyplot(input,istart,nin,x,y,z,n,
     &                      itype,jth,jsym,plot_symbol,hp,*1001)
                end if
c               
            end if
c           
        else if (c1.eq.'q') then
c
            if (c2.eq.'u') then
c
c               QU(IET): suppress most output (like '*' toggle).
c
                iprompt = 0
c
            else
c           
c               Q: Exit graphics/mcdraw.
c
                call devoff
                if (c2.ne.'1') then
                    call mcquit
                else
                    call mcquit1
                end if
c
                if (c2.ne.'c') go to 99999
                idevset = 0
c
            end if
c           
        else if (c1.eq.'r') then
c           
            if (c2.eq.'a') then
c               
c               RA: Put a ramp from 1 to n in the specified array.
c               
                call sdecode(input(istart:nin),2,iarg,*1001)
                if (iarg(1).le.0) iarg(1) = 1
c               
                if (narr(iarg(1)).le.0.and.iarg(2).le.0) then
                    if (iprompt.eq.1) then
                        write(6,*)'Array length unknown.'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c               
c               Specifying a second argument sets the value of narr.
c               
                if (iarg(2).gt.0) narr(iarg(1)) = iarg(2)
c               
                do 340 i=1,narr(iarg(1))
                    arr(i,iarg(1)) = i
340             continue
c               
            else if (c2.eq.'n') then
c               
c               RN: Put random numbers in the specified array.
c               
                call sdecode(input(istart:nin),2,iarg,*1001)
                if (iarg(1).le.0) iarg(1) = 1
c               
                if (narr(iarg(1)).le.0.and.iarg(2).le.0) then
                    if (iprompt.eq.1) then
                        write(6,*)'Array length unknown.'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c               
c               Specifying a second argument sets the value of n.
c               
                if (iarg(2).gt.0) narr(iarg(1)) = iarg(2)
c               
                do 341 i=1,narr(iarg(1))
                    arr(i,iarg(1)) = random(0)
341             continue
c               
            else if (c2.eq.'e') then
c               
                if (c3.eq.'b') then
c                   
c                   REB: Rebin the specified array, via a C routine
c                   to allocate space.
c                   
                    call crebin(input(istart:nin),arr,nmax,narr,
     &                      istat,iprompt)
                    if (istat.ne.0) go to 1001
c                   
                else if (c3.eq.'s') then
c                   
c                   RES: Reset most plotting parameters to their initial states.
c                   Note that the device and the plot offsets are untouched.
c                   
                    xmin = 0.
                    xmax = 1.
                    modex = 2
                    xttl = ' '
                    ymin = 0.
                    ymax = 1.
                    modey = 2
                    yttl = ' '
c                   
                    nzoom = 0
c                   
                    close(10)
                    iopen = 0
                    nx = 0
                    ny = 0
                    nz = 0
                    nhead0 = 0
                    nhead = 0
c                   
                    itype = 0
                    ibox = 0
                    jth = 0
                    jsym = 0
                    plot_symbol = ' '
                    npl = 1
c                   
                    rgin = 0.
                    sgin = 0.
                    nstring = 0
                    icolor = 1
                    call color(icolor)
c                   
                    call newbox('=-0',3,ierbox,*1001)
                    hp = .1
                    hn = .2
                    hs = .25
                    call sethts(hs,hn)
                    iplain = 0
                    call setpln
                    call nobounds
c                   
                    offlabel = .5
                    rloc = 0.
                    sloc = 0.
c                   
                    cwd = '.'
                    ncwd = 1
c                   
                else if (c3.eq.'p') then
c                   
c                   REP: Replay commands, with modifications and omissions.
c                   
                    if (ireplay.eq.0) then
                        if (istart.le.nin) then
                            call rdecode(input(istart:nin),nhist,
     &                              i1,i2,*1001)
                            i1=max(1,min(nhist,i1))
                            i2=max(i1,min(nhist,i2))
                        else
                            i1=1
                            i2=nhist
                        end if
                        ireplay = i1
                        jreplay = i2
                        iscolon = nl + 1
                    end if
                end if
c               
            else if (c2.eq.'p') then
c               
c               RP: Print register.
c               
                write(6,*)'Register = ',register
c               
            else if (c2.eq.'r') then
c               
c               RR: Specify read ranges (i1 to i2, after header).
c               
                call gettokens(input(istart:nin),token,nt)
                i1 = 1
                if (nt.gt.0) call readitoken(token(1),i1,1)
                i2 = 10000000
                if (nt.gt.1) call readitoken(token(2),i2,10000000)
                nhead = nhead0 + i1 - 1
                nrange = i2 - i1 + 1
c               
            else if (c2.eq.'s') then
c               
c               RS: Set register.
c               
                call readrq(input(istart:nin),1,
     &                      register,dum,dum,dum,*1001)
c               
            else if (c2.eq.'x'.or.c2.eq.'y'.or.c2.eq.'z') then
c               
c               RX,Y,Z: Set register from x,y,z.
c               
                call readiq(input(istart:nin),1,
     &                      ireg,idum,idum,idum,*1001)
c
                iarr = ichar(c2) - ichar('x') + 1
                register = arr(max(1,min(n,ireg)),iarr)
c               
            end if
c           
        else if (c1.eq.'s') then
c           
            if (c2.eq.' '.or.(c2.eq.'t'.and.c3.eq.'r')) then
c               
c               S: Plot string at current location.
c               
                if (istart.le.nin) then
                    nstring=nin-istart+1
                    stringsto(1:nstring)=input(istart:nin)
c
c                   Check for font modifications:
c
                    ntmp = nstring
                    tmp(1:ntmp)=stringsto(1:nstring)
                    if (iplain.eq.1) call plain_sim(tmp, ntmp)
c
                    call boxsim(rloc,sloc,hs,
     &                      tmp(1:ntmp)//'%%',angle,999)
c
                else if (nstring.gt.0) then
c
c                   Check for font modifications:
c
                    ntmp = nstring
                    tmp(1:ntmp)=stringsto(1:nstring)
                    if (iplain.eq.1) call plain_sim(tmp, ntmp)
c
                    call boxsim(rloc,sloc,hs,
     &                      tmp(1:ntmp)//'%%',angle,999)
                end if
c               
            else if (c2.eq.'a') then
c               
c               SA: Save x, y and z in specified file.
c               
                if (max(nx,ny,nz).le.0) go to 1001
c               
                open(1,file=input(istart:nin),status='new',
     &                  form='formatted',iostat=io)
                if (io.ne.0) then
                    write(6,*)'Error opening new ',input(istart:nin)
                    go to 1001
                end if
                do i=1,max(nx,ny,nz)
                    if (i.gt.nx) x(i) = x(nx)
                    if (i.gt.ny) y(i) = y(ny)
                    if (i.gt.nz) z(i) = z(nz)
                    write(1,*)i,x(i),y(i),z(i)
                end do
                close(1)
c               
            else if (c2.eq.'b') then
c               
c               SB:  Specify string-box parameters.
c               SB?: Print string-box parameters.
c
                if (c3.eq.'?') then
                    if (iprompt.ne.0) then
                        call sboxset(ib,ie,fr)
                        write(6,'(a,i1,a,i1,a,f7.3)')
     &                          ' box = ',ib,'  erase = ',ie,
     &                          '  fraction = ',fr
                    end if
                else
                    call readrq(input(istart:nin),3,
     &                          xib,xie,fr,dum,*1001)
                    ib = nint(xib)
                    if (ib.ne.0) ib = 1
                    ie = nint(xie)
                    if (ie.ne.0) ie = 1
                    if (fr.lt.0.) fr = 0.
                    call sboxset(ib,ie,fr)
                end if
c               
            else if (c2.eq.'c') then
c
c               SC: set data scaling.
c
                facx = 1.
                facy = 1.
                facz = 1.
                call readrq(input(istart:nin),3,
     &                      facx,facy,facz,dum,*1001)
c
            else if (c2.eq.'l') then
c               
c               SL: Sleep for a specified number of seconds
c                   (useful for pausing in loops).
c               
                call readiq(input(istart:nin),1,
     &                      isl,idum,idum,idum,*1001)
                call uwait(1000000*isl)
c               
            else if (c2.eq.'m') then
                if (c3.eq.'o') then
c                   
c                   SMO: Smooth y --> z.
c                   
                    if (n.le.0.or.istart.gt.nin) go to 1001
c                   
                    call readrq(input(istart:nin),3,
     &                          dtsmooth,xioption,xiwindow,dum,*1001)
c
                    ioption = nint(xioption)
                    iwindow = nint(xiwindow)
                    call smooth(x,y,z,n,dtsmooth,ioption,iwindow)
                    nz = n
                else
c                   
c                   SM: Make box smaller.
c                   
                    red = REDUCE
                    call readrq(input(istart:nin),1,
     &                          red,dum,dum,dum,*1001)
                    xlen = red*xlen
                    ylen = red*ylen
                end if
c               
            else if (c2.eq.'o') then
c               
                if (c3.eq.' '.or.c3.eq.'-') then
c                   
c                   SO: Sort an array.
c                   
                    call sdecode(input(istart:nin),1,iarg,*1001)
                    if (iarg(1).le.0) iarg(1) = 1
                    if (narr(iarg(1)).le.0) go to 1001
c                   
                    if (c3.eq.'1')
     &                      call negate(narr(iarg(1)),arr(1,iarg(1)))
c                   
                    call sort(narr(iarg(1)),arr(1,iarg(1)))
c                   
                    if (c3.eq.'1')
     &                      call negate(narr(iarg(1)),arr(1,iarg(1)))
c                   
                else if (c3.eq.'2') then
c                   
c                   SO2: Sort an array, carrying another along with it.
c                   
                    call sdecode(input(istart:nin),2,iarg,*1001)
                    if (iarg(1).le.0) iarg(1) = 1
                    if (narr(iarg(1)).le.0) go to 1001
                    if (iarg(2).le.0) iarg(2) = iarg(1)
c                   
                    call sort2(narr(iarg(1)),arr(1,iarg(1)),
     &                      arr(1,iarg(2)))
c                   
                else
                    go to 1001
                end if
c               
            else if (c2.eq.'s') then
c               
c               SS: Plot an X/SUN/PS string at the current point.
c               
                if (istart.le.nin) then
                    call symbl(rloc,sloc,hs,
     &                      input(istart:nin)//'%%',angle,999)
                    nstring=nin-istart+1
                    stringsto(1:nstring)=input(istart:nin)
                else if (nstring.gt.0) then
                    call symbl(rloc,sloc,hs,
     &                      stringsto(1:nstring)//'%%',angle,999)
                end if
c               
            else if (c2.eq.'t') then
c               
c               STAT: Dump status of mcdraw.
c
                if (iprompt.eq.1) then
c
                    call devoff
c
c                   Directory, etc:
c
                    write(6,'(/a)')' working directory = "'//
     &                            cwd(1:ncwd)//'"'
                    call mygetenv('PWD',pwd)
                    do i=len(pwd),1,-1
                        if (pwd(i:i).gt.' ') go to 400
                    end do
                    i = 1
                    pwd(1:1) = '?'
400                 write(6,*)' (. = ',pwd(1:i),')'


                    write(6,'(a)')' data file name    = "'//
     &                            datafile(1:ndata)//'"'
c
c                   Data input.
c
                    write(6,'(a,i1,a)')' data input mode   = ',
     &                      inpmode,' ('//modename(inpmode)//')'

                    if (inpmode.eq.1) then
                        if (nsep.eq.2) then
                            write(6,'(a)')
     &                              ' commas ignored in input data'
                        else
                            write(6,'(a)')' commas accepted as'//
     &                                    ' delimiters in input data'
                        end if
                        if (icomment.ne.0)
     &                          write(6,'(a)')
     &                          ' input comment character is "'//
     &                          comment//'"'
                    end if
c
c                   Graphics device:
c
                    nd = len(device)
                    do while (nd.gt.0.and.device(nd:nd).le.' ')
                        nd = nd - 1
                    end do
                    if (nd.le.0) nd = 1

                    write(6,'(/a)')' graphics device = "'//
     &                             device(1:nd)//'"'
c
                    write(6,'(a,f4.1,a,f4.1,a,i4,a,i4)')
     &                      ' xsize = ',xsize,'  ysize = ',ysize,
     &                      '  nxpix = ',nxpix,'  nypix = ',nypix
                    write(6,'(a,f5.2,a,f5.2)')
     &                      ' device aspect ratio = ',aspect,
     &                      '  mcdraw aspect ratio = ',aspect1
c
                    write(6,'(/a,i4,a,i4)')' colors = ',ncolors,
     &                                     '  current color = ',icolor
                    if (idev.eq.15.or.idev.eq.17) then
                        ncmap = 0
                        do i=len(colormapfile),1,-1
                            if (colormapfile(i:i).gt.' '
     $                              .and.ncmap.eq.0) then
                                ncmap = i
                                write(6,'(a)')' color map file = '
     $                                  //colormapfile(1:i)
                            end if
                        end do
                    end if
c
                    nwin = 0
                    do id = 1,20
                        nwin = nwin + num_win(id)
                    end do
                    write(6,'(a,i5,a,i5,a)')
     $                      ' current window = ',curr_win(),
     $                      '  (total number = ',nwin,')'
                    write(6,'(a,l1)')' audio = ',audio
c
c                   Plotting:
c
                    write(6,'(/a,f6.2,a,f6.2,a)')
     &                      ' box origin = (',ro,', ',so,')'
                    write(6,'(a,f5.2,a,f5.2,a,i2,a,i2)')
     &                      ' xlen = ',xlen,'  ylen = ',ylen,
     &                      '  modex = ',modex,'  modey = ',modey
                    write(6,'(3(a,f7.3),a,i1)')
     &                      ' hs = ',hs,'  hn = ',hn,'  hp = ',hp,
     &                      '  plain = ',iplain

                    write(6,'(3(a,i5),a)')' weight = ',iweight,
     &                      '  jth = ',jth,
     $                      '  ngons: ',jsym,'   '//plot_symbol
                    write(6,'(a,4i4)')' pattern: ',ipat
c
                    write(6,*)
                    write(6,'(3(a,i2)/)')' ix = ',ix,'  iy = ',iy,
     &                                   '  iz = ',iz
                    call iminmax(x,narr(1),ix1,ix2)
                    call iminmax(y,narr(2),iy1,iy2)
                    call iminmax(z,narr(3),iz1,iz2)
                    write(6,196)'x',nx,x(max(1,ix1)),ix1,
     &                           x(max(1,ix2)),ix2,
     &                           xmin,xmax
                    write(6,196)'y',ny,y(max(1,iy1)),iy1,
     &                           y(max(1,iy2)),iy2,
     &                           ymin,ymax
                    write(6,196)'z',nz,z(max(1,iz1)),iz1,
     &                           z(max(1,iz2)),iz2
c
c                   Miscellaneous:
c
                    write(6,'(/a,5i3)')' ierbox = ',ierbox
                    if (nbox.gt.0) write(6,'(a,3(a,f5.2))')' box:  ',
     &                      'xbox = ',xbox(nbox),
     &                      '  ybox = ',ybox(nbox),
     &                      '  scale = ',bscale(nbox)
c
                    write(6,'(a,f5.2)')' label = "'//label(1:nlabel)//
     &                                 '"  offset = ',offlabel
                    write(6,'(a)')' string = "'//
     &                            stringsto(1:nstring)//'"'
                    write(6,'(a,f7.2,a,2f5.2)')' angle = ',angle,
     &                                       '  offsets: ',offx,offy
                    write(6,'(a,1p,e14.6)')' register = ',register
c
                    write(6,*)
c
                end if
c
            else if (c2.eq.'w') then
c               
c               SW: Swap two arrays.
c               
                if (n.le.0)go to 1001
                call sdecode(input(istart:nin),2,iarg,*1001)
                if (iarg(1).le.0) iarg(1) = 1
                if (iarg(2).le.0) iarg(2) = 2
c               
                call swap(arr(1,iarg(1)),arr(1,iarg(2)),n)
c               
            else if (c2.eq.'y') then
c               
c               SY: Execute a UNIX command in the current ("pwd") directory..
c               
                if (ireplay.eq.0.and.istart.le.nin) then
                    iret = mysystem('cd '//cwd(1:ncwd)//'; '//
     &                      input(istart:nin)//'; exit')
                    if (iprompt.eq.1) write(6,*)'Return status = ',iret
                end if
c               
            else
c               
c               S+(-,*,/): Perform scalar arithmetic on an array.
c               
                call sarith(input,istart,nin,arr,nmax,narr,register,
     &                  iprompt,*1001)
c               
            end if
c           
        else if (c1.eq.'t') then
c
            if (c2.eq.'x'.or.c2.eq.'y'.or.c2.eq.'z') then
c               
c               TX,Y,Z: Type first few elements of x,y,z.
c
                if (istart.le.nin) then
                    call readiq(input(istart:nin),1,
     &                          ntype,idum,idum,idum,*1001)
                else
                    ntype = 10
                end if
c
                iarr = ichar(c2) - ichar('x') + 1
                call devoff
                write(6,'(1p,6e12.4)')(arr(i,iarr),i=1,ntype)
c
            else if (c2.eq.'i') then
c
c               TI: Specify tick level on frame.
c
                call readiq(input(istart:nin),1,
     &                      itik,idum,idum,idum,*1001)
                if (itik.ne.2.and.itik.ne.3) itik = 1
                call settiks(itik)
c
            else if (c2.eq.' ') then
c               
c               T: Specify line type.
c           
                if (istart.gt.nin) then
                    itype = 0
                    ipat(1) = 0
                    ipat(2) = 0
                    ipat(3) = 0
                    ipat(4) = 0
                else
                    call readiq(input(istart:nin),4,
     &                          i1,i2,i3,i4,*1001)
                    itype = 1
                    call setpat(i1,i2,i3,i4)
                    ipat(1) = i1
                    ipat(2) = i2
                    ipat(3) = i3
                    ipat(4) = i4
                end if
            end if
c           
        else if (c1.eq.'u') then
c           
c           U: Unplot the previous plot (if possible).
c           
            if (ibox.eq.0) go to 1001
c           
            call invert
            call unplot(x,y,z,n,itype,jth,jsym,plot_symbol,hp)
            call invert
c           
        else if (c1.eq.'v') then
c           
            if (c2.eq.'l') then
c               
c               VL: Vertical y-axis label.
c               
                call getmod(i1,i2,i3,i4,i5)
                call setmod(i1,i2,i3,i4,1)
c               
            else if (c2.eq.' ') then
c               
c               V: Invert subsequent plot colors.
c               
                if (idevset.eq.0)call setup(' ','v')
                call devon
                call invert
c               
            else
                go to 1001
            end if
c           
        else if (c1.eq.'w') then
c
            if (c2.eq.'-'.or.c2.eq.' ') then
c
                if (c2.eq.'-') then
c
c                   W-: restore line weight.
c
                    iwt = iwtsto
                else
c           
c                   W: Set line weight.
c           
                    call readiq(input(istart:nin),1,
     &                          iwt,idum,idum,idum,*1001)
c
                end if
c
                iwtsto = iweight
                iweight = iwt
                call weight(iweight)
c
            else if (c2.eq.'i') then
c
c               WI(N): Specify which X window to use.
c
                if (idev.ne.17) then
                    if (num_win(17).le.0) then
                        if (iprompt.ne.0) write(6,'(a)')'No X windows!'
                    else
                        if (istart.gt.nin) then
                            if (iprompt.ne.0)
     $                              write(6,'(a)')'X not selected.'
                        else
c
c                           Select an existing X window.
c
                            call readiq(input(istart:nin),1,
     &                                  iwin,idum,idum,idum,*1001)
c
                            no_new_xwin = iwin
                            call setup('x','win')
                            no_new_xwin = -1
                            call set_win(iwin, ierr)
                        end if
                    end if
                else
                    if (istart.gt.nin) then
c
c                       Just print the current window ID.
c
                        if (iprompt.ne.0) then
                            iwin = curr_win()
                            if (iwin.ge.0) then
                                write(6,'(a,i4)')
     $                             'Current X output window is #',iwin
                            else
                                write(6,'(a)')'No window open.'
                            end if
                        end if
                    else
                        call readiq(input(istart:nin),1,
     &                              iwin,idum,idum,idum,*1001)
c
c                       NOTE:  Context switching must be performed here if
c                       we want different windows/devices to be independent.
c
                        icurrwin = curr_win()
                        if (iwin.ne.icurrwin) then
                            call save_context(idev,icurrwin)
                            call set_win(iwin, iret)
                            if (iret.ge.0)
     $                              call load_context(idev,iwin)
                        else
                            call set_win(iwin, iret)
                        end if
                    end if
                end if
            end if
c           
        else if (c1.eq.'x') then
c           
            if (c2.eq.' ') then
c               
c               X: Read the x array.
c               
                if (iopen.eq.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No file open'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c               
                call readiq(input(istart:nin),1,
     &                      ix,idum,idum,idum,*1001)
                if (ix.le.0)go to 1001
c               
                call readcols(10,nhead,nrange,ix,ix,x,y,z,nx,nmax,
     &                  c1,0,iprompt,*1001)
c               
            else if (c2.eq.'l') then
c               
c               XL: Specify x label.
c               
                xttl=' '
                if (istart.le.nin)xttl=input(istart:nin)
                if (c3.eq.'p'.and.ibox.eq.1) call labels(xttl,' ')
c               
            else if (c2.eq.'o') then
c               
c               XO: Set x-axis mode (XON, XOFF).
c               
                if (c3.eq.'f') then
                    call setxaxis(0)
                else if (c3.eq.'n') then
                    call setxaxis(1)
                end if
            else if (c2.eq.'w') then
c
c               XW: Print current X-window.
c
                if (idev.ne.17) then
                    if (iprompt.ne.0) write(6,'(a)')'Not using X!'
                else if (iprompt.ne.0) then
                    iwin = curr_win()
                    if (iwin.ge.0) then
                        write(6,'(a,i4)')
     $                          'Current X output window is #',iwin
                    else
                        write(6,'(a)')'No window open.'
                    end if
                end if
c
            else if (c2.eq.'x') then
c               
c               XX: Read the x array, without worrying about columns.
c               
                if (iopen.eq.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No file open'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c               
                call readall(10,input,istart,nin,x,nx,nhead,nmax,
     $                       iprompt,*1001)
c               
            end if
c           
        else if (c1.eq.'y') then
c           
            if (c2.eq.' ') then
c               
c               Y: Read the y array.
c               
                if (iopen.eq.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No file open'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c               
                call readiq(input(istart:nin),1,
     &                      iy,idum,idum,idum,*1001)
                if (iy.le.0)go to 1001
c               
                call readcols(10,nhead,nrange,iy,iy,x,y,z,ny,nmax,
     &                  c1,0,iprompt,*1001)
c               
            else if (c2.eq.'y') then
c               
c               YY: Read the y array, without worrying about columns.
c               
                if (iopen.eq.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No file open'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c               
                call readall(10,input,istart,nin,y,ny,nhead,nmax,
     $                       iprompt,*1001)
c               
            else if (c2.eq.'h') then
c
c               YH: Toggle y label size following x.
c
                call getyfollowsx(iy)
                call setyfollowsx(1 - iy)
c
            else if (c2.eq.'l') then
c               
c               YL: Specify y label.
c               
                yttl=' '
                if (istart.le.nin)yttl=input(istart:nin)
                if (c3.eq.'p'.and.ibox.eq.1) call labels(' ',yttl)
c               
            else if (c2.eq.'a') then
c               
c               YA: Set y-axis mode.
c               
                call readiq(input(istart:nin),1,
     &                      ly,idum,idum,idum,*1001)
                call setyaxis(ly)
c               
            else if (c2.eq.'o') then
c               
c               YO: Set y-axis (YON, YOFF).
c               
                if (c3.eq.'f') then
                    call setyaxis(-1)
                else if (c3.eq.'n') then
                    call setyaxis(0)
                end if
c               
            end if
c           
        else if (c1.eq.'z') then
c           
            if (c2.eq.'s') then
c               
c               ZS: Display the zoom stack.
c               
                if (ireplay.eq.0.and.iprompt.eq.1) then
                    do 600 iz=1,nzoom
                        write(6,*)iz,xlzoom(iz),xrzoom(iz),
     &                          ybzoom(iz),ytzoom(iz)
600                 continue
                end if
c               
            else if (c2.eq.'o'.or.c2.eq.'+'.or.c2.eq.'-'.or.c2.eq.'0')
     &                  then
c               
c               ZO(+,-,0): Zoom in/out on a particular plot region.
c               
                if (ibox.eq.0) go to 1001
c               
c               Work with izoom, rather than nzoom, throughout.
c               
                if (ireplay.eq.0) izoom = nzoom
c               
                izm = 0
                if (c2.eq.'-'.or.c2.eq.'0') izm = 1
c               
                if (izm.eq.0.and.izoom.ge.NZMAX) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'Zoom stack full !'
                    end if
                    go to 1001
                else if (izm.eq.1.and.izoom.le.0) then
                    go to 1001
                end if
c               
                if (izm.eq.0) then
                    if (ireplay.eq.0) then
c                       
c                       Get graphic input.
c                       
                        if (iprompt.eq.1) then
                            call devoff
                            write(6,*)
     &                          'Indicate bottom left, ',
     &                          'top right corners'
                        end if
                    end if
c                   
                    call getgfx(r1,s1)
                    call getgfx(r2,s2)
c                   
                    izoom = izoom + 1
c                   
                    if (ireplay.eq.0) then
c                       
c                       Save old limits on the zoom stack.
c                       
                        nzoom = nzoom + 1
c                       
                        xlzoom(nzoom) = xmin
                        xrzoom(nzoom) = xmax
                        ybzoom(nzoom) = ymin
                        ytzoom(nzoom) = ymax
                    end if
c                   
c                   Convert new limits to user units.
c                   
                    call fr users(r1,s1,xz1,yz1)
                    call fr users(r2,s2,xz2,yz2)
c                   
                    xmin = min(xz1,xz2)
                    ymin = min(yz1,yz2)
                    xmax = max(xz1,xz2)
                    ymax = max(yz1,yz2)
c                   
                    if (modex.lt.0) then
                        xmin = 10.**xmin
                        xmax = 10.**xmax
                    end if
c                   
                    if (modey.lt.0) then
                        ymin = 10.**ymin
                        ymax = 10.**ymax
                    end if
c                   
                else
c                   
c                   Retrieve from the zoom stack.
c                   
                    if (c2.eq.'-'.and.istart.le.nin) then
                        call readiq(input(istart:nin),1,
     &                              ndown,idum,idum,idum,*1001)
c
                        if (ndown.gt.izoom) then
                            if (iprompt.eq.1)
     &                          write(6,*)'Zoom stack underflow'
                            ndown = izoom
                        end if
                    else
                        if (c2.eq.'-') then
                            ndown = 1
                        else
                            ndown = izoom
                        end if
                    end if
c                   
                    izoom = izoom - ndown + 1
c                   
                    xmin = xlzoom(izoom)
                    xmax = xrzoom(izoom)
                    ymin = ybzoom(izoom)
                    ymax = ytzoom(izoom)
c                   
                    izoom = izoom - 1
                    if (ireplay.eq.0) nzoom = izoom
                end if
c               
c               Erase the appropriate segment of the screen.
c               
                if (jbox.eq.0) then
                    call clear
                else
                    call erase(-roff,.5*xsize-roff,
     &                         -soff,.5*ysize-soff)
                end if
c               
c               Draw the new box.
c               
                call eframe(xmin,xmax,xlen,modex,xttl,
     &                  ymin,ymax,ylen,modey,yttl)
c               
                if (iprompt.eq.1) then
                    write(6,*)'Zoom level = ',izoom
                    call devoff
                end if
c               
            else if (c2.eq.'z') then
c               
c               ZZ: Read the z array, without worrying about columns.
c               
                if (iopen.eq.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No file open'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c               
                call readall(10,input,istart,nin,z,nz,nhead,nmax,
     $                       iprompt,*1001)
c               
            else
c               
c               Z: Read the z array.
c               
                if (iopen.eq.0) then
                    if (iprompt.eq.1) then
                        call devoff
                        write(6,*)'No file open'
                        go to 10
                    else
                        go to 1001
                    end if
                end if
c               
                call readiq(input(istart:nin),1,
     &                      iz,idum,idum,idum,*1001)
                if (iz.le.0)go to 1001
c               
                call readcols(10,nhead,nrange,iz,iz,x,y,z,nz,nmax,
     &                  c1,0,iprompt,*1001)
c               
            end if
        end if
c       
        go to 10
c
c       (Very) rudimentary error handling:
c
1001    if (inunit.eq.5.and.nino.eq.nl) then
            if (iprompt.eq.1) then
                if (audio.and.idev.eq.15) call message
                call devoff
                write(6,*)'Input error'
            end if
        else
            if (iprompt.eq.1) then
                if (audio.and.idev.eq.15) call message
                call devoff
                write(6,*)'Input error ('//
     &                  input(1:max(1,nino-1))//')'
            end if
        end if
c
c       Terminate looping and force a new line to be read:
c
        if (loop.ge.0) then
            loop = -1
            iscolon = nl + 1
        end if
c
        go to 10
c       
99999   end
        
        
        subroutine modifyh(hlast,h,c3,input,*)
        save
        character*1 c3
        character*(*) input
c       
c       Modify h according to the contents of c3 and inout.
c       
c       
        if (c3.eq.'-') then
            h = hlast
        else if (c3.eq.'*') then
            hlast = h
            call readrq(input,1,
     &                  fac,dum,dum,dum,*1001)
            if (fac.ge.0.) h = h * fac
        else if (c3.eq.'/') then
            hlast = h
            call readrq(input,1,
     &                  fac,dum,dum,dum,*1001)
            if (fac.gt.0.) h = h / fac
        else
            hlast = h
            call readrq(input,1,
     &                  h,dum,dum,dum,*1001)
        end if
c       
        return
1001    return 1
c       
        end
        
        
        subroutine message
        save
c       
c       Output an audio message, if possible.
c       
        character*80 value,file
c       
c       call mygetenv('SHELL',value)
c       if (index(value,'/bin/csh').eq.0) return
c       
        call mygetenv('WINDOW_PARENT',value)
        if (index(value,'/dev/win').eq.0) return
c       
        file = '/home/zonker_export/steve/sorrydave.au'
        idum = mysystem(
     &          '(((cat '//file//' > /dev/audio) >& /dev/null )&)')
c       
        write(6,*)
     &      '(To turn off the annoying audio message, type "h-")'
c       
        end


        subroutine readrq(string,n,a1,a2,a3,a4, *)
c
c       Read a specified number of real arguments from a string.
c
        character*(*) string
        real a1,a2,a3,a4
        character*50 temp(20)
        common /read_token_stat/ io
c
        call gettokens(string,temp,nt)
c
        if (nt.eq.1) then
            call readrtoken(temp(1),a1,a1)
        else if (nt.eq.2) then
            call readrtoken(temp(1),a1,a1)
            if (io.ne.0) return 1
            call readrtoken(temp(2),a2,a2)
        else if (nt.eq.3) then
            call readrtoken(temp(1),a1,a1)
            if (io.ne.0) return 1
            call readrtoken(temp(2),a2,a2)
            if (io.ne.0) return 1
            call readrtoken(temp(3),a3,a3)
        else if (nt.eq.4) then
            call readrtoken(temp(1),a1,a1)
            if (io.ne.0) return 1
            call readrtoken(temp(2),a2,a2)
            if (io.ne.0) return 1
            call readrtoken(temp(3),a3,a3)
            if (io.ne.0) return 1
            call readrtoken(temp(4),a4,a4)
        end if
        if (io.ne.0) return 1
c
        end


        subroutine readiq(string,n,a1,a2,a3,a4,*)
c
c       Read a specified number of integer arguments from a string.
c
        character*(*) string
        integer a1,a2,a3,a4
        character*50 temp(20)
        common /read_token_stat/ io
c
        call gettokens(string,temp,nt)
c
        if (nt.eq.1) then
            call readitoken(temp(1),a1,a1)
        else if (nt.eq.2) then
            call readitoken(temp(1),a1,a1)
            if (io.ne.0) return 1
            call readitoken(temp(2),a2,a2)
        else if (nt.eq.3) then
            call readitoken(temp(1),a1,a1)
            if (io.ne.0) return 1
            call readitoken(temp(2),a2,a2)
            if (io.ne.0) return 1
            call readitoken(temp(3),a3,a3)
        else if (nt.eq.4) then
            call readitoken(temp(1),a1,a1)
            if (io.ne.0) return 1
            call readitoken(temp(2),a2,a2)
            if (io.ne.0) return 1
            call readitoken(temp(3),a3,a3)
            if (io.ne.0) return 1
            call readitoken(temp(4),a4,a4)
        end if
        if (io.ne.0) return 1
c
        end


      subroutine plain_sim(string, nstring)
c
c     Replace numbers in string by "plain-font" versions.
c
      character*(*) string
      integer nstring
c
      if (nstring.le.0) return
c
      ndigit = 0
      do i=1,nstring
          if (string(i:i).ge.'0'.and.string(i:i).le.'9'
     $            .and.(i.eq.1.or.(string(i-1:i-1).ne.'@'
     $                             .and.string(i-1:i-1).ne.'%')))
     $            ndigit = ndigit + 1
      end do
c
c     String is pure numeric.  Replace digit by "@digit".
c
      nstring = nstring+ndigit
      i = nstring
      do while (i.gt.0)
          string(i:i) = string(i-ndigit:i-ndigit)
          ind = i-ndigit-1
          if (string(i:i).ge.'0'.and.string(i:i).le.'9'
     $            .and.(i.eq.1.or.(string(ind:ind).ne.'@'
     $                             .and.string(ind:ind).ne.'%'))) then
              i = i - 1
              string(i:i) = '@'
              ndigit = ndigit - 1
          end if
          i = i - 1
      end do
c
      if (ndigit.ne.0) write(6,*)'Error in plain_sim!!!!!!!'
c
      end
