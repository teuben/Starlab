
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

c------------------------------------------------------------------------
c     
c     Simple on-line help routines.
c     
c     Contents:    help
c     wr
c     fullhelp
c     
c------------------------------------------------------------------------

      subroutine help(keyword)
c     
      character*(*) keyword
c     
c     Describe available options.
c     
      call wrinit(keyword)
c     
      call wr('!{range}: reexecute {historical }command(s)')
      call wr('!n^x"y": {reexecute n and }subs{titute} y for x')
      call wr('?: {type }short help{ info}')
      call wr('={#}: specify box{ location and size}')
      call wr('*: suppress most output')
      call wr('2d{ array}: 2-D plot{ of the data in the'//
     &        ' specified array}')
      call wr('3d{ array}: 3-D plot{ of the data in the'//
     &        ' specified array}')
      call wr('a+(-*/^){ array1 array2 [array3]}: array1{ }'//
     &        '+{ }array2 (etc){ --> array2 [array3]}')
      call wr('a>{ array1 array2 [array3]}: max(array1,{ }'//
     &        'array2){ --> array2 [array3]}')
      call wr('a<{ array1 array2 [array3]}: min(array1,{ }'//
     &        'array2){ --> array2 [array3]}')
      call wr('a={ array1 array2}: array assignment'//
     &        '{ (array2 = array1)}')
      call wr('a\\{ array}: reverse array{ [array(1) <-->'//
     &        ' array(n), etc.]}')
      call wr('a!{ array nskip}: array compress'//
     &        '{ (select every (nskip+1)-th element)')
      call wr('aa{ array}: {array --> }abs(array)')
      call wr('ac{ array min max}: array clip'//
     &        '{ to min and max}')
      call wr('ai{ array}: invert array{ [array(i) --> 1/array(i)]}')
      call wr('an{ angle}: specify angle{ for strings}')
      call wr('an-: restore{ previous} angle{ for strings}')
      call wr('ao{ array offset}: offset array{ by specified'//
     &        ' amount [array(i+offset) = array(i)]}')
      call wr('ar{ xt yt xh yh}: draw an arrow{ from (xt,yt)'//
     &        ' to (xh,yh)}')
      call wr('arh{ size}: specify arrowhead size{ in screen units}')
      call wr('as{ ratio}: define aspect ratio')
      call wr('au{ array}: autocorrelation{ of specified array}')
      call wr('b: draw box{ (all parameters settable)}')
      call wr('b-{ n flag}: restore earlier box{ parameters}')
      call wr('bb: draw unadorned box')
      call wr('bc: {toggle attempt to }make "=n" frames closer')
      call wr('bi: increase box size')
      call wr('bn: scale, without {actually drawing }a box')
      call wr('bo{ x y scale}: offset{ and scale} box')
      call wr('bs: list box stack')
      call wr('c{ col#1 col#2}: input columns{ to x and y}')
      call wr('cb{ color}: set background{ color,'//
     &        ' if applicable}')
      call wr('cd{ directory}: change{ working} directory')
      call wr('cm{ file}: show {(or specify) }colormap'//
     &        '{ (if filename supplied)}')
      call wr('co{ color}: specify color{ of lines and text}')
      call wr('co-: restore{ previous} color{ of lines and text}')
      call wr('cu{ array}: {replace array by its }cumulative sum')
      call wr('cx(y,z){ filename}: read image{ file contining x(y,z)}')
      call wr('d {x y}: draw line{ to (x y)}')
      call wr('da {x y}: draw line{ to absolute (x y)}')
      call wr('dd: draw dashed line{ (set pattern with "t")}')
      call wr('de{ device}: specify {graphics }device{ or new window}')
      call wr('diff{ a1 a2 a3}: differentiate array{ a2 with'//
     &        ' respect to a1, place in a3 (def: a2)}')
      call wr('dx(y,z): {x(y,z) --> }10**x(y,z)')
      call wr('e[a]: erase {entire }screen')
      call wr('eb: erase {interior of }box')
      call wr('ec{ string}: echo input')
      call wr('el: erase label')
      call wr('er{ [xmin xmax ymin ymax]}: erase rectangle'//
     &        '{ [no argument ==> erase all]}')
      call wr('erc{ size}: {specify }error-bar cap{ size}')
      call wr('erp{ direction #sides}: {as for }"err", {but }limited'//
     $        ' by {the current }point size}')
      call wr('err{ direction #sides}: draw error bars{ from z}')
      call wr('ex(y,z): {x(y,z) --> }exp(x[y,z])')
      call wr('f{ filename [\\nh]}: specify filename'//
     &        '{ [and header length] for x, y and z input}')
      call wr('fc: close{ currently open data} file')
      call wr('f1[p]{ [xl xu]}: {linear }fit [+ plot]'//
     &        '{ in [xl, xu]}')
      call wr('f2[p]{ [xl xu]}: {linear }fit [+ plot]'//
     &        '{ in [xl, xu], weighted by z}')
      call wr('fo: toggle{ between plain and fancy}'//
     &        ' numeric labels')
      call wr('fp[p]{ m [xl xu]}: {m-th order polynomial }fit'//
     &        ' [+ plot]{ in [xl, xu]}')
      call wr('ft[p]{ m [xl xu]}: {m-th order trigonometric '//
     &        'polynomial }fit [+ plot]{ in [xl, xu]}')
      call wr('fz[p]{ m [xl xu]}: {m-th order polynomial }fit'//
     &        ' [+ plot]{ in [xl, xu], weighted by z}')
      call wr('g: get graphic input{ (and save location)}')
      call wr('gc: get color{ from color bar}')
      call wr('gd(m): g and draw (move)')
      call wr('gs{ prompt}: {prompt and }get string from stdin')
      call wr('h: {type }long help{ information}')
      call wr('h?: report {current }height{ setting}s{ (hn, hs, hp)}')
      call wr('hg{[n] del [ref]}: draw{ [normalized]} z histogram')
      call wr('hg1{[n] del [ref]}: draw{ another} z histogram')
      call wr('hg2: draw{ another} z histogram')
      call wr('hgxy: draw histogram{ from x and y}')
      call wr('hge: toggle{/set histogram} error-bars')
      call wr('hgm: toggle{/set} histo{gram} mode')
      call wr('hgs: toggle{/set} histo{gram} save{ mode}')
      call wr('hh: {type }complete help{ information}')
      call wr('hi{ range}: list {range of }history')
      call wr('hk{ string}: keyword help{ on the specified string}')
      call wr('hl: horizontal y-label')
      call wr('hn{ height}: set number heights')
      call wr('hn*{ factor}: increase hn{ by specified factor}')
      call wr('hn/{ factor}: reduce hn{ by specified factor}')
      call wr('hn-: restore{ previous} hn')
      call wr('hp{ height}: set point heights')
      call wr('hp*{ factor}: increase hp{ by specified factor}')
      call wr('hp/{ factor}: reduce hp{ by specified factor}')
      call wr('hp-: restore{ previous} hp')
      call wr('hs{ height}: set symbol heights')
      call wr('hs*{ factor}: increase hs{ by specified factor}')
      call wr('hs/{ factor}: reduce hs{ by specified factor}')
      call wr('hs-: restore{ previous} hs')
      call wr('i{ [array1]}: info on array(s)'//
     &        '{ [no argument ==> x and y]}')
      call wr('ic{ [comment-char]}: {toggle }comments in data '//
     &        'file{ [and specify comment character]}')
      call wr('id: idle{ (wait for graphic input)}')
      call wr('im: toggle input mode{ to allow non-numeric columns}')
      call wr('in{ filename}: take input from {specified }file')
      call wr('int{ a1 a2 a3}: integrate array{ a2 with respect'//
     &        ' to a1, place in a3 (def = a2)}')
      call wr('iw{ window}: iconify{ specified output} window')
      call wr('j{ jth}: set jth')
      call wr('ki{ window}: kill {specified X-}window')
      call wr('kw{ window}: kill {specified X-}window')
      call wr('l{ [xmin xmax ymin ymax]}: get x,y limits'//
     &        '{ [arguments ==> limits forced]}')
      call wr('l=: {choose }box from {loop }counter')
      call wr('la{ label}: draw overall label')
      call wr('lc: set color from loop{ counter}')
      call wr('lgx(y,z): {x(y,z) --> }log10(x[y,z])')
      call wr('ln: set n-gons from loop{ counter}')
      call wr('lnx(y,z): {x(y,z) --> }ln(x[y,z])')
      call wr('lo{ offset}: {specify }label offset'//
     &        '{ above top of box}')
      call wr('loop{ lower upper inc}: loop {as specified }over'//
     &        ' {the remaining }commands')
      call wr('lp: print loop{ counter}')
      call wr('lw: set weight from loop{ counter}')
      call wr('ly: read y from loop{ (counter gives column)}')
      call wr('m{ x y}: move pointer{ to (x,y)}')
      call wr('ma{ x y}: move pointer{ to absolute (x,y)}')
      call wr('mo{ modex modey}: specify {frame }modes')
      call wr('mv{ array1 array2}: move array')
      call wr('n{ #}: specify ngons')
      call wr('nc: toggle{ use of} comma {as delimiter in "smart" '//
     &        'input mode}')
      call wr('np{ x y}: plot an ngon{ at (x,y)}')
      call wr('nx{ value}: set x dimension')
      call wr('nxy{ value}: set x & y{ dimensions}')
      call wr('nxyz{ value}: set x, y, & z{ dimensions}')
      call wr('ny{ value}: set y dimension')
      call wr('nz{ value}: set z dimension')
      call wr('o{ frx fry}: {specify }strpos offsets')
      call wr('o-: restore{ previous strpos} offsets')
      call wr('of{ dx dy dz}: set offsets{ for x, y, z data}')
      call wr('p{ i1 i2}: plot graph{ from i1 to i2 (line '//
     &        'type set by "t")}')
      call wr('pa: new page')
      call wr('pc{ i1 i2}: plot graph, z-colored')
      call wr('pf{ color}: fill polygons{ with specified color}')
      call wr('po: plot points only{ (same as "j -1")}')
      call wr('pr: PostScript print{ and reinitialize output file}')
      call wr('prc: PostScript print{ and close current output file}')
      call wr('ps: power spectrum of y(x){, replace x by frequency}')
      call wr('psc: {toggle use of }color {in }Postscript')
      call wr('psq: close {current }Postscript {output }file')
      call wr('psr: toggle {forced }rewriting of {long }Postscript'//
     $        '{ files}')
      call wr('pwd: print {working }directory')
      call wr('px{ x1 x2}: plot graph{ from x1 to x2}')
      call wr('pz {z1 z2}: plot, subject to z'//
     &        '{ (z1 .le. z .le. z2)}')
      call wr('q: quit')
      call wr('q1: wait {(e.g. for graphic input) }and quit')
      call wr('qc: terminate graphics{, but don''t quit}')
      call wr('quiet: suppress {most }output')
      call wr('ra{ array}: assign ramp{ to array}')
      call wr('reb{ array interval}: rebin array'//
     &        '{ to specified interval}')
      call wr('rep{ range}: replay {specified }command range')
      call wr('res: reset graphics')
      call wr('rn{ array}: {assign }random numbers{ in [0,1) to array}')
      call wr('rp: register print')
      call wr('rr{ line1 line2}: specify read range{ for file input}')
      call wr('rs{ value}: register set')
      call wr('rx{[y,z] i}: save x{[y,z](i)} in register ')
      call wr('s[s]{ string}: draw a {[X/PostScript/SUN]}string')
      call wr('s+(-*/^){ array1 scalar [array2]}: scalar'//
     &        ' add, etc{ [result to array2]}')
      call wr('s={ array scalar}: array set{ to scalar value}')
      call wr('sa{ filename}: save x,y and z{ data in filename}')
      call wr('sb{ ib ie fr}: set symbol-box{ parameters'//
     &        ' (border, erase, fraction)}')
      call wr('sc{ sx sy sz}: set scaling{ for x, y, z data}')
      call wr('sl{ time}: sleep{ for a specified number of seconds}')
      call wr('sm{ factor}: decrease box size'//
     &        '{ by specified amount}')
      call wr('smo{ dx iopt1 iopt2}: smooth y(x) --> z')
      call wr('so{ array}: array {ascending }sort')
      call wr('so-{ array}: array {descending }sort')
      call wr('so2{ array1 array2}: array sort{ with baggage}')
      call wr('stat: print status {of graphics system}')
      call wr('sw{ array1 array2}: swap arrays')
      call wr('sy{ command}: {execute a }UNIX command')
      call wr('t{ [i1 i2 i3 i4]}: specify line type'//
     &        '{ [no argument ==> solid]}')
      call wr('ti{ [level]}: {specify }tick-mark level')
      call wr('tx{[y,z] n}: type {first n elements of'//
     &        ' specified }array')
      call wr('u: unplot graph')
      call wr('v: invert plot mode')
      call wr('vl: vertical y-label')
      call wr('w{ width}: {specify }line widths')
      call wr('w-: restore{ previous}line widths')
      call wr('wi{ window}: {specify }output window')
      call wr('x{ col#}: load x{ from specified column}')
      call wr('xl(p){ label}: specify {(and plot) }x-label')
      call wr('xoff: omit x-axis')
      call wr('xon: draw x-axis')
      call wr('xw: {print }current window')
      call wr('xx{ [#points]}: load x{, ignoring column structure}')
      call wr('y{ col#}: load y{ from specified column}')
      call wr('ya{ option}: specify y-axis mode')
      call wr('yh:{ (toggle) height of} y{-axis} numbers'//
     $        ' follow{s} x{-axis numbers}')
      call wr('yoff: omit y-axis')
      call wr('yon: draw y-axis')
      call wr('yl(p){ label}: specify {(and plot) }y-label')
      call wr('yy{ [#points]}: load y{, ignoring column structure}')
      call wr('z{ col#}: load z{ from specified column}')
      call wr('z-{ nout}: zoom out{ by specified number of frames}')
      call wr('z0: return to top {(unzoomed) }level')
      call wr('zoom: zoom in{ on cursor-specified region}')
      call wr('zs: display{ the} zoom stack')
      call wr('zz{ [#points]}: load z{, ignoring column structure}')
c     
      call wrend
c     
      end

      
      subroutine wr(string)
      save
c     
      character*(*) string
      character in*100,out*80,key*80,keyword*80
      common /hlevel/ ihelp
      character*2 hindx
      common /hindex/ hindx
c     
      data key/' '/
c     
      if (key(1:1).gt.' ') then
          if (index(string,key(1:nk)).le.0) return
      else
          if (hindx(1:1).gt.' '.and.hindx(1:1).ne.string(1:1)) return
          if (hindx(2:2).gt.' '.and.hindx(2:2).ne.string(2:2)) return
      end if
c     
      l = len(string)
      in(1:l) = string
      loff = 0
      iex = 0
      do 10 i = 1,l
          if(iex.eq.0)in(i-loff:i-loff) = in(i:i)
          if(ihelp.eq.2)then
              if(in(i:i).eq.'{'.or.in(i:i).eq.'}')loff = loff+1
          else
              if(in(i:i).eq.'{')iex = 1
              if(iex.eq.1)loff = loff+1
              if(in(i:i).eq.'}')iex = 0
          end if
10    continue
      l = l-loff
c     
      if(ihelp.eq.2)then
          write(6,*)in(1:l)
      else
          out(ioff+1:ioff+l) = in(1:l)
          if (l.ge.26) out(ioff+26:ioff+l) = ' '
          kount = kount+1
          if(kount.eq.3)then
              write(6,*)out(1:min(78,ioff+l))
              kount = 0
              ioff = 0
              out = ' '
          else
              ioff = ioff+26
          end if
      end if
      return
c     
      entry wrinit(keyword)
c     
      key = keyword
c     
      kount = 0
      ioff = 0
      out = ' '
      write(6,*)
      if (key(1:1).gt.' ') then
          nk = len(key)
          do while (key(nk:nk).le.' ')
              nk = nk - 1
          end do
          if (nk.le.0) then
              nk = 1
              key(1:1) = ' '
          end if
      end if
      return
c     
      entry wrend
c     
      if (kount.gt.0) write(6,*)out(1:ioff+l)
      write(6,*)
c     
      end


      subroutine fullhelp
      save
c     
c     Type out the contents of the mcdraw primer...
c     
      parameter (NDIR = 10)
      character*120 directory(NDIR)
      dimension ldir(NDIR)
c     
      character*100 line
c     
      logical UNIX
      parameter (UNIX = .true.)
c     
c     Look in the same places as for the simbol fonts.
c     
      nd = NDIR
      call listdir(directory,ldir,nd,iunit)
      if (iunit.lt.0) return
c     
      do 100 idir=1,nd
          if (ldir(idir).gt.0) then
              open(iunit,file=directory(idir)(1:ldir(idir))
     &                //'/READ_ME',
     &                status='old',form='formatted',err=100)
c             
c             Print out the contents of the help file.
c             
              write(6,*)'Found help file ',
     &                directory(idir)(1:ldir(idir))//'/READ_ME'
c             
              if (UNIX) then
c                 
c                 Suppress this, because unit 0 isn't always stderr!
c                 
c                 write(0,*)'(Spawning a subshell',
c                 &                    ' to display the help file...)'
                  call system('more '//directory(idir)(1:ldir(idir))
     &                    //'/READ_ME')
                  return
              else
c                 
c                 Brain-dead systems...
c                 
30                read(iunit,'(a)',err=200,end=200)line
c                 
                  do 50 i=100,1,-1
                      if (line(i:i).gt.' ') then
                          write(6,'(a)')line(1:i)
                          go to 30
                      end if
50                continue
c                 
                  write(6,*)
                  go to 30
              end if
          end if
100   continue
c     
200   return
      end
