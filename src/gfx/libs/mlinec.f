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
        subroutine mlinec(xarray,yarray,zarray,n,jth,jsymbl,htsym)
        save
c       
c	Plots the n points yarray(xarray), with identifying
c	symbols, of height |htsym|, at every |jth|-th point
c	for nonzero jth.
c       
c	This version differs from mline in that the color of the
c	plot can be continuously adjusted, according to zarray.
c       
c	if jth = 0, the points are joined by solid lines only.
c	if jth > 0, both lines and symbols are drawn.
c	if jth < 0, only the symbols are drawn.
c	
c	If jsymbl >0 or =0, the symbol is drawn by subroutine ngon.
c       
c	If jsymbl < 0, the symbol is a centered "simbol" symbol.
c	Defining iascii = |jsymbl|, the first font is obtained
c	for iascii.le.127, the second for 128.le.iascii.le.223
c	and the third for iascii.ge.224, so, for example, 58
c	becomes a colon (":"), 58+96=154 is an integral sign ("@:")
c	and 58+192=250 is a gothic "z" ("%:").
c       
c	No lines are drawn outside the box produced by "eframe"
c	if mode is nonzero.
c	If htsym < 0., the "simbol" symbols will not be centered:
c	the symbol will be drawn with the array point at the
c	bottom left-hand corner.
c       
        dimension xarray(1),yarray(1),zarray(1)
c       
        logical in0
        character sim*3
c       
        common /scales/ xminim,xmax,dxinch,
     &                  yminim,ymax,dyinch,rlen,slen
        common /fr bnds/ mode
        common /mline on/ imline
c       
        common /dash/ dpatrn(10),dpat,npatrn,ipat,lpen
c       
        common /ngon stars/ istar
        common /mcpak_colormap/ ncolors
c       
        cinch(x,x0,dxi) = (x-x0)*dxi
c       
        call routine id('mlinec')
        idline = 0
        go to 1
c       
        entry dlinec(xarray,yarray,zarray,n,jth,jsymbl,htsym)
c       
c       As for mlinec, but plot dashed lines.
c       
        call routine id('dlinec')
        idline = 1
c       
c       Reinitialize the dash pattern.
c       
        dpat = dpatrn(1)
        ipat = 1
        lpen = 2
c       
1       call minmax(zarray,n,zmin,zmax)
        if (zmax.le.zmin) then
            z0 = 1.
            zfac = 1.
        else
            z0 = .05*ncolors
            zfac = .9*(ncolors-1.)/(zmax-zmin)
        end if
c       
        if (jth.ne.0) then
c           
c           *** plot symbols ***
c           
            hite = abs(htsym)
c           
            if(jsymbl.lt.0)then
                iascii = -jsymbl
                if(iascii.le.127)then
                    nsym = -1
                    sim = char(iascii)
                else if(iascii.le.223)then
                    nsym = -2
                    sim = '@'//char(iascii-96)
                else
                    nsym = -2
                    sim = '%'//char(iascii-192)
                end if
                if(htsym.lt.0.)nsym = -nsym
            end if
c           
            int = abs(jth)
            do i=1,n,int
                x1 = cinch(xarray(i),xminim,dxinch)
                y1 = cinch(yarray(i),yminim,dyinch)
                call color(nint(z0+(zarray(i)-zmin)*zfac))
                if (mode.eq.0
     &                  .or.(x1.ge.0..and.x1.le.rlen
     &			.and.y1.ge.0..and.y1.le.slen)) then
                    if(jsymbl.ge.0)then
                        if(istar.eq.0)then
                            call ngon(x1,y1,.5*htsym,jsymbl,0.)
                        else
                            call ngon(x1,y1,.5*htsym,-jsymbl,0.)
                        end if
                    else
                        if (iascii.le.1) then
c                           
c                           Encode the number of the point in sim.
c                           
                            ii = i
                            do while (ii.gt.61)
                                ii = ii - 61
                            end do
                            if (ii.le.9) then
                                sim = char(48+ii)
                            else if (ii.le.35) then
                                sim = char(87+ii)
                            else
                                sim = char(29+ii)
                            end if
                        end if
                        call simbol(x1,y1,hite,sim,0.,nsym)
                    end if
                end if
            end do
c           
            if (jth.lt.0) return
        end if
c       
c       *** plot line ***
c       
        x0 = cinch(xarray(1),xminim,dxinch)
        y0 = cinch(yarray(1),yminim,dyinch)
        in0 = (x0.ge.0..and.x0.le.rlen
     &           .and.y0.ge.0..and.y0.le.slen)
c       
c       The following distinction is necessary to set up "plotin" internally:
c       
        if (mode.eq.0) then
            call plot(x0,y0,3)
        else
            call plotin(x0,y0,3)
            imline = 1
        end if
c       
        do 20 i=2,n
            x1 = cinch(xarray(i),xminim,dxinch)
            y1 = cinch(yarray(i),yminim,dyinch)
c           
            call color(nint(z0+(zarray(i)-zmin)*zfac))
c           
            if (mode.eq.0) then
                if (idline.eq.0) then
                    call plot(x1,y1,2)
                else
                    call dplot(x1,y1,2)
                end if
            else
c               
c               Is there anything to plot? (Problems with slow graphics devices.)
c               
                if (max(x0,x1).ge.0..and.min(x0,x1).le.rlen
     &                  .and.max(y0,y1).ge.0.
     &                  .and.min(y0,y1).le.slen) then
c                   
                    if (.not.in0) call plotin(x0,y0,3)
c                   
                    if (idline.eq.0) then
                        call plotin(x1,y1,2)
                    else
                        call dplotin(x1,y1,2)
                    end if
c                   
                    in0 = .true.
                else
                    in0 = .false.
                end if
                x0 = x1
                y0 = y1
            end if
20      continue
c       
c       Make sure all pointers are set to the end of the array.
c       
        if (mode.eq.0) then
            call plot(x1,y1,3)
        else
            call plotin(x1,y1,3)
            imline = 0
        end if
c       
        end
