
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
c

        block data ps setup
        save
        common /ps color set/ ipscolor
	common /ps write/ iwrite
        data ipscolor/0/iwrite/0/
        end


        subroutine set ps write
        save
c
c	Toggle forced rewriting of long files.
c
        common /ps write/ iwrite

        iwrite = 1 - iwrite

        end


	subroutine ps stroke
        save
	common /ps strokes/ nstroke,nstrpage,nstroketot
c
	if (nstroke.gt.0) then
	    write(42,'(''stroke'')')
	    nstroke = 0
	end if
c
	end


        subroutine set ps color
        save
c
c       Toggle between greyscale and color output.
c
        common /ps color set/ ipscolor
        common /mcpak_colormap/ ncolors
c
        parameter (NXCOLORS = 16)
c
        ipscolor = 1 - ipscolor
c
        if (ipscolor.eq.0) then
c
c           Greyscale, 256 levels:
c
            ncolors = 256
c
        else
c
c           Color, NXCOLORS levels:
c
            ncolors = NXCOLORS
c
        end if
c
        return
c
        entry get ps color(ipsc)
        ipsc = ipscolor
c
        end


        subroutine ps color(icolor)
        save
        common /ps color set/ ipscolor
        common /mcpak_colormap/ ncolors
c
c       Set up PostScript color/grey level.
c
c       The color map mimics the current X-window colors (see "mcdxsubs.c"):
c
c           0 = "white"
c           1 = "black"
c           2 = "blue"
c           3 = "purple"
c           4 = "violet"
c           5 = "magenta"
c           6 = "light blue"
c           7 = "gray"
c           8 = "cyan"
c           9 = "green"
c           10 = "lime green"
c           11 = "yellow"
c           12 = "orange"
c           13 = "brown"
c           14 = "red"
c           15 = "pink"

c	Note that color 0 is white here, whereas for greyscale, 0 is black!

        parameter (NXCOLORS = 16)
        integer red(0:NXCOLORS-1),green(0:NXCOLORS-1),blue(0:NXCOLORS-1)
c
        data   red/ 255,   0,   0, 160, 238, 255, 191, 192,   0,   0,
     &               50, 255, 255, 165, 255, 255/
        data green/ 255,   0,   0,  32, 130,   0, 239, 192, 255, 255,
     &              205, 255, 165,  42,   0, 192/
        data  blue/ 255,   0, 255, 240, 238, 255, 255, 192, 255,   0,
     &               50,   0,   0,  42,   0, 203/
c
        call ps stroke
c
        if (ipscolor.eq.0) then
c
c           Greyscale:
c
            ic = max(0,icolor)
            do while (ic.ge.ncolors)
                ic = ic - ncolors
            end do

c           Force 0 and 1 to be black for a large color range.

            if (ncolors.gt.16 .and. ic.eq.1) ic = 0
            write(42,'(f5.3,'' setgray'')')ic/(ncolors-1.)
c
        else
c
c           Color:
c
            ic = max(0,icolor)
            do while (ic.ge.ncolors)
                ic = ic - ncolors
            end do
            write(42,'(3f6.3,'' setrgbcolor'')')
     &                   red(ic)/255.,green(ic)/255.,blue(ic)/255.
c
        end if
c
        end


	subroutine ps page(iopt)
        save
	common /ps strokes/ nstroke,nstrpage,nstroketot
	common /ps pages/ npage
c
	if (nstrpage.eq.0) return
	call printid
c	write(42,'(''gsave showpage grestore newpath'')')
	write(42,'(''showpage newpath'')')
	if(iopt.eq.0)return
c
	npage = npage+1
	write(42,'(/''%%Page: '',i4)')npage
	nstrpage = 0
c
	end
