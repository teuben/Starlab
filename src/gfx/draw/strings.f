
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

c--------------------------------------------------------------------------
c
c       Contents:  chkhead    - Check for leading character in string
c                  chktail    - Check for trailing character in string
c                  chksubs    - Check for substitutions in string
c                  stripbl    - Remove non-significant blanks
c                  cleanup    - Beautify a string 
c                  shiftstr   - Left shift a string.
c                  locsubstr  - Locate a substring with specific delimiters.
c                  locchar    - Locate a character in a string
c                  repsubstr  - Replace a substring by location
c                  substitute - Repeatedly replace a substring by name.
c
c--------------------------------------------------------------------------

	subroutine chkhead(line,nl,char,iflag,istrip)
c
c       Check for presence of a specified leading character
c
c       Set iflag = 1 if the character is found, and optionally
c       strip the character.
c
	character*(*) line
        character*1 char
c
        iflag = 0
c
c       Search for the character.
c
        do i=1,nl
            if (line(i:i).eq.char) iflag = 1
            if (line(i:i).gt.' ') go to 100
        end do
	i = nl + 1
c
c       Strip the character, if desired.
c
100     if (iflag.eq.1.and.istrip.ne.0) then
            do ii = i+1,nl
                line(ii-i:ii-i) = line(ii:ii)
            end do
            nl = nl - i
        end if
c
	end


	subroutine chktail(line,nl,char,iflag,istrip)
c
c       Check for presence of a specified trailing character
c
c       Set iflag = 1 if the character is found, and optionally
c       strip the character.
c
	character*(*) line
        character*1 char
c
        iflag = 0
c
c       Search for the character.
c
        do i = nl,1,-1
            if (line(i:i).eq.char) iflag = 1
            if (line(i:i).gt.' ') go to 100
        end do
c
c       Strip the trailing character, if requested.
c
100     if (iflag.eq.1.and.istrip.ne.0) nl = i-1
c
	end

        
	subroutine chksubs(line,nl,nhist,iexpand,subso,no,subsn,nn)
c       
c       Check for substitutions (.....^xxx"yyy" substitutes yyy for xxx).
c	Return the old and new strings in the arrays subso and subsn, and
c       truncate the line.  Retain the trailing '"' in the convention to
c       allow blanks in the substitute string.  Each invocation of this
c       routine will lead to the final ^xxx"yyy" substitution being flagged.
c
c       On return, no > 0 if a substitution is due.
c
c       The substitutions are not made here (history expansion must occur
c       first), but the ^xxx"yyy" piece is removed from the string.
c       If only the "substitute" piece is found, insert "!!" as the string.
c       
	character*(*) line,subso,subsn
c
        isubs = 0
        iexpand = 0
        no = 0
        nn = 0
c
        i = nl
        do while (i.gt.0) 
            if (line(i:i).eq.'"') then
                if (line(i-1:i-1).ne.'"') then
                    isubs = isubs+1
                    if (isubs.eq.1) then
                        icar = i
                    else if (isubs.gt.2) then
                        i = 0
                    else
                        nn = icar - 1 - i
                        subsn(1:nn) = line(i+1:icar-1)
                        icar = i
                    end if
                else
                    i = i - 1
                end if
            else if (line(i:i).eq.'^') then
                if (isubs.ge.2) then
                    if (i.gt.1.and.line(i-1:i-1).eq.'^') then
                        i = i - 1
                    else
                        no = icar - 1 - i
                        subso(1:no) = line(i+1:icar-1)
                        nl = i - 1
                        iexpand = 1
                        i = 0
                        if (nl.le.0) then
                            write(line(1:6),'(''!'',i5)')nhist
                            nl = 6
                        end if
                    end if
                end if
            else if (isubs.eq.0.and.line(i:i).ne.' ') then
                i = 0
            end if
            i = i - 1
        end do
c
	end


	subroutine stripbl(line,nl,*,*)
c
c	Strip off trailing blanks and non-significant semicolons
c       from the input line.
c
	character*(*) line
c
        ns = 0
        do i = nl,1,-1
            if (line(i:i).ne.' ') then
                if (line(i:i).eq.';') then
                    ns = 1
                    do j = i-1,1,-1
                        if (line(j:j).ne.';') go to 100
                        ns = ns+1
                    end do
                    return 1
                end if
                go to 100
            end if
        end do
        return 2
c
100     if (2*(ns/2).eq.ns) then
            nl = i
        else
            nl = i - 1
        end if
c
	end


	subroutine cleanup(input,nin,istart,*)
c
c	Beautify the command string.
c
	character*(*) input
c
c       Strip leading blanks.
c
        do i1 = 1,nin
            if (input(i1:i1).gt.' ') go to 50
        end do
        return 1
c
50      input(1:nin-i1+1) = input(i1:nin)
        nin = nin-i1+1
c
c	Locate the first blank.
c
        ib = 0
        nn = 0
        nino = nin
        do i = 1,nin
            if (input(i:i).eq.' ') then
                if (nn.gt.0.and.ib.eq.0) then
                    ib = i
                    go to 100
                end if
            else
                nn = 1
            end if
        end do
c
c       Add a trailing blank if there are no others.
c
100     if (ib.eq.0) then
            ib = nin+1
            input(ib:ib) = ' '
            nin = ib
        end if
c
c       Convert to lowercase.
c
        do i = 1,ib-1
            if (input(i:i).ge.'A'.and.input(i:i).le.'Z')
     &              input(i:i) = char(ichar(input(i:i))+32)
        end do
c
c       Locate the start of the argument list.
c
        do i = ib+1,nin
            if (input(i:i).ne.' ') then
                istart = i
                return
            end if
        end do
        istart = ib+1
c
	end


        subroutine shiftstr(string,n,ishift)
c
c       Shift the string left by the specified amount.
c
        character*(*) string
c       
        if (ishift.le.0) return
c
        do i=1,n-ishift
            string(i:i) = string(i+ishift:i+ishift)
        end do
        n = n - ishift
c       
        end


        subroutine locsubstr(string,n,c1,c2,i1,i2,iend)
c
c       Locate the substring of string delimited by the characters
c       c1 and c2 (skip double characters), beginning the search at
c       location i2 + 1.  The end of the string is regarded as a
c       delimiter of type c2 if iend is nonzero.
c
        character*(*) string
        character*1 c1,c2
c
        i1 = 0
        if (i2.lt.0) return
c
        i = i2 + 1
        call locchar(string,n,i,c1)
        if (i.ge.n) return
c
        i1save = i
c
        i = i + 1
        call locchar(string,n,i,c2)
        if (i.gt.n.and.iend.eq.0) return
c
        i1 = i1save
        i2 = i
c
        end


        subroutine locchar(string,n,i,c)
c
c       Find the next location of the single character c in string,
c       beginning the search at location i
c
        character*(*) string
        character*1 c
c
10      do while(i.le.n.and.string(i:i).ne.c)
            i = i + 1
        end do
c
        if (i.ge.n) return
c
        if (string(i+1:i+1).eq.c) then
            i = i + 2
            go to 10
        end if
c
        end


        subroutine repsubstr(string,n,i1in,i2in,substr,nsub)
c
c       Replace the portion of string between i1in and i2in (inclusive)
c       by substr, and adjust the string length n accordingly.
c
c       Handle the special case of insertion at the start of the
c       string by i2in = 0.
c
        character*(*) string,substr
c
        i1 = i1in
        i2 = i2in
c
c       Make room for the new string.
c
        if (i2.le.0) then
            i1 = 1
            i2 = 0
            joff = nsub
        else
            joff = nsub - (i2 - i1 + 1)
        end if
c
        if (joff.ne.0) then
            if (joff.gt.0) then
                j1 = n
                j2 = i2 + 1
                jinc = -1
            else
                j1 = i2 + 1
                j2 = n
                jinc = 1
            end if
c
            do j=j1,j2,jinc
                string(j+joff:j+joff) = string(j:j)
            end do
c
        end if
c
        string(i1:i1+nsub-1) = substr(1:nsub)
        n = n + joff
c
        end

        
        subroutine substitute(line,nl,old,no,new,nn)
c       
c       Repeatedly substitute new for old in the input line.
c
        character*(*)line,old,new
c       
        idel = no - nn
        if (idel.ne.0) inc = sign(1,idel)
c
        i1 = 1
100     if (i1.gt.nl) return
c
        i2 = i1 - 1 + index(line(i1:nl),old(1:no))
        if (i2.lt.i1) return
c
        if (inc.gt.0) then
            if = i2 + no
            il = nl
        else
            if = nl
            il = i2 + no
        end if
        if (idel.ne.0) then
            do  i = if,il,inc
                line(i-idel:i-idel) = line(i:i)
            end do
        end if
        line(i2:i2+nn-1) = new(1:nn)
        nl = nl - idel
        i1 = i2 + nn
c
c       Continue the search.
c
        go to 100
c       
        end
