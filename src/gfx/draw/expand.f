
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

c----------------------------------------------------------------------
c
c       Perform (historical) expansions of the input line.
c
c       expandr   - replace !! by !nhist
c       expandh   - insert a range of history items
c       expandl   - expand historical references
c       expands   - expand historical references of the form "!string"
c       applysubs - apply substitutions, after expansion
c
c----------------------------------------------------------------------

	subroutine expandr(line,nl)
        save
c       
c	Expand repeat references within a line (!! --> !nhist).
c
	character*(*) line
c
        parameter (NHMAX = 500)
c
        character*200 history(NHMAX)
        common /histchars/ history
        common /histnums/ lhist(NHMAX),nhist
c
        if (nl.le.1) return
c
c       *** Special case -- line may start with ^..."...", meaning
c           that the previous line is to be modified.  If the first
c           nonblank character is "^", insert "!!" at the start of
c           the line.
c
        i = 1
        do while (line(i:i).le.' ')
            i = i + 1
        end do
        if (line(i:i).eq.'^') call repsubstr(line,nl,0,0,'!!',2)
c
	is = 1
	ie = 0
c
        i = 0
        do while (i.lt.nl)
            i = i+1
            if (is.eq.0) then
                if (line(i:i).ne.';') go to 100
                if (line(i+1:i+1).eq.';') then
                    i = i+1
                    go to 100
                end if
                is = 1
                ie = 0
            else
                if (line(i:i).le.' ') go to 100
                if (line(i:i).ne.'!') then
                    is = 0
                else if (ie.eq.0) then
                    ie = 1
                else
c
c                   Perform the substitution.
c
                    do j = nl,i+1,-1
                        line(j+4:j+4) = line(j:j)
                    end do
                    nl = nl + 4
                    write(line(i:i+4),'(i5)')nhist
                    i = i + 4
                    is = 0
                end if
            end if
100     end do
c
	end


	subroutine expandh(line,nl,ic0,ihmin,ihmax,lextra,iret)
        save
c
c       Expand the specified history range (ihmin to ihmax) into the line.
c       The insertion is to occur immediately after character #ic0 in
c       the line.  The total number of characters added is returned as
c       lextra.  No trailing semicolon is added at the end.
c
c       Neither ic0 nor nl are altered.  (Bookkeeping is to be performed 
c       by the calling routine.)
c
c       The status flag iret returns 1 iff an error occurs.
c
	character*(*) line
c
        parameter (NHMAX = 500)
c       
        character*200 history(NHMAX)
        common /histchars/ history
        common /histnums/ lhist(NHMAX),nhist
c
        iret = 0
        lextra = 0
c
        do ih = ihmin,ihmax
c
            jh = ih
            do while (jh.gt.NHMAX)
                jh = jh - NHMAX
            end do
c
            if (nl+lhist(jh).ge.len(line)) then
                write(6,*)'Expanded line too long.'//
     &                    ' Check for possible recursion.'
                iret = 1
                return
            end if
c
            l0 = lextra
            lextra = lextra + lhist(jh)
c
            line(ic0+l0+1:ic0+lextra) = history(jh)(1:lhist(jh))
c
c           Add a delimiting ";", except at the end.
c
            if (ih.lt.ihmax) then
                lextra = lextra + 1
                line(ic0+lextra:ic0+lextra) = ';'
            end if
c
        end do
c
	end

        
        subroutine expandl(line,nl,isuppr,*)
        save
c
c       Repeatedly expand historical references in line.  This routine
c       serves a dual function, so the substring to be expanded can be
c       of the form "!  n;" or "!  n^".  This allows use during normal
c       line expansion, and during substitutions.
c
c       The expansion is recursive unless isuppr is nonzero.
c       
        character*(*) line
c
        parameter (NHMAX = 500)
c       
        character*200 history(NHMAX),lsto
        common/histchars/history
        common/histnums/lhist(NHMAX),nhist
c       
        character*80 temp
        character*1 which
c
        if (nl.le.0) return 1
c       
c       The logic here is similar to that used in mcdraw when
c       interpreting individual commands.
c
        iscolon = 0
10      if (iscolon.ge.nl) return
c
c       Find the next substring.  If may be terminated by ";" or "^".
c       Note that if "^" is the terminator, we will have to search
c       for the next ";" before continuing.
c
c       Note that the two expected expansion patterns are "!....;"
c       and "!...^xxx"yyy";".
c
        ic0 = iscolon
        which = ';'
c
15      do ic1 = iscolon+1,nl
            which = line(ic1:ic1)
            if (which.eq.';'.or.which.eq.'^') go to 20
        end do
        ic1 = nl + 1
20      if (ic1.lt.nl.and.line(ic1+1:ic1+1).eq.which) then
            iscolon = ic1 + 1
            go to 15
        end if
c
c       The piece to be expanded runs from ic0+1 to ic1-1.
c       The terminating character is "which", at location ic1.
c
        nin = ic1 - ic0 - 1
        temp(1:nin) = line(ic0+1:ic1-1)
        iscolon = ic1
c       
c       Expand any historical references.
c       
        do i = 1,nin
            if (temp(i:i).gt.' ') then
                if (temp(i:i).ne.'!') go to 100
                if (i.ge.nin) then
                    write(6,*)'warning: possible input error...'
                    go to 100
                end if
                call rdecode(temp(i+1:nin),nhist,ihmin,ihmax,*1001)
                go to 32
            end if
        end do
c
c       Continue expanding the line.
c
        go to 100
c       
c       Save the remainder of the string before expansion.
c
32      nlsto = 0
        if (ic1.lt.nl) then
c
c           Note that the trailing ";" or "^" is saved.
c
            nlsto = nl - ic1 + 1
            lsto(1:nlsto) = line(ic1:nl)
        end if
c
c       Expand the history list (no trailing delimiter added).
c
        call expandh(line,nl,ic0,ihmin,ihmax,lextra,iret)
        if (iret.ne.0) return 1
c
c       Restore the rest of the line.  (Note: this would be a good
c       place for a consistency check, as the string lengths are not
c       forced to agree...)
c       
        nl = nl - nin + lextra
        if (nlsto.gt.0) then
            line(ic0+lextra+1:nl) = lsto(1:nlsto)
        end if
c
c       Continue recursively at ic0, or non-recursively, depending on
c       input option isuppr.
c
100     if (isuppr.eq.0.and.which.eq.';') then
            iscolon = ic0
        else
c
c           Start at the next ';'.
c
            iscolon = ic0 + lextra + 1
            if (which.eq.'^') then
                do while (iscolon.le.nl
     &                  .and.line(iscolon:iscolon).ne.';')
                    iscolon = iscolon + 1
                end do
            end if
        end if
c
        go to 10
c       
1001    return 1
        end


        subroutine expands(line,nl)
        save
c
c       Perform symbolic history expansion on line, replacing !string
c       by !ihist.
c
        character*(*)line
        character*6 temp
c
        parameter (NHMAX = 500)
c       
        character*200 history(NHMAX)
        common/histchars/history
        common/histnums/lhist(NHMAX),nhist
c
        i2 = 0
10      if (i2+1.ge.nl) return
c
        call locsubstr(line,nl,'!',';',i1,i2,1)
        if (i1.le.0) return
c
c       Strip leading and trailing blanks.
c
        i1save = i1
        i1 = i1 + 1
        do while (line(i1:i1).le.' ')
            i1 = i1 + 1
        end do
c
c       Note that i2 is either the location of a ";" or nl + 1.
c
        i2save = i2
        i2 = i2 - 1
        do while (line(i2:i2).le.' ')
            i2 = i2 - 1
        end do
c
        if (i1.le.i2) then
c
c           Now the string to be searched for lies in line(i1:i2).
c           Only try to substitute if the entire string is non-numeric.
c
            do i=i1,i2
                if (line(i:i).lt.'0'.or.line(i:i).gt.'9') go to 50
            end do
            go to 100
c
50          do ihist = nhist,max(1,nhist-NHMAX),-1
c
                jh = ihist
                do while (jh.gt.NHMAX)
                    jh = jh - NHMAX
                end do
c
                if (index(history(jh)(1:lhist(jh)),line(i1:i2)).eq.1)
     &                  then
                    write(temp,'(''!'',i5)')ihist
                    call repsubstr(line,nl,i1save,i2save-1,temp,6)
c
c                   Note that we keep the trailing ";", if one exists.
c
                    go to 100
                end if
c
            end do
        end if
c
c       Continue the search.
c
100     i2 = i2save
        go to 10
c
        end


        subroutine applysubs(line,nl,*)
        save
c
c       Apply substitutions to the input line.  Expand each subcommand
c       as necessary before making the changes.
c
        character*(*) line
        character temp*200,subsn*80,subso*80
c
        parameter (NHMAX = 500)
c
        character*200 history(NHMAX)
        common /histchars/ history
        common /histnums/ lhist(NHMAX),nhist
c
c       Locate and process each reference.
c
        i2 = 0
c
100     call locsubstr(line,nl,'!',';',i1,i2,1)
        if (i1.eq.0) return
c
        nt = i2 - i1
        temp(1:nt) = line(i1:i2-1)
c
c       Note that temp includes the leading "!" but not the trailing ";".
c
        kount = 0
        no = 1
        do while (no.gt.0.and.kount.lt.100)
c
            call chksubs(temp,nt,nhist,iexpand,subso,no,subsn,nn)
c
c           Perform expansion, non-recursively, prior to substitution.
c
            if (iexpand.gt.0) call expandl(temp,nt,1,*99999)
            if (no.gt.0) call substitute(temp,nt,subso,no,subsn,nn)
            kount = kount + 1
c
        end do
c
        call repsubstr(line,nl,i1,i2-1,temp,nt)
c
c       Continue the search.
c
        i2 = i1 + nt
        go to 100
c
99999   return 1
c
        end
