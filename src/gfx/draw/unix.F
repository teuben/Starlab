
c     ***********************************************************************
c     *									    *
c     *  	UNIX.F -- LOOKALIKES FOR UNIX BUILT-INS:		    *
c     *									    *
c     *				nparams = iargc				    *
c     *				getparam = getarg			    *
c     *				fdate					    *
c     *				flush					    *
c     *				random					    *
c     *									    *
c     ***********************************************************************


#include "translate.h"


C     ***********************************************************************
C     *                                                                     *
C     *  THE FOLLOWING SHOULD BE STANDARD UNIX, BUT SOMETIMES ARENT:	    *
C     *                                                                     *
C     *         IARGC -->  NPARAMS here          	                    *
C     *         GETARG --> GETPARAM here             	                    *
C     *                                                                     *
C     ***********************************************************************


      INTEGER FUNCTION NPARAMS()
C
C     Return number of command-line arguments (excluding the command name).
c     Just return 0 for non-UNIX systems.
C
      nparams = 0

#if defined(sun) || defined(cray) || defined(linux)

      nparams = iargc()

#elif defined(hpux)

$HP9000_800 INTRINSICS ON
      nparams = iargc() - 1	/* Command name is argument #1 in HPUX. */

#endif

      END


      SUBROUTINE GETPARAM(I,ARG)
      CHARACTER*(*) ARG
C
C     Return i-th command-line argument (0 = command name).
c     Just return null for non-UNIX systems.
C
      ARG = ' '

#if defined(sun) || defined(cray) || defined(linux)

      call getarg(i,arg)

#elif defined(hpux)

$HP9000_800 INTRINSICS ON
      call getarg(i+1,arg)	/* Command name is argument #1 in HPUX. */

#endif

      END


C     ***********************************************************************
C     *                                                                     *
C     *  THE FOLLOWING APPEARS TO BE SPECIFIC TO SUN UNIX:		    *
C     *                                                                     *
C     *         FDATE             			                    *
C     *                                                                     *
C     ***********************************************************************

#ifndef sun
      CHARACTER*24 FUNCTION FDATE()
C
C     RETURN CURRENT TIME OF DAY AND DATE AS A 24-CHARACTER STRING.
C
#ifdef cray
      CHARACTER*8 TIM,DAT
#endif
#ifdef hpux
$NOSTANDARD SYSTEM
      CHARACTER TIM*8,DAT*9
#endif
#ifdef SX2
      CHARACTER*8 TIM,DAT
#endif
C
      FDATE = ' '
#ifdef cray
      CALL CLOCK(TIM)
      CALL DATE(DAT)
      FDATE = ' '//TIM//'   '//DAT
#endif
#ifdef hpux
      CALL TIME(TIM)
      CALL DATE(DAT)
      FDATE = ' '//TIM//'   '//DAT
#endif
#ifdef SX2
      CALL DATIM(DAT,TIM,1)
      FDATE = ' '//TIM//'   '//DAT
#endif
C
      END
#endif


C     ***********************************************************************
C     *                                                                     *
C     *  THE FOLLOWING EXISTS IN SUN UNIX AND UNICOS (AT LEAST):	    *
C     *                                                                     *
C     *         FLUSH             			                    *
C     *                                                                     *
C     ***********************************************************************
C
C
C
#ifndef sun
#ifndef cray
#ifndef hpux
#ifndef linux
      SUBROUTINE FLUSH(LUNIT)
C
C     Flush the buffer associated with the specified logical unit.
C     Do nothing in non-unix cases.
C
      END
#endif
#endif
#endif
#endif

#ifdef hpux
      subroutine flush(lunit)
c
c     Special treatment for unit 6 only on hpux...
c
      if (lunit.eq.6) call cflush()
      end
#endif

      subroutine sx2dum
      end


C     ***********************************************************************
C     *                                                                     *
C     *  THERE DOESNT SEEM TO BE A STANDARD VERSION OF THIS...		    *
C     *                                                                     *
C     *         RANDOM             			                    *
C     *                                                                     *
C     ***********************************************************************
C
C
C
#ifdef hpux
$HP9000_800 INTRINSICS ON
#endif

      FUNCTION RANDOM(ISEED)
C
C     Return a random number in the range [0,1).
C
C     Convention:  ISEED = 0 ==> return next number in the seqruence.
C		   ISEED > 0 ==> initialize the random generator.
C
C     Use a built-in random number generator if one exists.  Otherwise, use
C     the congruential random number generator CONGRNO below.
C
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8 RANDOM
c
#ifdef sun
      real*8 DRAND
#endif

#ifndef hpux
#ifndef sun
#ifndef cray
#ifndef vms
c
c     Declare portable generator for systems with no known built-in.
c
      real*8 congrno
#endif
#endif
#endif
#endif
c
      SAVE IRSET
      DATA IRSET/0/
C
      IF (IRSET.EQ.0.OR.ISEED.GT.0) THEN
C
C         ADDITIONAL SETUP MAY BE REQUIRED ON SOME SYSTEMS:
C
#ifdef hpux
	  call srand(iseed)
#endif
#ifdef cray
          CALL RANSET(ISEED)
#endif
c
#ifndef hpux
#ifndef cray
#ifndef sun
#ifndef vms
C
C         NO SYSTEM BUILT-IN RANDOM NUMBER GENERATOR IS KNOWN:
C
C          dum = FLOAT(ISEED)
C          CALL CONGSET(dum)
c
        dum = congrno(iseed)
C
#endif
#endif
#endif
#endif
C
          IRSET = 1
      END IF
C
C     RETURN THE NEXT RANDOM NUMBER IN THE SEQUENCE
C     ---------------------------------------------
C
#ifdef hpux
      random = rand()
#endif
#ifdef sun
      RANDOM = DRAND(ISEED)
#endif
#ifdef vms
      RANDOM = RAN(ISEED)
#endif
#ifdef cray
      RANDOM = RANF()
#endif
c
#ifndef hpux
#ifndef sun
#ifndef cray
#ifndef vms
C
C     NO SYSTEM RANDOM BUILT-IN NUMBER GENERATOR IS KNOWN:
C
C      RANDOM = CONGRNO1()
c
      RANDOM = CONGRNO(0)
C
#endif
#endif
#endif
#endif
c
      END


#ifndef hpux
#ifndef sun
#ifndef cray
#ifndef vms
        
        
        function congrno(iseed)
c       
c       Random number generator (Press et al. p. 195).
c       ---------------------------------------------
c       
        implicit real*8 (a-h,o-z)
        real*8 congrno

        PARAMETER (M=714025, IA=1366, IC=150889, RM=1.d0/M)

        integer ir(97) 
        save iy,iff,ir,idum
        data iff /0/

        if (iseed.lt.0.or.iff.eq.0) then
            idum = -abs(iseed)
            iff = 1
            idum = mod(IC-idum,M)
            do 11 j = 1,97
                idum = mod(IA*idum+IC,M)
                ir(j) = idum
11          continue
            idum = mod(IA*idum+IC,M)
            iy = idum
        end if

        j = 1 + (97*iy)/m
        if (j.gt.97.or.j.lt.1) write (6,12) j, idum
12      format (/' RAN2:  j, idum = ',2i12)

        iy = ir(j)
        congrno = iy*RM
        idum = mod(IA*idum+IC,M)
        ir(j) = idum

        end 
        
c
c       Alternative version:
c
        FUNCTION CONGRNO1()
        IMPLICIT REAL*8 (A-H, O-Z)
        REAL*8 CONGRNO1
        SAVE
C       
C       Generate a pseudo-random number in the range [0,1).  The "magic"
C       numbers (AM, B1, B2) are intended for use with a 64-bit REAL datum.
C       
        REAL*8 DUM
        COMMON /RNUMS/ X1,X2
C       
        DATA X10,X20/2*1234567.0/
        DATA AM,B1,B2/7.0368743D13,1.4348907D7,1953125.0/
C       
        X1 = B1 * X1
        X2 = B2 * X2
        X1 = MOD (X1, AM)
        X2 = MOD (X2, AM)
        Y = X1 + X2
        Y = MOD (Y, AM)
        CONGRNO1 = Y / AM
        RETURN
C       
C       Initialization:
C       
        ENTRY CONGSET(DUM)
        X1 = X10 + DUM
        X2 = X20
        RETURN
C       
C       Restoration:
C       
        ENTRY CRESET(X1O,X2O)
        X1=X1O
        X2=X2O
C       
        END
#endif
#endif
#endif
#endif
