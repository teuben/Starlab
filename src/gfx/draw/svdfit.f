
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
c	Least-squares fitting by singular-value decomposition.
c	From "Numerical Recipes," section 14.3, pp. 515-519, with
c	minor modifications.
c
c	Perform a fit to the NDATA points Y(X), with weight 1/SIG, returning
c	the MA coefficients A of the fitting function, in terms of the basis
c	functions defined by the user-supplied function FUNCS.  For now, the
c	input arrays U, V, and W, of dimension (MP,NP), (NP,NP), and NP,
c	respectively, are for workspace.  We require MP >= NDATA, NP >= MA.
c

c
c	Polynomial basis functions:
c	--------------------------
c
	subroutine poly(x,f,n)
        save
	dimension f(n)
c
	f(1) = 1.
	do 10 i=2,n
10	f(i) = x*f(i-1)
c
	return
	end

c
c	Trigonometric basis function1s:
c	-----------------------------
c
	subroutine trig(x,f,n)
        save
	dimension f(n)
c
	do 10 i=1,n,2
10	f(i) = cos(((i-1)/2)*x)
	do 20 i=2,n-1,2
20	f(i) = sin((i/2)*x)
c
	return
	end



        SUBROUTINE SVDFIT(X,Y,SIG,NDATA,A,MA,U,V,W,B,MP,NP,CHISQ,FUNCS)
        save
c
        PARAMETER(MMAX=21,TOL=1.E-5)
        DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),V(NP,NP),
     *            U(MP,NP),W(NP),B(MP),AFUNC(MMAX)
c
        external funcs
c
        DO 95000 I=1,NDATA
                CALL FUNCS(X(I),AFUNC,MA)
                TMP=1./SIG(I)
                DO 95001 J=1,MA
                        U(I,J)=AFUNC(J)*TMP
95001           CONTINUE
                B(I)=Y(I)*TMP
95000   CONTINUE
c
c	Pad U with blank rows, if necessary (see NR, p. 60)...
c
	do 100 i=ndata+1,ma
	do 100 j=1,ma
100	u(i,j) = 0.
c
        CALL SVDCMP(U,max(ma,NDATA),MA,MP,NP,W,V)
        WMAX=0.
        DO 95002 J=1,MA
                IF(W(J).GT.WMAX)WMAX=W(J)
95002   CONTINUE
        THRESH=TOL*WMAX
        DO 95003 J=1,MA
                IF(W(J).LT.THRESH)W(J)=0.
95003   CONTINUE
        CALL SVBKSB(U,W,V,NDATA,MA,MP,NP,B,A)
        CHISQ=0.
        DO 95004 I=1,NDATA
                CALL FUNCS(X(I),AFUNC,MA)
                SUM=0.
                DO 95005 J=1,MA
                        SUM=SUM+A(J)*AFUNC(J)
95005           CONTINUE
                CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
95004   CONTINUE
        RETURN
        END



        SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
        save
c
        PARAMETER (NMAX=1000)
        DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
c
        G=0.0
        SCALE=0.0
        ANORM=0.0
        DO 95000 I=1,N
                L=I+1
                RV1(I)=SCALE*G
                G=0.0
                S=0.0
                SCALE=0.0
                IF (I.LE.M) THEN
                        DO 95001 K=I,M
                                SCALE=SCALE+ABS(A(K,I))
95001                   CONTINUE
                        IF (SCALE.NE.0.0) THEN
                                DO 95002 K=I,M
                                        A(K,I)=A(K,I)/SCALE
                                        S=S+A(K,I)*A(K,I)
95002                           CONTINUE
                                F=A(I,I)
                                G=-SIGN(SQRT(S),F)
                                H=F*G-S
                                A(I,I)=F-G
                                IF (I.NE.N) THEN
                                        DO 95003 J=L,N
                                                S=0.0
                                                DO 95004 K=I,M
                                                        S=S+A(K,I)*A(K,J
     +                                                  )
95004                                           CONTINUE
                                                F=S/H
                                                DO 95005 K=I,M
                                                        A(K,J)=A(K,J)+F*
     +                                                  A(K,I)
95005                                           CONTINUE
95003                                   CONTINUE
                                ENDIF
                                DO 95006 K= I,M
                                        A(K,I)=SCALE*A(K,I)
95006                           CONTINUE
                        ENDIF
                ENDIF
                W(I)=SCALE *G
                G=0.0
                S=0.0
                SCALE=0.0
                IF ((I.LE.M).AND.(I.NE.N)) THEN
                        DO 95007 K=L,N
                                SCALE=SCALE+ABS(A(I,K))
95007                   CONTINUE
                        IF (SCALE.NE.0.0) THEN
                                DO 95008 K=L,N
                                        A(I,K)=A(I,K)/SCALE
                                        S=S+A(I,K)*A(I,K)
95008                           CONTINUE
                                F=A(I,L)
                                G=-SIGN(SQRT(S),F)
                                H=F*G-S
                                A(I,L)=F-G
                                DO 95009 K=L,N
                                        RV1(K)=A(I,K)/H
95009                           CONTINUE
                                IF (I.NE.M) THEN
                                        DO 95010 J=L,M
                                                S=0.0
                                                DO 95011 K=L,N
                                                        S=S+A(J,K)*A(I,K
     +                                                  )
95011                                           CONTINUE
                                                DO 95012 K=L,N
                                                        A(J,K)=A(J,K)+S*
     +                                                  RV1(K)
95012                                           CONTINUE
95010                                   CONTINUE
                                ENDIF
                                DO 95013 K=L,N
                                        A(I,K)=SCALE*A(I,K)
95013                           CONTINUE
                        ENDIF
                ENDIF
                ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
95000   CONTINUE
        DO 95014 I=N,1,-1
                IF (I.LT.N) THEN
                        IF (G.NE.0.0) THEN
                                DO 95015 J=L,N
                                        V(J,I)=(A(I,J)/A(I,L))/G
95015                           CONTINUE
                                DO 95016 J=L,N
                                        S=0.0
                                        DO 95017 K=L,N
                                                S=S+A(I,K)*V(K,J)
95017                                   CONTINUE
                                        DO 95018 K=L,N
                                                V(K,J)=V(K,J)+S*V(K,I)
95018                                   CONTINUE
95016                           CONTINUE
                        ENDIF
                        DO 95019 J=L,N
                                V(I,J)=0.0
                                V(J,I)=0.0
95019                   CONTINUE
                ENDIF
                V(I,I)=1.0
                G=RV1(I)
                L=I
95014   CONTINUE
        DO 95020 I=N,1,-1
                L=I+1
                G=W(I)
                IF (I.LT.N) THEN
                        DO 95021 J=L,N
                                A(I,J)=0.0
95021                   CONTINUE
                ENDIF
                IF (G.NE.0.0) THEN
                        G=1.0/G
                        IF (I.NE.N) THEN
                                DO 95022 J=L,N
                                        S=0.0
                                        DO 95023 K=L,M
                                                S=S+A(K,I)*A(K,J)
95023                                   CONTINUE
                                        F=(S/A(I,I))*G
                                        DO 95024 K=I,M
                                                A(K,J)=A(K,J)+F*A(K,I)
95024                                   CONTINUE
95022                           CONTINUE
                        ENDIF
                        DO 95025 J=I,M
                                A(J,I)=A(J,I)*G
95025                   CONTINUE
                ELSE
                        DO 95026 J= I,M
                                A(J,I)=0.0
95026                   CONTINUE
                ENDIF
                A(I,I)=A(I,I)+1.0
95020   CONTINUE
        DO 95027 K=N,1,-1
                DO 95028 ITS=1,30
                        DO 95029 L=K,1,-1
                                NM=L-1
                                IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO
     +                          2
                                IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 
     +                          1
95029                   CONTINUE
1                       C=0.0
                        S=1.0
                        DO 95030 I=L,K
                                F=S*RV1(I)
                                IF ((ABS(F)+ANORM).NE.ANORM) THEN
                                        G=W(I)
                                        H=SQRT(F*F+G*G)
                                        W(I)=H
                                        H=1.0/H
                                        C= (G*H)
                                        S=-(F*H)
                                        DO 95031 J=1,M
                                                Y=A(J,NM)
                                                Z=A(J,I)
                                                A(J,NM)=(Y*C)+(Z*S)
                                                A(J,I)=-(Y*S)+(Z*C)
95031                                   CONTINUE
                                ENDIF
95030                   CONTINUE
2                       Z=W(K)
                        IF (L.EQ.K) THEN
                                IF (Z.LT.0.0) THEN
                                        W(K)=-Z
                                        DO 95032 J=1,N
                                                V(J,K)=-V(J,K)
95032                                   CONTINUE
                                ENDIF
                                GO TO 3
                        ENDIF
                        IF (ITS.EQ.30)
     &			    stop 'No convergence in 30 iterations'
                        X=W(L)
                        NM=K-1
                        Y=W(NM)
                        G=RV1(NM)
                        H=RV1(K)
                        F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
                        G=SQRT(F*F+1.0)
                        F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
                        C=1.0
                        S=1.0
                        DO 95033 J=L,NM
                                I=J+1
                                G=RV1(I)
                                Y=W(I)
                                H=S*G
                                G=C*G
                                Z=SQRT(F*F+H*H)
                                RV1(J)=Z
                                C=F/Z
                                S=H/Z
                                F= (X*C)+(G*S)
                                G=-(X*S)+(G*C)
                                H=Y*S
                                Y=Y*C
                                DO 95034 NM=1,N
                                        X=V(NM,J)
                                        Z=V(NM,I)
                                        V(NM,J)= (X*C)+(Z*S)
                                        V(NM,I)=-(X*S)+(Z*C)
95034                           CONTINUE
                                Z=SQRT(F*F+H*H)
                                W(J)=Z
                                IF (Z.NE.0.0) THEN
                                        Z=1.0/Z
                                        C=F*Z
                                        S=H*Z
                                ENDIF
                                F= (C*G)+(S*Y)
                                X=-(S*G)+(C*Y)
                                DO 95035 NM=1,M
                                        Y=A(NM,J)
                                        Z=A(NM,I)
                                        A(NM,J)= (Y*C)+(Z*S)
                                        A(NM,I)=-(Y*S)+(Z*C)
95035                           CONTINUE
95033                   CONTINUE
                        RV1(L)=0.0
                        RV1(K)=F
                        W(K)=X
95028           CONTINUE
3               CONTINUE
95027   CONTINUE
        RETURN
        END



        SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
        save
c
        PARAMETER (NMAX=1000)
        DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
c
        DO 95000 J=1,N
                S=0.
                IF(W(J).NE.0.)THEN
                        DO 95001 I=1,M
                                S=S+U(I,J)*B(I)
95001                   CONTINUE
                        S=S/W(J)
                ENDIF
                TMP(J)=S
95000   CONTINUE
        DO 95002 J=1,N
                S=0.
                DO 95003 JJ=1,N
                        S=S+V(J,JJ)*TMP(JJ)
95003           CONTINUE
                X(J)=S
95002   CONTINUE
        RETURN
        END
