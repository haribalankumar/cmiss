      SUBROUTINE ARCSCA(IDO,JDER,JEST,JSCA,NBJ,NEL,nl,NPL,NPNE,NVJL,
     '  DL,TOL,XP,ERROR,*)

C#### Subroutine: ARCSCA
C###  Description:
C###    ARCSCA calculates arc length DL(3) of global line segment nl,
C###    using Gaussian quadrature.

C**** Newton iteration is used to calculate arc length.
C**** TOL sets convergence tolerance.
C**** JDER=0,1 as nodal derivs are not,are updated, respec.
C**** JEST=0,1 as initial estimates of DL(3) are not,are given,respec
C**** JSCA=0,1 as values of global arc-length derivs DL(1) & DL(2)
C**** are not, are given, respectively. If JSCA=1 no iter.n is needed.
C**** Global derivs in XP are assumed to be wrt arc-length if JSCA=0.
C**** Note: SUM1 is square of deriv of arclength wrt Xi
C****       SUM2 is 0.5* second deriv of arclength wrt Xi
C****       SUM3 is arclength estimate
C****       SUM4 is used is Newton method

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),JDER,JEST,JSCA,NBJ(NJM,NEM),
     '  NEL(0:NELM),nl,NPL(5,0:3),NPNE(NNM,NBFM,NEM),NVJL(4,NJM)
      REAL*8 DL(3,NLM),TOL,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IG(4),k,it,IT_count,ITMAX,n,ne,nb,ng,NGA,ni,ni2,ni3,
     '  NITB,nj,nk,nn,np,nt,nv
      REAL*8 CL,CM,DA,PH3,PL1,PL2,PL2S1,PL2S3,PL3,SL,SM,SUM1,SUM2,SUM3,
     '  SUM4,W,WG_LOCAL(10),XA_LOCAL(4,3),XI,XIGG(10),XN_LOCAL(2,3,4)
      LOGICAL FOUND,LINEAR

      DATA XIGG/0.5000000000000D0,0.2113248654051D0,0.7886751345948D0,
     '  0.1127016653792D0,0.5000000000000D0,0.8872983346207D0,
     '  0.0694318442029D0,0.3300094782075D0,0.6699905217924D0,
     '  0.9305681557970D0/
      DATA WG_LOCAL/1.0000000000000D0,0.5000000000000D0,0.50000000000D0,
     '  0.2777777777778D0,0.4444444444444D0,0.2777777777778D0,
     '  0.1739274225687D0,0.3260725774313D0,0.3260725774313D0,
     '  0.1739274225687D0/
      DATA   IG/0,1,3,6/,NGA/4/

      CALL ENTERS('ARCSCA',*9999)

C     Check if line is mapped to another line
      IF(JTYP2B.EQ.1.AND.NPL(4,0).LT.0) GOTO 9998

C     initialise DL for current line nl
      DO k=1,3
        DL(k,nl)=0.0d0
      ENDDO !k
      IF(DOP.AND.JDER.EQ.1) THEN
        WRITE(OP_STRING,'('' Update nodal derivatives'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(JSCA.EQ.0) THEN
        ITMAX=20
      ELSE IF(JSCA.EQ.1) THEN
        ITMAX=1
      ENDIF
C KAT: Line may still have length as it may loop back to the same node
C      IF(NPL(2,1).EQ.NPL(3,1)) GO TO 80
      ni=NPL(1,0)
      nt=2
      IF(NPL(1,1).EQ.2) nt=3
      IF(NPL(1,1).EQ.3) nt=4
      DO n=1,nt
        np=NPL(n+1,1)
        DO nj=1,NJT
          nv=NVJL(N,nj)
          XN_LOCAL(1,nj,n)=XP(1,nv,nj,np)
          IF((NPL(1,nj).EQ.4).OR.(NPL(1,nj).EQ.6).OR.(NPL(1,nj).EQ.7))
     '      THEN
            ne=NEL(1)
            nb=NBJ(nj,ne)
            NITB=NIT(nb)
            IF(NITB.EQ.2) THEN
              ni2=1+MOD(ni,2)
            ELSE IF(NITB.EQ.3) THEN
              ni2=1+MOD(ni,3)
              ni3=1+MOD(ni2,3)
            ENDIF
!           !Final local node nn of global node np on element ne
            nn=1
            FOUND=.FALSE.
            DO WHILE((nn.LE.NNT(nb)).AND.(.NOT.FOUND))
              IF(np.EQ.NPNE(nn,nb,ne))THEN
                FOUND=.TRUE.
              ELSE
                nn=nn+1
              ENDIF
            ENDDO
            IF(.NOT.FOUND)THEN
              ERROR='Could not find local node in ARCSCA'
              GOTO 9999
            ENDIF
            DO nk=2,NKT(nn,nb)
              IF(NITB.EQ.1) THEN
                IF(IDO(nk,nn,ni,nb).EQ.2) THEN
                  XN_LOCAL(2,nj,n)=XP(NPL(n+3,1),nv,nj,np)
                ENDIF
              ELSE IF(NITB.EQ.2) THEN
                IF(IDO(nk,nn,ni,nb).EQ.2.AND.IDO(nk,nn,ni2,nb).EQ.1)
     '            THEN
                  XN_LOCAL(2,nj,n)=XP(nk,nv,nj,np)
                ENDIF
              ELSE IF(NITB.EQ.3) THEN
                IF(IDO(nk,nn,ni,nb).EQ.2.AND.IDO(nk,nn,ni2,nb).EQ.1.
     '            AND.IDO(nk,nn,ni3,nb).EQ.1) THEN
                  XN_LOCAL(2,nj,n)=XP(nk,nv,nj,np)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' XN_LOCAL(K,'',I1,'','',I1,'')='','
     '        //'2E12.3)')nj,n,(XN_LOCAL(k,nj,n),k=1,2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDDO
      ENDDO
      IF(ITYP10(1).GT.1) THEN
        DO n=2,nt
          IF(ITYP10(1).LT.4) THEN
            IF(ni.EQ.1.AND.XN_LOCAL(1,2,n).LE.XN_LOCAL(1,2,1))
     '        XN_LOCAL(1,2,N)=XN_LOCAL(1,2,N)+2.0d0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,2,n)-XN_LOCAL(1,2,1)).GT.PI)
     '        XN_LOCAL(1,2,1)=XN_LOCAL(1,2,1)+2.0d0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,2,1)-XN_LOCAL(1,2,N)).GT.PI)
     '        XN_LOCAL(1,2,n)=XN_LOCAL(1,2,n)+2.0d0*PI
          ELSE IF(ITYP10(1).EQ.4) THEN
            IF(ni.EQ.1.AND.XN_LOCAL(1,3,n).GE.XN_LOCAL(1,3,1))
     '        XN_LOCAL(1,3,1)=XN_LOCAL(1,3,1)+2.0d0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,3,n)-XN_LOCAL(1,3,1)).GT.PI)
     '        XN_LOCAL(1,3,1)=XN_LOCAL(1,3,1)+2.0d0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,3,1)-XN_LOCAL(1,3,n)).GT.PI)
     '        XN_LOCAL(1,3,n)=XN_LOCAL(1,3,n)+2.0d0*PI
          ENDIF
        ENDDO
        IF(ITYP10(1).eq.2.and.ni.EQ.2) THEN
          IF(DABS(XN_LOCAL(1,1,1)).LT.1.0d-6) THEN
            XN_LOCAL(1,2,1)=XN_LOCAL(1,2,2)
          ELSE IF(DABS(XN_LOCAL(1,1,2)).LT.1.0d-6) THEN
            XN_LOCAL(1,2,2)=XN_LOCAL(1,2,1)
          ENDIF
        ELSE IF(ITYP10(1).eq.4.and.ni.EQ.2) THEN
          IF(DABS(XN_LOCAL(1,2,1)).LT.1.0d-6) THEN
            XN_LOCAL(1,3,1)=XN_LOCAL(1,3,2)
          ELSE IF(DABS(XN_LOCAL(1,2,2)).LT.1.0d-6) THEN
            XN_LOCAL(1,3,2)=XN_LOCAL(1,3,1)
          ENDIF
        ENDIF
      ENDIF

      IF(JEST.EQ.0) THEN
        SUM2=0.0d0
        DO ng=1,NGA
          XI=XIGG(IG(NGA)+ng)
          W=WG_LOCAL(IG(NGA)+ng)
          DO nj=1,3
            DO k=1,2
              XA_LOCAL(k,nj)=0.0d0
              IF(nj.LE.NJT) XA_LOCAL(k,nj)=PL1(1,k,XI)*XN_LOCAL(1,nj,1)
     '                              +PL1(2,k,XI)*XN_LOCAL(1,nj,2)
            ENDDO
          ENDDO
          IF(ITYP10(1).EQ.1) THEN
            SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
          ELSE IF(ITYP10(1).EQ.2) THEN
            SUM1=XA_LOCAL(2,1)**2+(XA_LOCAL(1,1)*XA_LOCAL(2,2))**2
     '                      +XA_LOCAL(2,3)**2
          ELSE IF(ITYP10(1).EQ.3) THEN
            SUM1=XA_LOCAL(2,1)**2+(XA_LOCAL(1,1)*XA_LOCAL(2,2)
     '                  *DCOS(XA_LOCAL(1,3)))**2+
     '                  (XA_LOCAL(1,1)*XA_LOCAL(2,3))**2
          ELSE IF(ITYP10(1).EQ.4) THEN
            SL=DSINH(XA_LOCAL(1,1))
            SM=DSIN(XA_LOCAL(1,2))
            SUM1=FOCUS*FOCUS*((SL*SL+SM*SM)*
     '                       (XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2)
     '                       +(SL*SM*XA_LOCAL(2,3))**2)
          ELSE IF(ITYP10(1).EQ.5) THEN
          ENDIF
          SUM2=SUM2+W*DSQRT(SUM1)
        ENDDO
        IF(JSCA.EQ.0) THEN
          DL(1,nl)=SUM2
          DL(2,nl)=SUM2
        ENDIF
        DL(3,nl)=SUM2
      ENDIF

      LINEAR=.TRUE.
      IF(ITYP10(1).NE.1) LINEAR=.FALSE.
      DO nj=1,NJT
        IF(NPL(1,nj).GT.1) LINEAR=.FALSE.
      ENDDO
      IF(LINEAR) THEN
        DL(1,nl)=DL(3,nl)
        DL(2,nl)=DL(3,nl)
        GO TO 80
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Starting line length='',E12.3)') DL(3,nl)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      DO it=1,ITMAX
        IT_count=it
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(/'' Iteration '',I2)') it
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        SUM3=0.0d0
        SUM4=0.0d0
        DO ng=1,NGA
          XI=XIGG(IG(NGA)+ng)
          W=WG_LOCAL(IG(NGA)+ng)
          DO nj=1,NJT
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' nj='',I1)') nj
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            IF(NPL(1,nj).EQ.1) THEN
              DO k=1,2
                XA_LOCAL(k,nj)=0.0d0
                DO n=1,2
                  XA_LOCAL(k,nj)=XA_LOCAL(k,nj)+
     '                           PL1(n,k,XI)*XN_LOCAL(1,nj,n)
                ENDDO
              ENDDO
              XA_LOCAL(3,nj)=0.0d0
              XA_LOCAL(4,nj)=0.0d0
            ELSE IF(NPL(1,nj).EQ.2) THEN
              DO k=1,2
                XA_LOCAL(k,nj)=0.0d0
                DO n=1,3
                  XA_LOCAL(k,nj)=XA_LOCAL(k,nj)+
     '                           PL2(n,k,XI)*XN_LOCAL(1,nj,n)
                ENDDO
              ENDDO
              XA_LOCAL(3,nj)=0.0d0
              XA_LOCAL(4,nj)=0.0d0
            ELSE IF(NPL(1,nj).EQ.3) THEN
              DO k=1,2
                XA_LOCAL(k,nj)=0.0d0
                DO n=1,4
                  XA_LOCAL(k,nj)=XA_LOCAL(k,nj)+
     '                           PL3(n,k,XI)*XN_LOCAL(1,nj,n)
                ENDDO
              ENDDO
              XA_LOCAL(3,nj)=0.0d0
              XA_LOCAL(4,nj)=0.0d0
            ELSE IF(NPL(1,nj).EQ.4) THEN
              DO k=1,2
                XA_LOCAL(k,nj)=0.0d0
                DO n=1,2
                  XA_LOCAL(k,nj)=XA_LOCAL(k,nj)+
     '                              PH3(n,1,k,XI)*XN_LOCAL(1,nj,n)
     '                             +PH3(n,2,k,XI)*
     '                              XN_LOCAL(2,nj,n)*DL(n,nl)
                ENDDO
              ENDDO
              XA_LOCAL(3,nj)=0.0d0
C             second deriv wrt xi and total arc length
              DO n=1,2
                XA_LOCAL(3,nj)=XA_LOCAL(3,nj)+
     '            PH3(n,2,2,XI)*XN_LOCAL(2,nj,n)
              ENDDO
C             also need deriv wrt total arc length in polar coord systems
              XA_LOCAL(4,nj)=0.0d0
              DO n=1,2
                XA_LOCAL(4,nj)=XA_LOCAL(4,nj)+
     '            PH3(n,2,1,XI)*XN_LOCAL(2,nj,n)
              ENDDO
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' XA_LOCAL(k,nj)='',3E11.4)')
     '            (XA_LOCAL(k,nj),k=1,4)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ELSE IF(NPL(1,nj).EQ.6) THEN
!             !Added AJP 2-6-93
!             !Special Hermite simplex line Apex at node 1
              DO k=1,2
                XA_LOCAL(k,nj)=PL2S1(1,1,k,XI)*XN_LOCAL(1,nj,1)
     '                  +PL2S1(2,1,k,XI)*XN_LOCAL(1,nj,2)
     '                  +PL2S1(2,2,k,XI)*XN_LOCAL(2,nj,2)*DL(2,nl)
              ENDDO
C             second deriv wrt xi and total arc length
              XA_LOCAL(3,nj)=PL2S1(2,2,2,XI)*XN_LOCAL(2,nj,2)
C             also need deriv wrt total arclength in polar coord systems
              XA_LOCAL(4,nj)=PL2S1(2,2,1,XI)*XN_LOCAL(2,nj,2)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' XA_LOCAL(k,nj)='',3E11.4)')
     '            (XA_LOCAL(k,nj),k=1,4)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ELSE IF(NPL(1,nj).EQ.7) THEN
!             !Added AJP 2-6-93
!             !Special Hermite simplex line Apex at node 3
              DO k=1,2
                XA_LOCAL(k,nj)=PL2S3(1,1,k,XI)*XN_LOCAL(1,nj,1)
     '                  +PL2S3(1,2,k,XI)*XN_LOCAL(2,nj,1)*DL(1,nl)
     '                  +PL2S3(2,1,k,XI)*XN_LOCAL(1,nj,2)
              ENDDO
C             second deriv wrt xi and total arc length
              XA_LOCAL(3,nj)=PL2S3(1,2,2,XI)*XN_LOCAL(2,nj,1)
C             also need deriv wrt total arclength in polar coord systems
              XA_LOCAL(4,nj)=PL2S3(1,2,1,XI)*XN_LOCAL(2,nj,1)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' XA_LOCAL(k,nj)='',3E11.4)')
     '            (XA_LOCAL(k,nj),k=1,4)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ENDIF
          ENDDO
c cpb 2/8/95 ??? Why is the nj=3 part of XA_LOCAL accessed here with no
c check that njt=3. Have done the rectangular cartesian case as this
c is obvious. Need to check other cases.
          IF(ITYP10(1).EQ.1) THEN
            SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
            SUM2=0.0d0
            DO nj=1,NJT
              SUM2=SUM2+XA_LOCAL(2,nj)*XA_LOCAL(3,nj)
            ENDDO !nj
          ELSE IF(ITYP10(1).EQ.2) THEN
            SUM1=XA_LOCAL(2,1)**2+(XA_LOCAL(2,2)*XA_LOCAL(1,1))**2
            SUM2=XA_LOCAL(2,1)*XA_LOCAL(3,1)+XA_LOCAL(2,2)*
     '           XA_LOCAL(3,2)*XA_LOCAL(1,1)**2
     '          +XA_LOCAL(1,1)*XA_LOCAL(4,1)*XA_LOCAL(2,2)**2
            IF(NJT.EQ.3) THEN
              SUM1=SUM1+XA_LOCAL(2,3)**2
              SUM2=SUM2+XA_LOCAL(2,3)*XA_LOCAL(3,3)
            ENDIF
          ELSE IF(ITYP10(1).EQ.3) THEN
            SUM1=XA_LOCAL(2,1)**2+(XA_LOCAL(2,2)*
     '           XA_LOCAL(1,1)*DCOS(XA_LOCAL(1,3)))**2
     '          +(XA_LOCAL(2,3)*XA_LOCAL(1,1))**2
            SUM2=XA_LOCAL(2,1)*XA_LOCAL(3,1)
     '          +XA_LOCAL(2,2)*XA_LOCAL(3,2)*
     '           (XA_LOCAL(1,1)*DCOS(XA_LOCAL(1,3)))**2
     '          +XA_LOCAL(2,3)*XA_LOCAL(3,3)*XA_LOCAL(1,1)**2
     '          +XA_LOCAL(1,1)*XA_LOCAL(4,1)*(XA_LOCAL(2,2)*
     '          DCOS(XA_LOCAL(1,3)))**2
     '          -XA_LOCAL(4,3)*DCOS(XA_LOCAL(1,3))*
     '          DSIN(XA_LOCAL(1,3))*(XA_LOCAL(1,1)*
     '           XA_LOCAL(2,2))**2
     '          +XA_LOCAL(1,1)*XA_LOCAL(4,1)*XA_LOCAL(2,3)**2
          ELSE IF(ITYP10(1).EQ.4) THEN
            SL=DSINH(XA_LOCAL(1,1))
            SM=DSIN(XA_LOCAL(1,2))
            SUM1=FOCUS*FOCUS*((SL*SL+SM*SM)*(XA_LOCAL(2,1)**2+
     '                        XA_LOCAL(2,2)**2)
     '                       +(SL*SM*XA_LOCAL(2,3))**2)
            CL=DCOSH(XA_LOCAL(1,1))
            CM=COS (XA_LOCAL(1,2))
            SUM2=FOCUS**2*((SL*SL+SM*SM)*(XA_LOCAL(2,1)*XA_LOCAL(3,1)
     '                                   +XA_LOCAL(2,2)*XA_LOCAL(3,2))
     '                        +(SL*SM)**2*XA_LOCAL(2,3)*XA_LOCAL(3,3)
     '         +(SL*CL*XA_LOCAL(4,1)+SM*CM*XA_LOCAL(4,2))*
     '         (XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2)
     '         +(SL*SL*SM*CM*XA_LOCAL(4,2)+SL*CL*SM*SM*
     '           XA_LOCAL(4,1))*XA_LOCAL(2,3)**2)
          ELSE IF(ITYP10(1).EQ.5) THEN
          ENDIF
          SUM3=SUM3+W*DSQRT(SUM1)
          IF(SUM1.GT.1.0d-6) SUM4=SUM4+W*SUM2/DSQRT(SUM1)
        ENDDO !ng

        IF(JSCA.EQ.0) THEN
          !arc-length derivs DL(1) & DL(2) not given
          !& must be calculated with Newton
          !Note: see notes in FE02 folder for maths on this
          DA=-(DL(3,nl)-SUM3)/(1.0d0-SUM4)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' SUM1='',D11.4,'' SUM2='',D11.4,'' SUM3='',D11.4,'
     '        //''' SUM4='',D11.4,'' DA='',D11.4)')
     '        SUM1,SUM2,SUM3,SUM4,DA
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(DABS(DA).GT.1.0d6) THEN
            WRITE(OP_STRING,'('' Length of line nl='',I4,'
     '        //''' has not converged & is set to unity'')') nl
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DL(3,nl)=1.0d0
            GOTO 80
          ENDIF
          DL(3,nl)=DL(3,nl)+DA !is new arclength
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' DL(3,nl)='',E12.3)') DL(3,nl)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(JDER.EQ.1) THEN !update nodal derivatives
            DO n=1,2
              IF(DABS(XN_LOCAL(2,1,n)).LT.1.0d0) THEN
                XN_LOCAL(2,2,n)=DSQRT(1.0d0-XN_LOCAL(2,1,n)**2)
              ELSE
                XN_LOCAL(2,1,n)=1.0d0
                XN_LOCAL(2,2,n)=0.0d0
              ENDIF
! Note: XN_LOCAL is updated so as to keep dX_dXi fixed at its previous value
              DO nj=1,NJT
                XN_LOCAL(2,nj,n)=XN_LOCAL(2,nj,n)*DL(n,nl)/DL(3,nl) !new node deriv
              ENDDO                                     !wrt s
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' XN_LOCAL(2,nj,n='',I1,''): '','
     '             //'3E12.3)') n,(XN_LOCAL(2,nj,n),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ENDDO
          ENDIF
          DL(1,nl)=DL(3,nl)
          DL(2,nl)=DL(3,nl)
          IF(DABS(DA).LE.DL(3,nl)*TOL) GOTO 70
        ELSE IF(JSCA.EQ.1) THEN
          DL(3,nl)=SUM3
          GOTO 80
        ENDIF
      ENDDO !iteration
      IF(IT_count.EQ.ITMAX) THEN
        WRITE(OP_STRING,'('' >>WARNING!!! Iteration in ARCSCA has not'
     '    //' converged '',3F12.8)') DA,DL(3,nl),TOL
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

C *** Update derivs (transfer XN_LOCAL(2,nj,n) to XP(nk,nv,nj,np))
 70   IF(JDER.EQ.1) THEN
        ni=NPL(1,0)
        nt=2
        IF(NPL(1,1).EQ.2) nt=3
        IF(NPL(1,1).EQ.3) nt=4
        DO n=1,nt
          np=NPL(n+1,1)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          DO nj=1,NJT
            nv=NVJL(n,nj)
            IF((NPL(1,nj).EQ.4).OR.(NPL(1,nj).EQ.6).OR.
     '        (NPL(1,nj).EQ.7)) THEN
C AJP 3-6-93 IF(NPL(1,nj).EQ.4) THEN
              ne=NEL(1)
              nb=NBJ(nj,ne)
              NITB=NIT(nb)
              IF(NITB.EQ.2) THEN
                ni2=1+MOD(ni,2)
              ELSE IF(NITB.EQ.3) THEN
                ni2=1+MOD(ni,3)
                ni3=1+MOD(ni2,3)
              ENDIF
              DO nk=2,NKT(n,nb)
                IF(NITB.EQ.1) THEN
                  IF(IDO(nk,nn,ni,nb).EQ.2) THEN
                    XP(NPL(n+3,1),nv,nj,np)=XN_LOCAL(2,nj,n)
                  ENDIF
                ELSE IF(NITB.EQ.2) THEN
                  IF(IDO(nk,nn,ni,nb).EQ.2.AND.IDO(nk,nn,ni2,nb).EQ.1)
     '              THEN
                    XP(nk,nv,nj,np)=XN_LOCAL(2,nj,n)
                  ENDIF
                ELSE IF(NITB.EQ.3) THEN
                  IF(IDO(nk,nn,ni,nb).EQ.2.AND.IDO(nk,nn,ni2,nb).EQ.1.
     '              AND.IDO(nk,nn,ni3,nb).EQ.1) THEN
                    XP(nk,nv,nj,np)=XN_LOCAL(2,nj,n)
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

 80   CONTINUE

      IF(JTYP2B.EQ.1.AND.NPL(4,0).GT.0) THEN
        DL(1,NPL(4,0))=-1.0d0*DL(2,nl)
        DL(2,NPL(4,0))=-1.0d0*DL(1,nl)
        DL(3,NPL(4,0))=DL(3,nl)
      ENDIF

 9998 CALL EXITS('ARCSCA')
      RETURN
 9999 CALL ERRORS('ARCSCA',ERROR)
      CALL EXITS('ARCSCA')
      RETURN 1
      END


