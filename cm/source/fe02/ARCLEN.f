      SUBROUTINE ARCLEN(IDO,NBJ,NEL,nl,NPL,NPNE,NVJL,DL,XP,ERROR,*)

C#### Subroutine: ARCLEN
C###  Description:
C###    ARCLEN calculates arc length DL(3) of global line segment using
C###    Gaussian quadrature.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),NEL(0:NELM),
     '  nl,NPL(5,0:3),NPNE(NNM,NBFM,NEM),NVJL(4,NJM)
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IG(4),K,N,nb,ne,ng,NGA,ni,ni2,ni3,NITB,nj,nk,nn,np,nt,nv
      REAL*8 DERIV,SUM,W,WG_LOCAL(10),XI,XIGG(10),XN_LOCAL(2,3,4)
      LOGICAL FOUND

      DATA XIGG/0.5000000000000D0,0.2113248654051D0,0.7886751345948D0,
     '          0.1127016653792D0,0.5000000000000D0,0.8872983346207D0,
     '          0.0694318442029D0,0.3300094782075D0,0.6699905217924D0,
     '          0.9305681557970D0/
      DATA WG_LOCAL/1.0000000000000D0,0.5000000000000D0,0.50000000000D0,
     '          0.2777777777778D0,0.4444444444444D0,0.2777777777778D0,
     '          0.1739274225687D0,0.3260725774313D0,0.3260725774313D0,
     '          0.1739274225687D0/
      DATA   IG/0,1,3,6/,NGA/4/

      CALL ENTERS('ARCLEN',*9999)

      IF(JTYP2B.EQ.1.AND.NPL(4,0).LT.0) GOTO 9998
C *** Calculate array XN_LOCAL of coords at line nodes
      ni=NPL(1,0)
      nt=2
      IF(NPL(1,1).EQ.2) nt=3
      IF(NPL(1,1).EQ.3) nt=4
      DO N=1,nt
        np=NPL(N+1,1)
        DO nj=1,NJT
          nv=NVJL(N,nj)
          XN_LOCAL(1,nj,N)=XP(1,nv,nj,np)
C AJP 4-6-93   IF(NPL(1,nj).EQ.4) THEN
          IF((NPL(1,nj).EQ.4).OR.(NPL(1,nj).EQ.6)
     '      .OR.(NPL(1,nj).EQ.7)) THEN
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
              ERROR='Could not find local node in ARCLEN'
              GOTO 9999
            ENDIF
            DO nk=2,NKT(nn,nb)
              IF(NITB.EQ.1) THEN
                IF(IDO(nk,nn,ni,nb).EQ.2) THEN
                  XN_LOCAL(2,nj,N)=XP(NPL(N+3,1),nv,nj,np)
                ENDIF
              ELSE IF(NITB.EQ.2) THEN
                IF(IDO(nk,nn,ni,nb).EQ.2.AND.IDO(nk,nn,ni2,nb).EQ.1)
     '            THEN
                  XN_LOCAL(2,nj,N)=XP(nk,nv,nj,np)
                ENDIF
              ELSE IF(NITB.EQ.3) THEN
                IF(IDO(nk,nn,ni,nb).EQ.2.AND.IDO(nk,nn,ni2,nb).EQ.1
     '            .AND.IDO(nk,nn,ni3,nb).EQ.1) THEN
                  XN_LOCAL(2,nj,N)=XP(nk,nv,nj,np)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C           Critical section is not essential.
CC$            call mp_setlock()
            WRITE(OP_STRING,'('' XN_LOCAL(k,'',I1,'','',I1,'')='','
     '        //'2E12.3)') nj,n,(XN_LOCAL(k,nj,n),k=1,2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
          ENDIF
        ENDDO
      ENDDO

      IF(ITYP10(1).GT.1) THEN
        DO N=2,nt
          IF(ITYP10(1).LT.4) THEN
            IF(ni.EQ.1.AND.XN_LOCAL(1,2,N).LE.XN_LOCAL(1,2,1))
     '        XN_LOCAL(1,2,N)=XN_LOCAL(1,2,N)+2.0d0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,2,N)-XN_LOCAL(1,2,1)).GT.PI)
     '        XN_LOCAL(1,2,1)=XN_LOCAL(1,2,1)+2.0d0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,2,1)-XN_LOCAL(1,2,N)).GT.PI)
     '        XN_LOCAL(1,2,N)=XN_LOCAL(1,2,N)+2.0d0*PI
          ELSE IF(ITYP10(1).EQ.4) THEN
            IF(ni.EQ.1.AND.XN_LOCAL(1,3,N).GE.XN_LOCAL(1,3,1))
     '        XN_LOCAL(1,3,1)=XN_LOCAL(1,3,1)+2.0d0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,3,N)-XN_LOCAL(1,3,1)).GT.PI)
     '        XN_LOCAL(1,3,1)=XN_LOCAL(1,3,1)+2.0d0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,3,1)-XN_LOCAL(1,3,N)).GT.PI)
     '        XN_LOCAL(1,3,N)=XN_LOCAL(1,3,N)+2.0d0*PI
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

C ***   Calculate derivatives of arc-length/angle wrt Xi at Xi=0,1
C     DO nb=1,NBFT
C     IF(NBI(nb).EQ.1) THEN
C     ELSE IF(NBI(nb).EQ.4) THEN
C     CALL ANGSCA(NPL,DL(1,nl),XP,ERROR,*9999)
C     ELSE IF(NBI(nb).EQ.5) THEN
C     XI=0.0D0
C     CALL ARCDER(NPL,DERIV,DL(1,nl),XI,XN_LOCAL,ERROR,*9999)
C     DL(1,nl)=DERIV
C     XI=1.0D0
C     CALL ARCDER(NPL,DERIV,DL(1,nl),XI,XN_LOCAL,ERROR,*9999)
C     DL(2,nl)=DERIV
C     ENDIF
C     ENDDO

C *** Calculate total arclength by Gaussian quadrature
      SUM=0.0d0
      DO ng=1,NGA
        XI=XIGG(IG(NGA)+ng)
        W=WG_LOCAL(IG(NGA)+ng)
        CALL ARCDER(NPL,DERIV,DL(1,nl),XI,XN_LOCAL,ERROR,*9999)
        SUM=SUM+W*DERIV
      ENDDO
      DL(3,nl)=SUM
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C       Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' DL(3,nl)='',E13.5)') DL(3,nl)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(JTYP2B.EQ.1.AND.NPL(4,0).GT.0) THEN
        DL(3,NPL(4,0))=DL(3,nl)
      ENDIF

 9998 CALL EXITS('ARCLEN')
      RETURN
 9999 CALL ERRORS('ARCLEN',ERROR)
      CALL EXITS('ARCLEN')
      RETURN 1
      END


