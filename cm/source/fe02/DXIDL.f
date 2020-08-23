      SUBROUTINE DXIDL(IDO,NBJ,NEL,NPL,DL,DXI,XP,ERROR,*)

C#### Subroutine: DXIDL
C###  Description:
C###    DXIDL calculates arc lengths DL(3,nl) of global line segment
C###    nl, from derivatives wrt Xi given by DXI(n,nj), using Gaussian
C###    quadrature.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),NEL(0:NELM),
     '  NPL(5,0:3)
      REAL*8 DL(3),DXI(2,3),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,IG(4),K,N,nb,ne,ng,NGA,ni,ni2,ni3,NITB,nj,nk,np,nt
      REAL*8 PH3,PL1,PL2,PL3,SL,SM,SUM1,SUM3,W,WG_LOCAL(10),
     '  XA_LOCAL(3,3),XI,XIGG(10),XN_LOCAL(2,3,3)

      DATA XIGG/0.5000000000000D0,0.2113248654051D0,0.7886751345948D0,
     '          0.1127016653792D0,0.5000000000000D0,0.8872983346207D0,
     '          0.0694318442029D0,0.3300094782075D0,0.6699905217924D0,
     '          0.9305681557970D0/
      DATA WG_LOCAL/1.0000000000000D0,0.5000000000000D0,0.50000000000D0,
     '          0.2777777777778D0,0.4444444444444D0,0.2777777777778D0,
     '          0.1739274225687D0,0.3260725774313D0,0.3260725774313D0,
     '          0.1739274225687D0/
      DATA   IG/0,1,3,6/,NGA/4/

      CALL ENTERS('DXIDL',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' DL(3)='',E13.5)')
     '    DL(3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      IF(NPL(2,1).EQ.NPL(3,1)) GO TO 80

      ni=NPL(1,0)
      nt=2
      IF(NPL(1,1).EQ.2) nt=3
      IF(NPL(1,1).EQ.3) nt=4
      DO N=1,nt
        np=NPL(N+1,1)
        DO nj=1,NJT
          XN_LOCAL(1,nj,N)=XP(1,1,nj,np)
          IF(NPL(1,nj).EQ.4) THEN
            ne=NEL(1)
            nb=NBJ(nj,ne)
            NITB=NIT(nb)
            IF(NITB.EQ.2) THEN
              ni2=1+MOD(ni,2)
            ELSE IF(NITB.EQ.3) THEN
              ni2=1+MOD(ni,3)
              ni3=1+MOD(ni2,3)
            ENDIF
            DO nk=2,NKT(0,nb)
              IF(NITB.EQ.1) THEN
                IF(IDO(nk,1,ni,nb).EQ.2) THEN
                  XN_LOCAL(2,nj,N)=DXI(N,nj)
                ENDIF
              ELSE IF(NITB.EQ.2) THEN
                IF(IDO(nk,1,ni,nb).EQ.2.AND.IDO(nk,1,ni2,nb).EQ.1)
     '            THEN
                  XN_LOCAL(2,nj,N)=DXI(N,nj)
                ENDIF
              ELSE IF(NITB.EQ.3) THEN
                IF(IDO(nk,1,ni,nb).EQ.2.AND.IDO(nk,1,ni2,nb).EQ.1
     '            .AND.IDO(nk,1,ni3,nb).EQ.1) THEN
                  XN_LOCAL(2,nj,N)=DXI(N,nj)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' XN_LOCAL(K,'',I1,'','',I1,'')='','
     '        //'2E12.3)')
     '        nj,N,(XN_LOCAL(K,nj,N),K=1,2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDDO
      ENDDO
      IF(ITYP10(1).GT.1) THEN
        DO N=2,nt
          IF(ITYP10(1).LT.4) THEN
            IF(ni.EQ.1.AND.XN_LOCAL(1,2,N).LE.XN_LOCAL(1,2,1))
     '        XN_LOCAL(1,2,N)=XN_LOCAL(1,2,N)+2.0D0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,2,N)-XN_LOCAL(1,2,1)).GT.PI)
     '        XN_LOCAL(1,2,1)=XN_LOCAL(1,2,1)+2.0D0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,2,1)-XN_LOCAL(1,2,N)).GT.PI)
     '        XN_LOCAL(1,2,N)=XN_LOCAL(1,2,N)+2.0D0*PI
          ELSE IF(ITYP10(1).EQ.4) THEN
            IF(ni.EQ.1.AND.XN_LOCAL(1,3,N).GE.XN_LOCAL(1,3,1))
     '        XN_LOCAL(1,3,1)=XN_LOCAL(1,3,1)+2.0D0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,3,N)-XN_LOCAL(1,3,1)).GT.PI)
     '        XN_LOCAL(1,3,1)=XN_LOCAL(1,3,1)+2.0D0*PI
            IF(ni.NE.1.AND.(XN_LOCAL(1,3,1)-XN_LOCAL(1,3,N)).GT.PI)
     '        XN_LOCAL(1,3,N)=XN_LOCAL(1,3,N)+2.0D0*PI
          ENDIF
        ENDDO
        IF(ITYP10(1).eq.2.and.ni.EQ.2) THEN
          IF(DABS(XN_LOCAL(1,1,1)).LT.1.0D-6) THEN
            XN_LOCAL(1,2,1)=XN_LOCAL(1,2,2)
          ELSE IF(DABS(XN_LOCAL(1,1,2)).LT.1.0D-6) THEN
            XN_LOCAL(1,2,2)=XN_LOCAL(1,2,1)
          ENDIF
        ELSE IF(ITYP10(1).eq.4.and.ni.EQ.2) THEN
          IF(DABS(XN_LOCAL(1,2,1)).LT.1.0D-6) THEN
            XN_LOCAL(1,3,1)=XN_LOCAL(1,3,2)
          ELSE IF(DABS(XN_LOCAL(1,2,2)).LT.1.0D-6) THEN
            XN_LOCAL(1,3,2)=XN_LOCAL(1,3,1)
          ENDIF
        ENDIF
      ENDIF

      SUM3=0.0D0
      DO ng=1,NGA
        XI=XIGG(IG(NGA)+ng)
        W=WG_LOCAL(IG(NGA)+ng)
        DO nj=1,NJT
          IF(NPL(1,nj).EQ.1) THEN
            DO K=1,3
              XA_LOCAL(K,nj)=0.0D0
              DO N=1,2
                XA_LOCAL(K,nj)=XA_LOCAL(K,nj)+PL1(N,K,XI)*
     '            XN_LOCAL(1,nj,N)
              ENDDO
            ENDDO
          ELSE IF(NPL(1,nj).EQ.2) THEN
            DO K=1,3
              XA_LOCAL(K,nj)=0.0D0
              DO N=1,3
                XA_LOCAL(K,nj)=XA_LOCAL(K,nj)+PL2(N,K,XI)*
     '            XN_LOCAL(1,nj,N)
              ENDDO
            ENDDO
          ELSE IF(NPL(1,nj).EQ.3) THEN
            DO K=1,3
              XA_LOCAL(K,nj)=0.0D0
              DO N=1,4
                XA_LOCAL(K,nj)=XA_LOCAL(K,nj)+PL3(N,K,XI)*
     '            XN_LOCAL(1,nj,N)
              ENDDO
            ENDDO
          ELSE IF(NPL(1,nj).EQ.4) THEN
            DO K=1,3
              XA_LOCAL(K,nj)=0.0D0
              DO N=1,2
                XA_LOCAL(K,nj)=XA_LOCAL(K,nj)+PH3(N,1,K,XI)*
     '                            XN_LOCAL(1,nj,N)
     '                           +PH3(N,2,K,XI)*XN_LOCAL(2,nj,N)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'(1X,3E11.4,3I2,3E11.4)')
     '              (DL(I),I=1,3),nj,K,N,
     '              (XN_LOCAL(I,nj,N),I=1,2),XA_LOCAL(K,nj)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        IF(ITYP10(1).EQ.1) THEN
          IF(NJT.EQ.2) THEN
            SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2
          ELSE IF(NJT.EQ.3) THEN
            SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
          ENDIF
        ELSE IF(ITYP10(1).EQ.2) THEN
          SUM1=XA_LOCAL(2,1)**2+(XA_LOCAL(2,2)*XA_LOCAL(1,1))**2+
     '     XA_LOCAL(2,3)**2
        ELSE IF(ITYP10(1).EQ.3) THEN
          SUM1=XA_LOCAL(2,1)**2+(XA_LOCAL(2,2)*XA_LOCAL(1,1)*
     '         DCOS(XA_LOCAL(1,3)))**2
     '        +(XA_LOCAL(2,3)*XA_LOCAL(1,1))**2
        ELSE IF(ITYP10(1).EQ.4) THEN
          SL=DSINH(XA_LOCAL(1,1))
          SM=DSIN(XA_LOCAL(1,2))
          SUM1=FOCUS*FOCUS*((SL*SL+SM*SM)*(XA_LOCAL(2,1)**2+
     '                      XA_LOCAL(2,2)**2)
     '                     +(SL*SM*XA_LOCAL(2,3))**2)
        ELSE IF(ITYP10(1).EQ.5) THEN
        ENDIF
        SUM3=SUM3+W*DSQRT(SUM1)
      ENDDO
      DL(3)=SUM3

 80   CALL EXITS('DXIDL')
      RETURN
 9999 CALL ERRORS('DXIDL',ERROR)
      CALL EXITS('DXIDL')
      RETURN 1
      END


