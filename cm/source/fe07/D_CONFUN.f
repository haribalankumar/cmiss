      SUBROUTINE D_CONFUN(nocont,nl,IDO,NBJ,NEL,NONL,
     '  NONY,NPL,NPNE,nr,NVJP,NYNP,
     '  CJACM,CONJAC,DL,XP,PARAMTYPE,ERROR,*)

C#### Subroutine: D_CONFUN
C###  Description:
C###    D_CONFUN returns analytic derivatives of the contraints wrt
C###    the parameters to be optimised.

C**** NOTE : This only works for ITYP10(nr) = 1 (Rectangular Cartesian)
C**** Note : SUM1 is square of deriv of arclength wrt Xi
C****        SUM2 is 0.5* second deriv of arclength wrt Xi
C****        SUM3 is ??
C****        SUM4 is integral of the constraint jacobian wrt arclength
C****        SUM5 is integral of the constraint jacobian wrt nodal
C****             parameters

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),NEL(0:NELM),
     '  nl,NONL(NLM),NONY(0:NOYM,NYM,NRCM),nocont,NPL(5,0:3),
     '  NPNE(NNM,NBFM,NEM),nr,
     '  NVJP(NJM,NPM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CJACM(NCOM,*),CONJAC(NCOM,*),DL(3),XP(NKM,NVM,NJM,NPM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER IG(4),k,nb,ne,NGA,ng,n,nn,ni,NI2,NITB,nj,njj2,nk,
     '  noopti,noy,np,nt,nv,ny,KTOT
      REAL*8 PH3,PHIN(2,2),SUM1,SUM2,SUM3(2,4,3),SUM4,SUM5(2,4,3),W,
     '  WG_LOCAL(10),XA_LOCAL(4,3),XI,XIGG(10),XN_LOCAL(4,3,2),
     '  PL2S1,PL2S3
      LOGICAL FOUND

      DATA XIGG/0.5000000000000d0,0.2113248654051d0,0.7886751345948d0,
     '          0.1127016653792d0,0.5000000000000d0,0.8872983346207d0,
     '          0.0694318442029d0,0.3300094782075d0,0.6699905217924d0,
     '          0.9305681557970d0/
      DATA WG_LOCAL/1.0000000000000d0,0.5000000000000d0,0.50000000000d0,
     '          0.2777777777778d0,0.4444444444444d0,0.2777777777778d0,
     '          0.1739274225687d0,0.3260725774313d0,0.3260725774313d0,
     '          0.1739274225687d0/
      DATA   IG/0,1,3,6/,NGA/4/

      CALL ENTERS('D_CONFUN',*9999)

C      nr=1 !These are temporary and may need fixing
C      fixed RGB 7/4/98
      IF(PARAMTYPE(1:12).EQ.'DATA_FITTING') THEN
        IF(NPL(2,1).EQ.NPL(3,1)) GO TO 9998
        ni=NPL(1,0)
        nt=2
        IF(NPL(1,1).EQ.2) nt=3
        IF(NPL(1,1).EQ.3) nt=4
        DO n=1,nt
          np=NPL(n+1,1)
          DO nj=1,NJT
            DO nv=1,NVJP(nj,np)
              XN_LOCAL(1,nj,n)=XP(1,nv,nj,np)
              IF((NPL(1,nj).EQ.4).OR.(NPL(1,nj).EQ.6).OR.
     '          (NPL(1,nj).EQ.7)) THEN
                ne=NEL(1)
                nb=NBJ(nj,ne)
                NITB=NIT(nb)
                IF(NITB.EQ.2) THEN
                  NI2=1+MOD(ni,2)
                ELSE IF(NITB.EQ.3) THEN
                  ERROR='..NITB=3 Not implented'
                  GOTO 9999
                ENDIF
C               Find local node nn of global node np on element ne
                nn=1
                FOUND=.FALSE.
                DO WHILE((nn.LE.NNT(nb)).AND.(.NOT.FOUND))
                  IF(np.EQ.NPNE(nn,nb,ne)) THEN
                    FOUND=.TRUE.
                  ELSE
                    nn=nn+1
                  ENDIF
                ENDDO
                IF(.NOT.FOUND)THEN
                  ERROR='>>Could not find local node'
                  GOTO 9999
                ENDIF
                DO nk=2,NKT(nn,nb)
                  IF(NITB.EQ.1) THEN
                    IF(IDO(nk,nn,ni,nb).EQ.2) THEN
                      XN_LOCAL(2,nj,n)=XP(NPL(n+3,1),nv,nj,np)
                    ENDIF
                  ELSE IF(NITB.EQ.2) THEN
                    IF(IDO(nk,nn,ni,nb).EQ.2.AND.IDO(nk,nn,NI2,nb)
     '                .EQ.1) THEN
                      XN_LOCAL(2,nj,n)=XP(nk,nv,nj,np)
                    ENDIF
                  ELSE IF(NITB.EQ.3) THEN
                    ERROR='>>NITB=3 Not implented'
                    GOTO 9999
                  ENDIF
                ENDDO !nk
              ENDIF
              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP           CRITICAL(D_CONFUN_1)
                WRITE(OP_STRING,'('' XN_LOCAL(K,'',I1,'','',I1,'
     '            //''')='',2D12.3)') nj,n,(XN_LOCAL(k,nj,n),k=1,2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP           END CRITICAL(D_CONFUN_1)
              ENDIF
            ENDDO !nv
          ENDDO !nj
        ENDDO !n
        SUM4=0.0d0
        DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj2,nr)
          DO nk=1,2
            DO nn=1,2
              SUM5(nn,nk,nj)=0.0d0
            ENDDO !nn
          ENDDO !nk
        ENDDO !nj
        DO ng=1,NGA
          XI=XIGG(IG(NGA)+ng)
          W=WG_LOCAL(IG(NGA)+ng)
          DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj2,nr)
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(D_CONFUN_2)
              WRITE(OP_STRING,'('' nj='',I1)') nj
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP         END CRITICAL(D_CONFUN_2)
            ENDIF
            IF(NPL(1,nj).EQ.4) THEN
              DO k=1,2
                XA_LOCAL(k,nj)=0.0d0
                DO n=1,2
                  XA_LOCAL(k,nj)=XA_LOCAL(k,nj)+PH3(n,1,k,XI)*
     '                              XN_LOCAL(1,nj,n)
     '                             +PH3(n,2,k,XI)*XN_LOCAL(2,nj,n)*DL(n)
                ENDDO !n
              ENDDO !k
C             second deriv wrt xi and total arc length
              XA_LOCAL(3,nj)=0.0d0
              DO n=1,2
                XA_LOCAL(3,nj)=XA_LOCAL(3,nj)+PH3(n,2,2,XI)*
     '            XN_LOCAL(2,nj,n)
              ENDDO !n
C             calculate basis function values
              DO n=1,2
                PHIN(1,n)=PH3(n,1,2,XI)
                PHIN(2,n)=PH3(n,2,2,XI)*DL(n)
              ENDDO
            ELSE IF(NPL(1,nj).EQ.6) THEN
!             !Special Hermite simplex line Apex at node 1
              DO k=1,2
                XA_LOCAL(k,nj)=PL2S1(1,1,k,XI)*XN_LOCAL(1,nj,1)
     '                  +PL2S1(2,1,k,XI)*XN_LOCAL(1,nj,2)
     '                  +PL2S1(2,2,k,XI)*XN_LOCAL(2,nj,2)*DL(2)
              ENDDO
C             second deriv wrt xi and total arc length
              XA_LOCAL(3,nj)=PL2S1(2,2,2,XI)*XN_LOCAL(2,nj,2)
C             calculate basis function values
              PHIN(1,1)=PL2S1(1,1,2,XI)
              PHIN(2,1)=0.0d0
              PHIN(1,2)=PL2S1(2,1,2,XI)
              PHIN(2,2)=PL2S1(2,2,2,XI)*DL(2)
            ELSE IF(NPL(1,nj).EQ.7) THEN
!             !Special Hermite simplex line Apex at node 3
              DO k=1,2
                XA_LOCAL(k,nj)=PL2S3(1,1,k,XI)*XN_LOCAL(1,nj,1)
     '                  +PL2S3(1,2,k,XI)*XN_LOCAL(2,nj,1)*DL(1)
     '                  +PL2S3(2,1,k,XI)*XN_LOCAL(1,nj,2)
              ENDDO
C             second deriv wrt xi and total arc length
              XA_LOCAL(3,nj)=PL2S3(1,2,2,XI)*XN_LOCAL(2,nj,1)
C             calculate basis function values
              PHIN(1,1)=PL2S3(1,1,2,XI)
              PHIN(2,1)=PL2S3(1,2,2,XI)*DL(1)
              PHIN(1,2)=PL2S3(2,1,2,XI)
              PHIN(2,2)=0.0d0
            ELSE
              ERROR='>>Line type not implemented'
              GOTO 9999
            ENDIF
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP         CRITICAL(D_CONFUN_3)
              WRITE(OP_STRING,'('' XA_LOCAL(k,nj)='',3D11.4)')
     '          (XA_LOCAL(k,nj),K=1,4)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP         END CRITICAL(D_CONFUN_3)
            ENDIF
          ENDDO !nj
          IF(ITYP10(nr).EQ.1) THEN
            SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
            SUM2=XA_LOCAL(2,1)*XA_LOCAL(3,1)+
     '        XA_LOCAL(2,2)*XA_LOCAL(3,2)+XA_LOCAL(2,3)*XA_LOCAL(3,3)
            DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj2,nr)
              DO nk=1,2
                DO nn=1,2
                  SUM3(nn,nk,nj)=PHIN(nk,nn)*XA_LOCAL(2,nj)
                ENDDO !nn
              ENDDO !nk
            ENDDO !nj
          ELSE
            ERROR='>>Current ITYP10(nr) not implemented'
            GOTO 9999
          ENDIF
          IF(SUM1.GT.1.0d-6) THEN
            SUM4=SUM4+W*SUM2/DSQRT(SUM1)
            DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj2,nr)
              DO nk=1,2
                DO nn=1,2
                  SUM5(nn,nk,nj)=SUM5(nn,nk,nj)+W*SUM3(nn,nk,nj)/
     '              DSQRT(SUM1)
                ENDDO !nn
              ENDDO !nk
            ENDDO !nj
          ENDIF
        ENDDO !ng
      ENDIF
      noopti=NONL(nl)
      IF(KTYP29.EQ.1) THEN
        CONJAC(nocont,noopti)=SUM4-1.0d0
      ELSE IF(KTYP29.EQ.2) THEN
        CJACM(nocont,noopti)=SUM4-1.0d0
      ENDIF
      DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
        nj=NJ_LOC(NJL_GEOM,njj2,nr)
        DO nn=1,2
          np=NPL(1+nn,1)
          DO nv=1,NVJP(nj,np)
            IF(NPL(1,nj).EQ.4) THEN
              KTOT=2
            ELSE IF(NPL(1,nj).EQ.6) THEN
              IF(nn.EQ.1) THEN
                KTOT=1
              ELSE
                KTOT=2
              ENDIF
            ELSE IF(NPL(1,nj).EQ.7) THEN
              IF(nn.EQ.2) THEN
                KTOT=1
              ELSE
                KTOT=2
              ENDIF
            ENDIF
            DO nk=1,KTOT
              IF(nk.EQ.1) THEN
                ny=NYNP(nk,nv,nj,np,0,1,nr)
              ELSE
                ny=NYNP(NPL(3+nn,1),nv,nj,np,0,1,nr)
              ENDIF
              DO noy=1,NONY(0,ny,2)
                noopti=NONY(noy,ny,2)
                IF(KTYP29.EQ.1) THEN
                  CONJAC(nocont,noopti)=SUM5(nn,nk,nj)
                ELSE IF(KTYP29.EQ.2) THEN
                  CJACM(nocont,noopti)=SUM5(nn,nk,nj)
                ENDIF
              ENDDO !noy
            ENDDO !nk
          ENDDO !nv
        ENDDO !nn
      ENDDO !nj

 9998 CALL EXITS('D_CONFUN')
      RETURN
 9999 CALL ERRORS('D_CONFUN',ERROR)
      CALL EXITS('D_CONFUN')
      RETURN 1
      END


