      SUBROUTINE OPVARI(iy,NBH,NCLIST,NEELEM,NENQ,NHE,NHP,NKH,NPNODE,
     '  NP_INTERFACE,NQNP,NQS,NQXI,nr,nrg,NVHP,nx,nxg,NXQ,NYNE,NYNP,
     '  NYNR,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,XQ,YP,YQ,ZA,ZP,FIX,FIXQ,TYPE,
     '  ERROR,*)

C#### Subroutine: OPVARI
C###  Description:
C###    OPVARI outputs YP variables.

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'nqloc00.inc'

!     Parameter List
      INTEGER iy,NBH(NHM,NCM,NEM),NCLIST(0:4),NEELEM(0:NE_R_M,0:NRM),
     '  NENQ(0:8,NQM),NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NPNODE(0:NP_R_M,0:NRM),NP_INTERFACE(0:NPM,0:3),NQNP(NPM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),nr,nrg,NVHP(NHM,NPM,NCM),nx,nxg,
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     '  DXDXIQ2(3,3,NQM),XQ(NJM,NQM),
     '  YP(NYM,NIYM),YQ(NYQM,NIQM,NAM),ZA(NAM,NHM,NCM,NEM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL FIX(NYM,NIYFIXM),FIXQ(NYQM,NIYFIXM)
!     Local Variables
      INTEGER COUNT,i,na,nc,ne,nh,nhx,niqV,nk,no_nc,noelem,
     '  nonode,no_nynr,np,npp,nq,nv,ny,nyf,nyfd,nyp,nypd
      REAL*8 FLUX,FLUXD,FLUXDIFF,FLUXDIFFDS,
     '  MAXFLUXDIFF,MAXFLUXDIFFDS,MAXPOTEDIFF,MAXPOTEDIFFDS,
     '  POTE,POTED,POTEDIFF,POTEDIFFDS,SUMSQ,TOTFLUXDIFF,TOTFLUXDIFFDS,
     '  TOTPOTEDIFF,TOTPOTEDIFFDS
      LOGICAL BLOOD,DERIV,EXCLUDE,FOUND

      CALL ENTERS('OPVARI',*9999)

      IF(TYPE(1:4).NE.'GRID') THEN
        WRITE(OP_STRING,'(/'' Region #'',I1)') nr
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(TYPE(1:2).EQ.'YP') THEN
        DO no_nc=1,NCLIST(0)
          nc=NCLIST(no_nc)
          WRITE(OP_STRING,'(/'' nc= '',I1)') nc
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          SUMSQ=0.0d0
         DO no_nynr=1,NYNR(0,0,nc)
            ny=NYNR(no_nynr,0,nc)
            IF(iy.LE.NIYFIXM) THEN
            WRITE(OP_STRING,'('' YP(ny='',I6,'',iy='',I2,'
     '          //''')= '',D25.18,'' FIX(ny='',I6,'',iy='',I2,'
     '          //''')= '',L1)') ny,iy,YP(ny,iy),ny,IY,FIX(ny,iy)
            ELSE
              WRITE(OP_STRING,'('' YP(ny='',I6,'',iy='',I2,'
     '          //''')= '',D25.18)') ny,iy,YP(ny,iy)
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            SUMSQ=SUMSQ+YP(ny,iy)**2
          ENDDO !no_nynr (ny)
          WRITE(OP_STRING,'(/'' Sum of squares = '',E12.5)') SUMSQ
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !no_nc (nc)

      ELSE IF(TYPE(1:2).EQ.'ZP') THEN
        DO no_nc=1,NCLIST(0)
          nc=NCLIST(no_nc)
          WRITE(OP_STRING,'(/'' nc= '',I1)') nc
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL YPZP(iy,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            FORMAT='('' '')'
            WRITE(OP_STRING,FORMAT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nhx=1,NHP(np)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,nc)
                IF(iy.LE.NIYFIXM) THEN
                  FORMAT='('' ZP(nk,'',I1,'','',I1,'','',I3,'','',I1,'
     '              //'''): '',6(D12.4,'' ('',L1,'')''))'
                  WRITE(OP_STRING,FORMAT) nv,nh,np,nc,
     '              (ZP(nk,nv,nh,np,nc),FIX(NYNP(nk,nv,nh,np,0,nc,nr),
     '              iy),nk=1,NKH(nh,np,nc))
                ELSE
                  FORMAT='('' ZP(nk,'',I1,'','',I1,'','',I3,'','',I1,'
     '              //'''): '',6D12.4)'
                  WRITE(OP_STRING,FORMAT) nv,nh,np,nc,
     '              (ZP(nk,nv,nh,np,nc),nk=1,NKH(nh,np,nc))
                ENDIF
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !nv
            ENDDO !nh
          ENDDO !nonode (np)
        ENDDO !no_nc (nc)

        FORMAT='('' '')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            IF(NAT(NBH(nh,1,ne)).GT.0) THEN
              FORMAT='('' ZA(na,'',I1,'',1,'',I5,''): '',6D12.4)'
              WRITE(OP_STRING,FORMAT) nh,ne,
     '          (ZA(na,nh,1,ne),na=1,NAT(NBH(nh,1,ne)))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !nh
        ENDDO !noelem (ne)
      ELSE IF(TYPE(1:4).EQ.'GRID') THEN
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        IF(niqV.EQ.0) niqV=1
        CALL ASSERT(.NOT.UP_NQNP,'>>nq to np mapping not yet set up',
     '    ERROR,*9999)
        CALL ASSERT(.NOT.UP_NENQ,'>>nq to ne mapping not yet set up',
     '    ERROR,*9999)

        !initialise variables
        MAXFLUXDIFF=0.0d0
        MAXFLUXDIFFDS=0.0d0
        MAXPOTEDIFF=0.0d0
        MAXPOTEDIFFDS=0.0d0

        TOTFLUXDIFF=0.0d0
        TOTFLUXDIFFDS=0.0d0
        TOTPOTEDIFF=0.0d0
        TOTPOTEDIFFDS=0.0d0

        COUNT=0

        !calculate differences
        DO npp=1,NPNODE(0,nrg)
          np=NPNODE(npp,nrg)
          FOUND=.FALSE.
          DO i=1,NP_INTERFACE(np,0)
            IF(nr.EQ.NP_INTERFACE(np,i)) FOUND=.TRUE.
          ENDDO

          !Don't add excluded nodes
          EXCLUDE=.FALSE.
          DO i=1,CPLST(0,1)
            IF(np.EQ.CPLST(i,1)) EXCLUDE=.TRUE.
          ENDDO

          IF(FOUND) THEN
            IF(.NOT.EXCLUDE) COUNT=COUNT+1

            nh=NH_LOC(1,nx)
            IF(NKH(nh,np,1).GT.1) THEN
              DERIV=.TRUE.
            ELSE
              DERIV=.FALSE.
            ENDIF

            nyf=NYNP(1,1,nh,np,0,2,nr)
            IF(DERIV) nyfd=NYNP(2,1,nh,np,0,2,nr)
            nyp=NYNP(1,1,nh,np,0,1,nr)
            IF(DERIV) nypd=NYNP(2,1,nh,np,0,1,nr)

            nq=NQNP(np)
            BLOOD=.FALSE.
            IF(FIXQ(nq,3)) BLOOD=.TRUE.
            
            CALL GGRADPHIQDN(NENQ,niqV,nq,NQS,NQXI,NXQ,AQ,CQ(6,nq),
     &        DNUDXQ,DXDXIQ,DXDXIQ2,FLUX,YQ(1,1,1),ERROR,*9999)

            FLUXDIFF=FLUX+YP(nyf,1)

            POTE=YQ(nq,niqV,1)
            POTEDIFF=POTE-YP(nyp,1)

            IF(DERIV) THEN
              CALL NQDS(2,1,NENQ,1,niqV,nq,NQS,NQXI,NXQ,AQ,CQ(6,nq),
     &          DNUDXQ,FLUXD,DXDXIQ,DXDXIQ2,XQ,YQ,.TRUE.,ERROR,*9999)

              IF(BLOOD) THEN
                FLUXDIFFDS=FLUXD-YP(nyfd,1)
              ELSE
                FLUXDIFFDS=FLUXD+YP(nyfd,1)
              ENDIF

              CALL NQDS(2,1,NENQ,1,niqV,nq,NQS,NQXI,NXQ,AQ,CQ(6,nq),
     &          DNUDXQ,POTED,DXDXIQ,DXDXIQ2,XQ,YQ,.FALSE.,ERROR,*9999)

              IF(BLOOD) THEN
                POTEDIFFDS=POTED+YP(nypd,1)
              ELSE
                POTEDIFFDS=POTED-YP(nypd,1)
              ENDIF
            ENDIF

            IF(DOP) THEN
              IF(DERIV) THEN
                WRITE(OP_STRING,'(/'' Grid: '',I6,'' Node: '',I6)')
     '            nq,np
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' grid phi,dphids,dphidn,'
     '            //'d2phidnds: '',4F12.6)') POTE,POTED,FLUX,FLUXD
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)

                WRITE(OP_STRING,'('' bem  phi,dphids,dphidn,'
     '            //'d2phidnds: '',4F12.6)') YP(nyp,1),YP(nypd,1),
     '            YP(nyf,1),YP(nyfd,1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'(/'' Grid: '',I6,'' Node: '',I6)')
     '            nq,np
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' grid phi,dphids,dphidn,'
     '            //'d2phidnds: '',2F12.6)') POTE,FLUX
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)

                WRITE(OP_STRING,'('' bem  phi,dphids,dphidn,'
     '            //'d2phidnds: '',2F12.6)') YP(nyp,1),YP(nyf,1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
C              IF(BLOOD) THEN
C                WRITE(OP_STRING,'('' bem  phi,dphids,dphidn,'
C     '            //'d2phidnds: '',4F12.6)') YP(nyp,1),-YP(nypd,1),
C     '            -YP(nyf,1),YP(nyfd,1)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              ELSE
C                WRITE(OP_STRING,'('' bem  phi,dphids,dphidn,'
C     '            //'d2phidnds: '',4F12.6)') YP(nyp,1),YP(nypd,1),
C     '            -YP(nyf,1),-YP(nyfd,1)
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              ENDIF
            ENDIF

            IF(.NOT.EXCLUDE) THEN
              TOTFLUXDIFF=TOTFLUXDIFF+DABS(FLUXDIFF)
              TOTPOTEDIFF=TOTPOTEDIFF+DABS(POTEDIFF)
              IF(DERIV) THEN
                TOTFLUXDIFFDS=TOTFLUXDIFFDS+DABS(FLUXDIFFDS)
                TOTPOTEDIFFDS=TOTPOTEDIFFDS+DABS(POTEDIFFDS)
              ENDIF

              IF(DABS(FLUXDIFF).GE.DABS(MAXFLUXDIFF)) MAXFLUXDIFF=
     '          DABS(FLUXDIFF)
              IF(DABS(POTEDIFF).GE.DABS(MAXPOTEDIFF)) MAXPOTEDIFF=
     '          DABS(POTEDIFF)
              IF(DERIV) THEN
                IF(DABS(FLUXDIFFDS).GE.DABS(MAXFLUXDIFFDS))
     '            MAXFLUXDIFFDS=DABS(FLUXDIFFDS)
                IF(DABS(POTEDIFFDS).GE.DABS(MAXPOTEDIFFDS))
     '            MAXPOTEDIFFDS=DABS(POTEDIFFDS)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF(COUNT.GT.0) THEN
          TOTFLUXDIFF=TOTFLUXDIFF/DBLE(COUNT)
          TOTPOTEDIFF=TOTPOTEDIFF/DBLE(COUNT)
          IF(DERIV) THEN
            TOTFLUXDIFFDS=TOTFLUXDIFFDS/DBLE(COUNT)
            TOTPOTEDIFFDS=TOTPOTEDIFFDS/DBLE(COUNT)
          ENDIF
        ENDIF

        !write out the results
        WRITE(OP_STRING,'(/'' Grid region is: '',I2,''  Grid class '
     '    //'is: '',I2)') nrg,nxg
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' BE region is: '',I2,''  BE class '
     '    //'is: '',I2)') nr,nx
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        !potential information
        WRITE(OP_STRING,'(/'' Phi_e: '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average absolute difference: '',F12.6)')
     '    TOTPOTEDIFF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Maximum absolute difference: '',F12.6)')
     '    MAXPOTEDIFF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        !arc length derivative of potential
        IF(DERIV) THEN
          WRITE(OP_STRING,'(/'' dPhi_e/ds: '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Average absolute difference: '',F12.6)')
     '      TOTPOTEDIFFDS
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Maximum absolute difference: '',F12.6)')
     '      MAXPOTEDIFFDS
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        !flux information
        WRITE(OP_STRING,'(/'' dPhi_e/dn: '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average absolute difference: '',F12.6)')
     '    TOTFLUXDIFF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Maximum absolute difference: '',F12.6)')
     '    MAXFLUXDIFF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        !arc length derivative of fluxes
        IF(DERIV) THEN
          WRITE(OP_STRING,'(/'' d2Phi_e/dnds: '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Average absolute difference: '',F12.6)')
     '      TOTFLUXDIFFDS
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Maximum absolute difference: '',F12.6)')
     '      MAXFLUXDIFFDS
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('OPVARI')
      RETURN
 9999 CALL ERRORS('OPVARI',ERROR)
      CALL EXITS('OPVARI')
      RETURN 1
      END
