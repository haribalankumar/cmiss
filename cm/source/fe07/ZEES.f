      SUBROUTINE ZEES(IBT,IDO,INP,LGE,NAN,NBH,NBJ,NBJF,
     '  ne,NFF,NGAP,NHE,NKEF,NMNO,NNF,NPNE,NPNY,nr,NRE,
     '  NW,nx,NXI,NYNE,NYNP,
     '  CE,CG,CGE,CP,D_RE,D_RI3,D_TG,D_ZG,ES,FEXT,
     '  PG,RE1,RE2,RG,SE,WG,XE,XG,YG,ZE,ZE1,ZEA,ZG,ZG1,FIX,ERROR,*)

C#### Subroutine: ZEES
C###  Description:
C###    <HTML>
C###    ZEES calculates element tangent stiffness matrix ES from
C###    current dependent variable array ZE.
C###    <PRE>
C###    KTYP1D=1 : calculation is done algebraically.
C###    KTYP1D=2 : calculation is done by one-sided finite differences.
C###    KTYP1D=3 : calculation is done by central finite differences.
C###    </PRE>
C###    </HTML>

C**** NW=-1 : element is 'dry' (no contribution to global system).
C**** DELTA is global perturbation. DELTA is element pert'n.
C**** ES(mhs,nhs) has scaling factor correction.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'nonl00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),
     '  NBJF(NJM,NFM),ne,NFF(6),
     '  NGAP(NIM,NBM),NHE,NKEF(0:4,16,6,NBFM),NMNO(1:2,0:NOPM),
     '  NNF(0:17,6,NBFM),NPNE(NNM,NBFM,NEM),
     '  NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),NW,nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM),CG(NMM,NGM),CGE(NMM,NGM),
     '  CP(NMM,NPM),D_RE(NSM,NHM,NOPM),
     '  D_RI3(NHM*NSM),D_TG(3,3,NHM*NSM),D_ZG(NHM,NUM,NHM*NSM),
     '  ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM),PG(NSM,NUM,NGM,NBM),
     '  RE1(NSM,NHM),RE2(NSM,NHM),RG(NGM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XE(NSM,NJM),XG(NJM,NUM),YG(NIYGM,NGM),
     '  ZE(NSM,NHM),ZE1(NSM,NHM),ZEA(NSM,NHM),ZG(NHM,NUM),ZG1(NHM,NUM)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ik,GETNYR,mb,mh,mhs,mhx,mj,ms,MSTB,my,my1,nb,nb1,ng,
     '  nh,nhs,nhs1,nhx,nhx1,nh1,nj,njj,nk,nn,ns,ns1,
     &  NSTNAT,ny,ny1
      INTEGER*4 DRES_PTR
      REAL*8 COSHZ,CSS,D,DELTA,DES,HALFDELTA,R,RWG,SINHZ,SS,THETA,
     '  XE_STORE,Y,ZE_STORE,ZE1_STORE,ZEA_STORE
          
      CALL ENTERS('ZEES',*9999)

      IF(NW.EQ.-1) GO TO 999
      IF(KTYP1D.NE.1) THEN
C***    Finite difference calculation of tangent stiffness matrix

C!!! CS 17/8/2002 dropping index on NPNE
C!!! This appears to have been a major bug .....
        IF(ITYP1(nr,nx).EQ.4) THEN !fe40 problems
          CALL CPCG(NW,NBH(NH_LOC(1,nx)),NPNE(1,1,ne),
     '      nr,nx,CE,CG,CGE,CP,PG,
     '      ERROR,*9999)
        ELSE
          CALL CPCG(1,NBH(NH_LOC(1,nx)),NPNE(1,1,ne),
     '      nr,nx,CE,CG,CGE,CP,PG,
     '      ERROR,*9999)
        ENDIF

        IF(KTYP1D.EQ.2) THEN !one-sided finite differences
C         Residual is probably close to linear so a large perturbation
C           can be used (to minimize numerical error).
          DELTA=1.0d-4
C         Calculate reference residual
          IF(ITYP1(nr,nx).EQ.3) THEN !fe30 problems
            CALL ZERE30(NBH,NBJ,NHE,nr,nx,
     '        CG,PG,RE1,WG,XE,XG,YG,ZE,ZG,ERROR,*9999)
            CALL ASSERT(JTYP10.EQ.1,
     '        '>> Finite diffs not implemented for isochoric interp.',
     '        ERROR,*9999)
          ELSE IF(ITYP1(nr,nx).EQ.4) THEN !fe40 problems
            CALL ZERE40(NBH,NBJ,ne,NHE,NW,nx,CE,CG,PG,RE1,SE,
     '        WG,XE,XG,ZE,ZG,ERROR,*9999)
          ELSE IF(ITYP1(nr,nx).EQ.5) THEN !fe50 problems
            IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
              CALL ZERE55(INP,NBH,ne,NHE,nr,nx,
     '          CG,PG,RE1,WG,ZE,ZE1,ZG,ERROR,*9999)
            ELSE !finite elasticity stress analysis
              CALL ZERE50(IBT,IDO,INP,NAN,NBH,NBJ,NBJF,ne,
     '          NFF,NGAP,NHE,NKEF,NNF,NPNE,nr,NRE,NW,nx,NXI,
     '          CE,CG,CP,FEXT,PG,RE1,RG,SE,WG,
     '          XE,XG,YG,ZE,ZEA,ZG,ERROR,*9999)
            ENDIF
          ENDIF
        ELSE !central differences
C         Residual is probably not so linear so a small perturbation
C           must be used.
          HALFDELTA=1.0d-8
          DELTA=2*HALFDELTA
          CALL ASSERT(JTYP10.EQ.1,
     '      '>> Central diffs not implemented for isochoric interp.',
     '      ERROR,*9999)
        ENDIF

C news VJ 17Nov2004: Added if blocks to calculate element stiffness matrices with derivatives
C with respect to deformed (perturbing ZE) or with respect to undeformed variables (perturbing XE)
        nhs=0
        DO nhx=1,NH_LOC(0,nx)
          IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &         nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
            nh=NH_LOC(nhx,nx)
            nb=NBH(nh)
          ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
            nj=NJ_LOC(NJL_GEOM,nhx,nr)
            nb=NBJ(nj)
          ENDIF
          DO ns=1,NST(nb)+NAT(nb)
            nhs=nhs+1
            ny=IABS(LGE(nhs,2)) !local variable number
            ny1=GETNYR(1,NPNY,nr,0,2,ny,NYNE,NYNP) !global variable #
            IF(.NOT.FIX(ny1,1)) THEN

C ***         Store current solution
              IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &             nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                 ZE_STORE=ZE(ns,nhx)
              ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
C new ryas002: Added code to be able to calculate undeformed state when cavity elements are being used for constant
C volume constraint. ZE1 stores undeformed state of cavity. When solving for undeformed, parameters of ZE1 must be
C peturbed to evaluate element stiffness matrix.              
                IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                  ZE1_STORE=ZE1(ns,nj)
                ELSE !all other types of elements
                  XE_STORE=XE(ns,nj)
                ENDIF
              ENDIF
              IF (KTYP5I(nr).EQ.1) THEN !inertia
                ZEA_STORE=ZEA(ns,nhx)
              ENDIF !inertia

              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP           CRITICAL(ZEES_1)
                WRITE(OP_STRING,'(/'' ******* Element '',I5,'
     '            //'''   nh='',I2,'' ns='',I2,'' *******'')') ne,nh,ns
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP           END CRITICAL(ZEES_1)
              ENDIF

              IF(KTYP1D.NE.2) THEN !central differences
C ***           Perturb solution in one direction
C ***           and valuate perturbed residual
                IF(ITYP1(nr,nx).EQ.3) THEN !fe30 problems
                  IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                 nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                    ZE(ns,nhx)=ZE_STORE-HALFDELTA
                  ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                    IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                      ZE1(ns,nj)=ZE1_STORE-HALFDELTA
                    ELSE !all other types of elements
                      XE(ns,nj)=XE_STORE-HALFDELTA
                    ENDIF
                  ENDIF
                ELSE !other problem types
                  IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                 nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                    ZE(ns,nhx)=ZE_STORE-HALFDELTA*SE(ns,nb,ne)
                  ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                    IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                      ZE1(ns,nj)=ZE1_STORE-HALFDELTA*SE(ns,nb,ne)
                    ELSE
                      IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                        ZE1(ns,nj)=ZE1_STORE-HALFDELTA*SE(ns,nb,ne)
                      ELSE  
                        XE(ns,nj)=XE_STORE-HALFDELTA*SE(ns,nb,ne)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
                IF(ITYP1(nr,nx).EQ.3) THEN !fe30 problems
                  CALL ZERE30(NBH,NBJ,NHE,nr,nx,
     '              CG,PG,RE1,WG,XE,XG,YG,ZE,ZG,ERROR,*9999)
                ELSE IF(ITYP1(nr,nx).EQ.4) THEN !fe40 problems
                  CALL ZERE40(NBH,NBJ,ne,NHE,NW,nx,CE,CG,PG,RE1,SE,
     '              WG,XE,XG,ZE,ZG,ERROR,*9999)
                ELSE IF(ITYP1(nr,nx).EQ.5) THEN !fe50 problems
                  IF (KTYP5I(nr).EQ.1) THEN ! inertia
C!!!                Is this restored
                    ZEA(ns,nhx)=ZEA_STORE-(HALFDELTA*SE(ns,nb,ne)
     '                /T_inc)
                  ENDIF
                  IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                     CALL ZERE55(INP,NBH,ne,NHE,nr,nx,
     '                CG,PG,RE1,WG,ZE,ZE1,ZG,ERROR,*9999)
                  ELSE !finite elasticity stress analysis
                    CALL ZERE50(IBT,IDO,INP,NAN,NBH,NBJ,NBJF,ne,
     '                NFF,NGAP,NHE,NKEF,NNF,NPNE,nr,NRE,NW,nx,NXI,
     '                CE,CG,CP,FEXT,PG,RE1,RG,SE,WG,
     '                XE,XG,YG,ZE,ZEA,ZG,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDIF
              IF(nhx.EQ.1.AND.JTYP10.GE.2) THEN
C ***           Special perturbation for isochoric interpolation
C               IF(ITYP11(nr).LE.3) THEN
C                 R=ZE(ns,1)**(1.0D0/ITYP11(nr))
C                 ZE(ns,1)=(R+DELTA*SE(ns,nb,ne))**ITYP11(nr)
C               ELSE IF(ITYP11(nr).EQ.4) THEN
C                 IF(JTYP10.EQ.2) THEN
C                   SS=ZE(ns,1)/(FOCUS*FOCUS)
C                   SINHZ=DSQRT(SS)
C                   COSHZ=DSQRT(1.0D0+SS)
C                   R=DLOG(COSHZ+SINHZ)
C                   Y=DSINH(R+DELTA*SE(ns,nb,ne))
C                   ZE(ns,1)=(FOCUS*Y)*(FOCUS*Y)
C                 ELSE IF(JTYP10.EQ.3) THEN
C                   CSS=ZE(ns,1)/FOCUS**3
C                   DES=CSS*CSS-4.0D0/27.0D0
C                   IF(DES.GT.0.0D0) THEN
C                     D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
C                     COSHZ=D+1.0D0/(3.0D0*D)
C                   ELSE
C                     THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
C                     COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
C                   ENDIF
C                   SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
C                   R=DLOG(COSHZ+SINHZ)
C                   Y=DCOSH(R+DELTA*SE(ns,nb,ne))
C                   ZE(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
C                 ENDIF
C               ENDIF
                IF(NKT(0,nb).GT.1) THEN
                  DO ik=1,NKT(0,nb)
                    DO nn=1,NNT(nb)
                      IF(ns.EQ.(ik+(nn-1)*NKT(0,nb))) THEN
                        nk=ik
C GMH 10/7/96 ns1 is not used
C                       ns1=1+(nn-1)*NKT(0,nb)
                        GOTO 10
                      ENDIF
                    ENDDO
                  ENDDO
                ELSE
                  nk=1
                ENDIF
 10             CONTINUE
                IF(nk.EQ.1) THEN
                  IF(ITYP11(nr).LE.3) THEN
                    IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                      R=ZE(ns,1)**(1.0D0/ITYP11(nr))
                      ZE(ns,1)=(R+DELTA*SE(ns,nb,ne))**ITYP11(nr)
                    ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                      IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                        R=ZE1(ns,1)**(1.0D0/ITYP11(nr))
                        ZE1(ns,1)=(R+DELTA*SE(ns,nb,ne))**ITYP11(nr)                      
                      ELSE  
                        R=XE(ns,1)**(1.0D0/ITYP11(nr))
                        XE(ns,1)=(R+DELTA*SE(ns,nb,ne))**ITYP11(nr)
                      ENDIF
                    ENDIF
                  ELSE IF(ITYP11(nr).EQ.4) THEN
                    IF(JTYP10.EQ.2) THEN
                      IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                    nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                        SS=ZE(ns,1)/(FOCUS*FOCUS)
                      ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) 
     &                  THEN !const vol
                          SS=ZE1(ns,1)/(FOCUS*FOCUS)
                        ELSE  
                          SS=XE(ns,1)/(FOCUS*FOCUS)
                        ENDIF
                      ENDIF
                        SINHZ=DSQRT(SS)
                        COSHZ=DSQRT(1.0D0+SS)
                        R=DLOG(COSHZ+SINHZ)
                        Y=DSINH(R+DELTA*SE(ns,nb,ne))
                      IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                     nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                        ZE(ns,1)=(FOCUS*Y)*(FOCUS*Y)
                      ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) 
     &                  THEN !const vol
                          ZE1(ns,1)=(FOCUS*Y)*(FOCUS*Y)
                        ELSE
                          XE(ns,1)=(FOCUS*Y)*(FOCUS*Y)
                        ENDIF
                      ENDIF
                    ELSE IF(JTYP10.EQ.3) THEN
                      IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                     nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                        CSS=ZE(ns,1)/FOCUS**3
                      ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) 
     &                  THEN !const vol
                          CSS=ZE1(ns,1)/FOCUS**3
                        ELSE
                          CSS=XE(ns,1)/FOCUS**3
                        ENDIF
                      ENDIF
                      DES=CSS*CSS-4.0D0/27.0D0
                      IF(DES.GT.0.0D0) THEN
                        D=((CSS+DSQRT(DES))/2.0D0)**(1.0D0/3.0D0)
                          COSHZ=D+1.0D0/(3.0D0*D)
                      ELSE
                          THETA=DACOS(CSS*DSQRT(27.0D0)/2.0D0)
                          COSHZ=2.0D0/DSQRT(3.0D0)*DCOS(THETA/3.0D0)
                      ENDIF
                      SINHZ=DSQRT(COSHZ*COSHZ-1.0D0)
                      R=DLOG(COSHZ+SINHZ)
                      Y=DCOSH(R+DELTA*SE(ns,nb,ne))
                      IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                     nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                        ZE(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
                      ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) 
     &                  THEN !const vol
                          ZE1(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
                        ELSE
                          XE(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE IF(ITYP11(nr).EQ.5) THEN
                  ENDIF
                ELSE IF(nk.GT.1) THEN
C****           18-Feb-89: Further modification needed if nodal params
C****           include second deriv, in which case need IDO(nk,nn,0,nb)
                  IF(ITYP11(nr).LE.3) THEN
                    IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                   nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                      ZE(ns,1)=ZE(ns,1)+DELTA*SE(ns,nb,ne)
     '                *ITYP11(nr)*R**(ITYP11(nr)-1)
                    ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                      IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) 
     &                THEN !const vol
                        ZE1(ns,1)=ZE1(ns,1)+DELTA*SE(ns,nb,ne)
     '                  *ITYP11(nr)*R**(ITYP11(nr)-1)
                      ELSE
                        XE(ns,1)=XE(ns,1)+DELTA*SE(ns,nb,ne)
     '                  *ITYP11(nr)*R**(ITYP11(nr)-1)
                      ENDIF
                    ENDIF
                  ELSE IF(ITYP11(nr).EQ.4) THEN
                    IF(JTYP10.EQ.2) THEN
                      IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                     nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                        ZE(ns,1)=ZE(ns,1)+DELTA*SE(ns,nb,ne)
     '                  *2.0D0*FOCUS*FOCUS*SINHZ*COSHZ
                      ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) 
     &                  THEN !const vol
                          ZE1(ns,1)=ZE(ns,1)+DELTA*SE(ns,nb,ne)
     '                    *2.0D0*FOCUS*FOCUS*SINHZ*COSHZ
                        ELSE
                          XE(ns,1)=ZE(ns,1)+DELTA*SE(ns,nb,ne)
     '                    *2.0D0*FOCUS*FOCUS*SINHZ*COSHZ
                        ENDIF
                      ENDIF
                    ELSE IF(JTYP10.EQ.3) THEN
                      IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                     nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                        ZE(ns,1)=ZE(ns,1)+DELTA*SE(ns,nb,ne)
     '                    *FOCUS**3*SINHZ*(3.0D0*COSHZ*COSHZ-1.0D0)
                      ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                        IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) 
     &                  THEN !const vol
                          ZE1(ns,1)=ZE1(ns,1)+DELTA*SE(ns,nb,ne)
     '                      *FOCUS**3*SINHZ*(3.0D0*COSHZ*COSHZ-1.0D0)
                        ELSE
                          XE(ns,1)=XE(ns,1)+DELTA*SE(ns,nb,ne)
     '                      *FOCUS**3*SINHZ*(3.0D0*COSHZ*COSHZ-1.0D0)
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE IF(ITYP11(nr).EQ.5) THEN
                  ENDIF
                ENDIF

              ELSE IF(ITYP1(nr,nx).EQ.3) THEN !fe30 problems
C ***           Normal perturbation for non-isochoric interpolation
                IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &               nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                  ZE(ns,nhx)=ZE(ns,nhx)+DELTA
                ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                  IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                    ZE1(ns,nhx)=ZE1(ns,nhx)+DELTA    
                  ELSE
                    XE(ns,nhx)=XE(ns,nhx)+DELTA
                  ENDIF
                ENDIF  
              ELSE
C ***           Normal perturbation for non-isochoric interpolation
                IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &               nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                  ZE(ns,nhx)=ZE(ns,nhx)+DELTA*SE(ns,nb,ne)
                ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                  IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol
                    ZE1(ns,nhx)=ZE1(ns,nhx)+DELTA*SE(ns,nb,ne)
                  ELSE
                    XE(ns,nhx)=XE(ns,nhx)+DELTA*SE(ns,nb,ne)
                  ENDIF
                ENDIF  
                IF (KTYP5I(nr).EQ.1) THEN ! inertia
                  ZEA(ns,nhx)=ZEA(ns,nhx)+(DELTA*SE(ns,nb,ne)
     '              /T_inc)
                ENDIF
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP             CRITICAL(ZEES_2)
                    WRITE(OP_STRING,'('' ******* ZE(ns,nhx)='',D12.5,'
     '                //''' Perturbation='',D12.5)')
     '                ZE(ns,nhx),DELTA*SE(ns,nb,ne)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP               END CRITICAL(ZEES_2)
                  ENDIF
                ENDIF

C ***           Evaluate perturbed residual
                IF(ITYP1(nr,nx).EQ.3) THEN
                  CALL ZERE30(NBH,NBJ,NHE,nr,nx,
     '              CG,PG,RE2,WG,XE,XG,YG,ZE,ZG,ERROR,*9999)
                ELSE IF(ITYP1(nr,nx).EQ.4) THEN
                  CALL ZERE40(NBH,NBJ,ne,NHE,NW,nx,CE,CG,PG,RE2,SE,
     '              WG,XE,XG,ZE,ZG,ERROR,*9999)
                ELSE IF(ITYP1(nr,nx).EQ.5) THEN
                  IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !cnst vol
                    CALL ZERE55(INP,NBH,ne,NHE,nr,nx,
     '                CG,PG,RE2,WG,ZE,ZE1,ZG,ERROR,*9999)
                  ELSE
                    CALL ZERE50(IBT,IDO,INP,NAN,NBH,NBJ,NBJF,ne,
     '                NFF,NGAP,NHE,NKEF,NNF,NPNE,nr,NRE,NW,nx,NXI,
     '                CE,CG,CP,FEXT,PG,RE2,RG,SE,WG,
     '                XE,XG,YG,ZE,ZEA,ZG,ERROR,*9999)
                  ENDIF
                ENDIF

C ***           Return ZE to original value
                IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &               nhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                  ZE(ns,nhx)=ZE_STORE
                ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                  IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN !const vol                
                    ZE1(ns,nhx)=ZE1_STORE
                  ELSE
                    XE(ns,nhx)=XE_STORE
                  ENDIF
                ENDIF  
                IF (KTYP5I(nr).EQ.1) THEN !inertia
                  ZEA(ns,nhx)=ZEA_STORE
                ENDIF !inertia

                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP           CRITICAL(ZEES_3)
                  WRITE(OP_STRING,'('' ******* Diagonal term: '','
     '              //'''RE1(ns,nh)='',D20.10,'' RE2(ns,nh)='',D20.10)')
     '              RE1(ns,nh),RE2(ns,nh)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP           END CRITICAL(ZEES_3)
                ENDIF

C ***         Assemble element stiffness matrix
                mhs=0
                DO mhx=1,NH_LOC(0,nx)
                  IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                 mhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                    mh=NH_LOC(mhx,nx)
                    mb=NBH(mh)
                  ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
C                    mh=NH_LOC(mhx,nx)
                    mj=NJ_LOC(NJL_GEOM,mhx,nr)
                    mb=NBJ(mj)
                  ENDIF
                  DO ms=1,NST(mb)+NAT(mb)
                    mhs=mhs+1
                    my=IABS(LGE(mhs,1)) !local variable number
                    my1=GETNYR(1,NPNY,nr,0,2,my,NYNE,NYNP) !global var #
                    IF(.NOT.FIX(my1,1)) THEN
C ***                 Note: DELTA here is global delta
                      IF(KTYP5L.EQ.1.OR. !solving for deformed coordinates
     &                     mhx.GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !nhx not a geometric variable
                        ES(mhs,nhs)=(RE2(ms,mh)-RE1(ms,mh))/DELTA
                      ELSE IF(KTYP5L.EQ.2) THEN !Solving for undeformed coordinates
                        ES(mhs,nhs)=(RE2(ms,mj)-RE1(ms,mj))/DELTA
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO

              ENDIF !FIX
            ENDDO !ns
          ENDDO !nh
C newe VJ
      ELSE !IF(KTYP1D.EQ.1) THEN
C       algebraic calculation of tangent stiffness matrix

        IF(ITYP1(nr,nx).EQ.3) THEN
          CALL CPCG(1,NBH(NH_LOC(1,nx)),NPNE,nr,nx,CE,CG,CGE,CP,PG,
     '      ERROR,*9999)
          DRES_PTR=0
          IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3.
     '      AND.(ITYP15(nr,nx).EQ.1.OR.ITYP15(nr,nx).EQ.2)) THEN
            CALL ALLOCATE_MEMORY(NSM,1,DPTYPE,DRES_PTR,
     '        MEM_INIT,ERROR,*9999)
          ENDIF
          CALL ZEES30(NBH,NBJ,NHE,nr,nx,CG,%VAL(DRES_PTR),ES,PG,
     '      WG,XE,XG,YG,ZE,ZG,ERROR,*9993)
          IF(ITYP5(nr,nx).EQ.5.AND.ITYP2(nr,nx).EQ.3.
     '      AND.(ITYP15(nr,nx).EQ.1.OR.ITYP15(nr,nx).EQ.2)) THEN
            CALL FREE_MEMORY(DRES_PTR,ERROR,*9999)
          ENDIF
        ELSE IF(ITYP1(nr,nx).EQ.4) THEN
          ERROR=' Analytic derivs not implemented for lin. elasticity'
          GO TO 9999
        ELSE IF(ITYP1(nr,nx).EQ.5) THEN !fe50 problems
          CALL CPCG(1,NBH(NH_LOC(1,nx)),NPNE,nr,nx,CE,CG,CGE,CP,PG,
     '      ERROR,*9999)
          IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
            ERROR='>>Not implemented'
            GO TO 9999
          ELSE
            CALL D_ZERE50('GEOMETRIC_PARAMETERS',IBT,IDO,INP,
     '        NAN,NBH,NBJ,NBJF,NGAP,ne,NFF,NHE,NKEF,NMNO,
     '        NNF,NPNE,nr,NRE,NW,nx,NXI,CE,CG,CP,D_RE,D_RI3,
     '        D_TG,D_ZG,ES,FEXT,PG,RG,SE,WG,XE,XG,
     '        YG,ZE,ZE1,ZG,ZG1,ERROR,*9999)
          ENDIF
        ELSE
          ERROR=' Analytic derivs not implemented'
          GO TO 9999
        ENDIF
      ENDIF

      IF(ITYP1(nr,nx).EQ.3) THEN !fe30 problems
C       Apply scale factors if necessary
        nhs=0
        DO nhx=1,NHE !Melge uses NHE instead of NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh)
          NSTNAT=NST(nb)+NAT(nb)
          mhs=0
          DO mhx=1,NHE !Melge uses NHE instead of NH_LOC(0,nx)
            mh=NH_LOC(mhx,nx)
            mb=NBH(mh)
            MSTB=NST(mb)
            IF(NBI(mb).NE.1.OR.NBI(nb).NE.1) THEN !scale factors not unit
              DO ns=1,NSTNAT
                DO ms=1,MSTB
                  ES(mhs+ms,nhs+ns)=
     '              ES(mhs+ms,nhs+ns)*SE(ms,mb,ne)*SE(ns,nb,ne)
                ENDDO !ms
              ENDDO !ns
            ENDIF !non unit scale factors
            mhs=mhs+MSTB
          ENDDO !mhx
          nhs=nhs+NSTNAT
        ENDDO !nhx
      ENDIF

C     JWF 27/09/04      
      IF((KTYP5G(nr).GE.1).AND.(KTYP5K(nr).GE.1)) THEN !psuedo-viscosity for contact problems
        DO ng=1,NGT(NBH(NH_LOC(1,nx)))
          RWG=RG(ng)*WG(ng,NBH(NH_LOC(1,nx)))
          nhs=0
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NJ_LOC(NJL_GEOM,njj,nr)
            nb=NBH(nj)
            DO ns=1,NST(nb)+NAT(nb)
              nhs=nhs+1
              nhs1=0
              DO nhx1=1,NH_LOC(0,nx)
                nh1=NH_LOC(nhx1,nx)
                NB1=NBH(nh1)
                DO ns1=1,NST(NB1)+NAT(NB1)
                  nhs1=nhs1+1
                  ES(nhs,nhs1)=ES(nhs,nhs1)+VIS_WT*RWG
     &              *PG(ns,1,ng,nb)*PG(ns1,1,ng,nb)
     &              *SE(ns,nb,ne)*SE(ns1,nb,ne)   
                ENDDO !ns1
              ENDDO !nhx1
            ENDDO !ns
          ENDDO !njj                                 
        ENDDO !ng
      ENDIF !ktyp5g   
      

 999  CALL EXITS('ZEES')
      RETURN
 9993 IF(DRES_PTR.NE.0) CALL FREE_MEMORY(DRES_PTR,ERROR,*9999)
 9999 CALL ERRORS('ZEES',ERROR)
      CALL EXITS('ZEES')
      RETURN 1
      END


