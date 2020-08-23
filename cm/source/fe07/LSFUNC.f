      SUBROUTINE LSFUNC(CONY,IBT,IDO,
     & INP,ITER1,LGE,NAN,
     '  NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,NNF,NONY,
     '  NPF,NPNE,NPNODE,NPNY,nr,NRE,NRLIST,NSB,NVHE,NVHP,
     '  NVJE,NW,nx,NXI,
     '  NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
     '  CURVCORRECT,FEXT,FIX,GRR,
     '  PG,RE,RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,
     '  ZE,ZE1,ZP,ZP1,
     '  ALPHA,F,ERROR,*)

C#### Subroutine: LSFUNC
C###  Description:
C###    LSFUNC scales the current search vector with the input
C###    parameter alpha, adds it to the current solution vector, and
C###    calculates the global residual vector. Previously this routine
C###    returned the norm of the residual vector, however it now
C###    returns the ratio of unconstrained to constrained componants to
C###    be consistant with nonlin.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'opti00.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ITER1,LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM,0:NRM),NKB(2,2,2,NNM,NBFM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     & NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),NRLIST(0:NRM),
     '  NSB(NKM,NNM,NBFM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     & NVJE(NNM,NBFM,NJM,NEM),
     & NW(NEM,3),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     & NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM),Z_CONT_LIST(NDM,2,7)
      REAL*8 ALPHA,CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),
     '  CP(NMM,NPM),CURVCORRECT(2,2,NNM,NEM),
     '  F,FEXT(NIFEXTM,NGM,NEM),GRR(NOM),
     '  PG(NSM,NUM,NGM,NBM),RE(NSM,NHM),RG(NGM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YGF(NIYGFM,NGFM,NFM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZA1(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZE(NSM,NHM),ZE1(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM),ZP1(NKM,NVM,NHM,NPM,NCM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Work array pointers
      INTEGER*4 INTWORK_PTR,REALWORK_PTR
      SAVE INTWORK_PTR,REALWORK_PTR

!     Local Variables
      INTEGER IBEG,IEND,ne,nonrlist,no_nynr,ny,
     &  nr_cont
      REAL*8 ERRMAX,RATIO(0:6),YP1_TEMP(NYM),YP5_TEMP(NYM)
      CHARACTER STRING*(MXCH)
      LOGICAL CONTACT,END

      DATA ERRMAX /1.0d-8/

      CALL ENTERS('LSFUNC',*9999)

C *** XSL 18Aug2010 Delete the old code
C      DO no_nynr=1,NYNR(0,0,1) !loop over global variables
C        ny=NYNR(no_nynr,0,1) !is global varible number
C        IF(NPNY(0,ny,0).EQ.1) THEN
C          np=NPNY(4,ny,0)
CC GMH 8/1/97 Update cmgui link
C          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C        ENDIF
C        YP(ny,8) = YP(ny,1) + YP(ny,5)*ALPHA
C      ENDDO !no_nynr (ny)
C
C      CALL YPZP(8,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,
C     '  NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
C      CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
C     '  NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
C     '  NNF,NPF,NPNE,NPNODE,NPNY,nr,NRE,NSB,NVHE,NVHP,NVJE,NW,nx,NXI,
C     '  NYNE,NYNP,NYNR,Z_CONT_LIST,CE,CG,CGE,CP,
C     '  CURVCORRECT,FEXT,PG,RE,
C     '  RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,%VAL(0),Z_CONT,
C     '  ZE,ZE1,ZP,ZP1,
C     '  %VAL(0),FIX,ERROR,*9999)
C
C      RATIO(0)=0.0d0
C      DO no_resid=1,6
C        RSUM_CONSTRAINED(no_resid)=0.0d0
C        RSUM_UNCONSTRAIN(no_resid)=0.0d0
C        RATIO(no_resid)=0.0d0
C      ENDDO
C      DO no=1,NOT(1,1,nr,nx)
C        RESID_TEMP(no)=0.0d0
C      ENDDO                     !no
C      DO no_nynr=1,NYNR(0,1,1) !loop over rows
C        ny1=NYNR(no_nynr,1,1) !is row number
C        IF(NONY(0,ny1,1).GT.0) THEN !free dependent variable
C          ny2=GETNYR(2,NPNY,nr,0,1,ny1,NYNE,NYNP) !is RHS var #
C          DO noy=1,NONY(0,ny1,1) !loop on rows assoc with ny1
C            no=NONY(noy,ny1,1) !is row number for ny1
C            co=CONY(noy,ny1,1) !is coupling coeff for ny1
C            RESID_TEMP(no)=RESID_TEMP(no)+(YP(ny1,4)-YP(ny2,1))*co
C          ENDDO !noy
C        ELSE !bdry cond applied to dependent variable
C ! Stick in code from NONLIN to differentiate dimensional
C ! differences
C          IF(NPNY(0,ny1,0).EQ.1) THEN !node based
C            nk=NPNY(1,ny1,0)
C            nh=NPNY(3,ny1,0)
C            IF(nh.LE.3) THEN !dealing with nodal attributes
C              IF(nk.EQ.1) THEN !then dealing with nodal coordinates
C                RSUM_CONSTRAINED(1)=RSUM_CONSTRAINED(1)+
C     &            DABS(YP(ny1,4))
C              ELSE IF((nk.EQ.2).OR.(nk.EQ.3).OR.(nk.EQ.5)) THEN !dealing with 1st derivs
C                RSUM_CONSTRAINED(2)=RSUM_CONSTRAINED(2)+
C     &            DABS(YP(ny1,4))
C              ELSE IF((nk.EQ.4).OR.(nk.EQ.6).OR.(nk.EQ.7)) THEN
C                RSUM_CONSTRAINED(3)=RSUM_CONSTRAINED(3)+
C     &            DABS(YP(ny1,4))
C              ELSE IF(nk.EQ.8) THEN
C                RSUM_CONSTRAINED(4)=RSUM_CONSTRAINED(4)+
C     &            DABS(YP(ny1,4))
C              ENDIF
C            ELSE !dealing with hydrostatic pressures...
C              RSUM_CONSTRAINED(5)=RSUM_CONSTRAINED(5)+
C    &          DABS(YP(ny1,4))
C            ENDIF 
C          ELSE !dealing with element based params
C            RSUM_CONSTRAINED(6)=RSUM_CONSTRAINED(6)+
C     &        DABS(YP(ny1,4))
C          ENDIF
C !RSUM_CONSTRAIN=RSUM_CONSTRAIN+YP(ny1,4)*YP(ny1,4)
C        ENDIF
C      ENDDO!no_nynr
C      DO no=1,NOT(1,1,nr,nx) !loop over global soln rows
C        DO nyo=1,NYNO(0,no,1)
C          ny1=NYNO(nyo,no,1) !is row number
C          ny2=GETNYR(2,NPNY,nr,0,1,ny1,NYNE,NYNP) !is RHS var #
C          IF(.NOT.FIX(ny2,1)) THEN
C            ! Stick in code from NONLIN to differentiate dimensional
C            ! differences
C            IF(NPNY(0,ny1,0).EQ.1) THEN !node based
C              nk=NPNY(1,ny1,0)
C              nh=NPNY(3,ny1,0)
C              IF(nh.LE.3) THEN !dealing with nodal attributes
C                IF(nk.EQ.1) THEN !then dealing with nodal coordinates
C                  RSUM_UNCONSTRAIN(1)=RSUM_UNCONSTRAIN(1)+
C     &              DABS(RESID_TEMP(no))
C                ELSE IF((nk.EQ.2).OR.(nk.EQ.3).OR.(nk.EQ.5)) THEN !dealing with 1st derivs
C                  RSUM_UNCONSTRAIN(2)=RSUM_UNCONSTRAIN(2)+
C     &              DABS(RESID_TEMP(no))
C                ELSE IF((nk.EQ.4).OR.(nk.EQ.6).OR.(nk.EQ.7)) THEN
C                  RSUM_UNCONSTRAIN(3)=RSUM_UNCONSTRAIN(3)+
C     &              DABS(RESID_TEMP(no))
C                ELSE IF(nk.EQ.8) THEN
C                  RSUM_UNCONSTRAIN(4)=RSUM_UNCONSTRAIN(4)+
C     &              DABS(RESID_TEMP(no))
C                ENDIF
C              ELSE !dealing with hydrostatic pressures...
C                RSUM_UNCONSTRAIN(5)=RSUM_UNCONSTRAIN(5)+
C     &            DABS(RESID_TEMP(no))
C              ENDIF 
C            ELSE !dealing with element based params
C              RSUM_UNCONSTRAIN(6)=RSUM_UNCONSTRAIN(6)+
C     &          DABS(RESID_TEMP(no))
C            ENDIF
C            !RSUM_UNCONSTRAIN=RSUM_UNCONSTRAIN+DABS(RESID_TEMP(no))
C          ENDIF
C        ENDDO
C      ENDDO
C      !backward compatibility
C      RSUM_UNCONSTRAINED=0.0d0
C      DO no_resid=1,6
CC       only include the ratio of a particular residual type if some
CC       d.o.f are constrained
C        IF(RSUM_CONSTRAINED(no_resid).GT.ZERO_TOL) THEN
C          RATIO(no_resid)=
C     &      RSUM_UNCONSTRAIN(no_resid)/RSUM_CONSTRAINED(no_resid)
C          RATIO(0)=RATIO(0)+RATIO(no_resid)
C        ENDIF
C        !backward compatibility
C        RSUM_UNCONSTRAINED=RSUM_UNCONSTRAINED+RSUM_UNCONSTRAIN(no_resid)
C      ENDDO
CC      IF(RSUM_CONSTRAIN.LT.ZERO_TOL) THEN
CC        IF(RSUM_UNCONSTRAIN.LT.ZERO_TOL) THEN
CC           F=ZERO_TOL
CC        ELSE
CC           F=1.0d0/ZERO_TOL
CC        ENDIF
CC      ELSE
CC        F = RSUM_UNCONSTRAIN/RSUM_CONSTRAIN
CC      ENDIF
CCC      F = RSUM_UNCONSTRAIN
C
C      F = RATIO(0)
C
C      !backward compatibility - remove if you want ratio
C      F = RSUM_UNCONSTRAINED


C 14/07/08 XSL YP1 is used instead of YP8 because the projection 
C code for contact is performed on YP1
      IF (.NOT.ALPHA.EQ.0.0d0) THEN
        DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables
          ny=NYNR(no_nynr,0,1,nr) !is global variable number
C 27/11/08 XSL NEWS Update YP5 with alpha to keep consistant with NONLIN when calculating energy norm
C 28/11/08 XSL Save YP1 and YP5 into TEMP variables to avoid building up of numerical error during line search
          YP1_TEMP(ny)=YP(ny,1)
          YP5_TEMP(ny)=YP(ny,5)
          YP(ny,5)=ALPHA*YP(ny,5)
          YP(ny,1)=YP(ny,1)+YP(ny,5)
C 27/11/08 NEWE
        ENDDO !no_nynr (ny)
      ENDIF

C *** 26/06/08 XSL Redo projections for contact problems   

C 29/12/08 XSL NEWS Copied from NONLIN
C Check if it's a contact problem
C nr_solve is passed in as nr
      CONTACT=.FALSE.
      IF(nr.EQ.0) THEN ! coupled problem
        CONTACT=.TRUE.
        DO nonrlist=1,NRLIST(0)
          nr_cont=NRLIST(nonrlist)  
          IF(KTYP5G(nr_cont).GT.0.AND.
     &      CONTACT) THEN!contact
            CONTACT=.TRUE. ! coupled contact
          ELSE
            CONTACT=.FALSE. ! it is coupled problem but not contact
          ENDIF
        ENDDO
      ELSE IF(KTYP5G(nr).GT.0) THEN ! not coupled contact
        CONTACT=.TRUE.
      ENDIF
C 29/12/08 NEWE

C XSL 18Aug2010 NEWS Redo projection for contact problems
      IF(CONTACT) THEN !contact
C Get projection command from .com file
        CALL STRING_TRIM(COM_FILE,IBEG,IEND)
        CO(1)='FEM'
        CO(2)='READ'
        CO(3)='COMMAND'
        NTCO=3
        noco=1
        COQU(3,1)=COM_FILE(IBEG:IEND)
        NTCOQU(1)=0
        NTCOQU(2)=0
        NTCOQU(3)=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,%VAL(INTWORK_PTR),
     '    %VAL(REALWORK_PTR),ERROR,*9999) 
      ENDIF !contact
C XSL NEWE

C 04/07/08 XSL YPZP and ZPRP should be called for each region
      DO nonrlist=1,NRLIST(0)
        nr_cont=NRLIST(nonrlist)
        CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr_cont),
     &   NKH(1,1,1,nr_cont),NPNODE,nr_cont,
     '    NVHP(1,1,1,nr_cont),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        CALL ZPRP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '    NEELEM,NFF,NFFACE,NGAP,NHE,NHP(1,nr_cont),NKB,NKEF,
     &   NKH(1,1,1,nr_cont),NKHE,NKJE,NNB,
     '    NNF,NPF,NPNE,NPNODE,NPNY,nr_cont,NRE,NSB,NVHE,
     '    NVHP(1,1,1,nr_cont),NVJE,NW,nx,NXI,
     '    NYNE,NYNP,NYNR(0,0,1,nr_cont),Z_CONT_LIST,CE,CG,CGE,CP,
     '    CURVCORRECT,FEXT,PG,RE,
     '    RG,SE,WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,%VAL(0),Z_CONT,
     '    ZE,ZE1,ZP,ZP1,
     '    %VAL(0),FIX,ERROR,*9999)
      ENDDO !nonrlist

C 15/07/08 XSL Calculate residual
      IF(CONTACT) THEN !contact
C ITER1 is zero when CALC_CONV_RATIO is called in NONLIN for contact
C but has been incremented when LINMIN is called
C hence a value of zero is passed into CALC_VONC_RATIO,
C instead of ITER1 iteself
        CALL CALC_CONV_RATIO(0,IOOP,NBH,ne,NEELEM,
     &    NONY,NPNY,nr_cont,nr,nx,NYNO,NYNR,
     &    CONY,ERRMAX,GRR,RATIO,YG,YP,CONTACT,
     &    .FALSE.,ERROR,*9999)
      ELSE ! not contact
        IF (ALPHA.EQ.0.0d0) THEN
C To calculate previous residual at alpha=0, use GRR_TEMP (before GRR is updated in SOLVE5)
          CALL CALC_CONV_RATIO(ITER1,IOOP,NBH,ne,NEELEM,
     &      NONY,NPNY,nr_cont,nr,nx,NYNO,NYNR,
     &      CONY,ERRMAX,GRR,RATIO,YG,YP,CONTACT,
     &      .FALSE.,ERROR,*9999)
        ELSE ! use updated GRR
          CALL CALC_CONV_RATIO(ITER1,IOOP,NBH,ne,NEELEM,
     &      NONY,NPNY,nr_cont,nr,nx,NYNO,NYNR,
     &      CONY,ERRMAX,%VAL(GRR_PTR),
     &      RATIO,YG,YP,CONTACT,
     &      .FALSE.,ERROR,*9999)
        ENDIF !ALPHA
      ENDIF !CONTACT

      F = RATIO(0)
      !backward compatibility - remove if you want ratio
C      F = RSUM_UNCONSTRAINED

C 14/07/08 XSL Return YP1 and YP5 to original value
      IF (.NOT.ALPHA.EQ.0.0d0) THEN
        DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables
          ny=NYNR(no_nynr,0,1,nr) !is global variable number
C 27/11/08 XSL NEWS Restore YP1 and YP5
          YP(ny,1)=YP1_TEMP(ny)
          YP(ny,5)=YP5_TEMP(ny)
C 27/11/08 NEWE
        ENDDO !no_nynr (ny)
      ENDIF

      CALL EXITS('LSFUNC')
      RETURN
 9999 CALL ERRORS('LSFUNC',ERROR)
      CALL EXITS('LSFUNC')
      RETURN 1
      END


