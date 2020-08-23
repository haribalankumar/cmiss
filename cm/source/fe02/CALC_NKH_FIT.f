      SUBROUTINE CALC_NKH_FIT(NBH,NEELEM,NHP,NKH,
     '  NPNE,NPNODE,nr,nx,ERROR,*)

C#### Subroutine: CALC_NKH_FIT
C###  Description:
C###    CALC_NKH_FIT calculates the array NKH for fitting problems.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M),NHP(NPM),
     '  NKH(NHM,NPM,NCM,0:NRM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),nr,nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,ne_tot,nh,nhx,nn,noelem,nonode,np
C      INTEGER elem,elemf,NKJF(NKM,NNM,NJM)

      CALL ENTERS('CALC_NKH_FIT',*9999)

C     Initialise NKH for nh's in this nx in current region
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        DO nhx=1,NHP(np)
          nh=NH_LOC(nhx,nx)
          NKH(nh,np,1,nr)=0
        ENDDO !nh
      ENDDO !nonode

C PM 14Aug2002: No need to have separate parameters for face fitting.
C      IF(KTYP11.EQ.1) then
       ne_tot=NEELEM(0)
C      ELSE !face fitting
C        ne_tot=LN(0)
C      ENDIF
      DO noelem=1,ne_tot
C        IF(KTYP11.EQ.1) then
        ne=NEELEM(noelem)
C        ELSE !face fitting
C          ne=LN(noelem)
C        ENDIF
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh,1,ne)
          DO nn=1,NNT(nb)
C            IF(KTYP11.EQ.1) THEN ! elements
            np=NPNE(nn,nb,ne)
C            ELSE !faces
C              elem=NPF(6,ne)
C              elemf=NPF(8,ne)
C              CALL CALC_FACE_INFORMATION_IND(NBJ(1,elem),
C     '          NBJF(1,ne),elemf,
C     '          NKJE(1,1,1,elem),NKEF,NKJF,NNF,
C     '          NPNE(1,1,elem),NPNF,nr,
C     '          NVJE(1,1,1,elem),NVJF,SE(1,1,elem),
C     '          SF,ERROR,*9999)
C              np=NPNF(nn,nb)
C            ENDIF
            IF(NKT(nn,nb).GT.NKH(nh,np,1,nr))
     '        NKH(nh,np,1,nr)=NKT(nn,nb)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' nh='',I1,'' np='',I3,'
     '          //''' NKH(nh,np,nc,nr)='',I2)') nh,np,NKH(nh,np,1,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDDO !nn
        ENDDO !nh
      ENDDO !noelem

C     Find max NKH for zero'th nr location
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        DO nhx=1,NH_LOC(0,nx)
          nh=NH_LOC(nhx,nx)
          IF(NKH(nh,np,1,nr).GT.NKH(nh,np,1,0))
     '      NKH(nh,np,1,0)=NKH(nh,np,1,nr)
        ENDDO !nh
      ENDDO !nonode (np)

      CALL EXITS('CALC_NKH_FIT')
      RETURN
 9999 CALL ERRORS('CALC_NKH_FIT',ERROR)
      CALL EXITS('CALC_NKH_FIT')
      RETURN 1
      END



