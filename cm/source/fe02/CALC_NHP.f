      SUBROUTINE CALC_NHP(NBH,NEELEM,NHE,NHP,NPNE,NPNODE,nr,nx,
     '  ERROR,*)

C#### Subroutine: CALC_NHP
C###  Description:
C###    CALC_NHP calculates the array NHP.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM),
     '  NHP(NPM,0:NRM),nhx,
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,NBP,nc,ne,nn,noelem,nonode,np
      LOGICAL PRESS_NODAL

      CALL ENTERS('CALC_NHP',*9999)

C     Initialise NHP for current region
      DO np=1,NPM
        NHP(np,nr)=0
      ENDDO !np

C***  Set NHP for the problem.
      IF(ITYP1(nr,nx).EQ.5) THEN !FE50 problems
C       check if the hydrostatic press interpolation is nodally based
        PRESS_NODAL=.FALSE.
        IF(NH_LOC(0,nx).GT.NJ_LOC(NJL_GEOM,0,nr)) THEN !non-geom vars included
          noelem=0
          DO WHILE(noelem.LT.NEELEM(0,nr).AND..NOT.PRESS_NODAL)
            noelem=noelem+1
            ne=NEELEM(noelem,nr)
            IF(KTYP52(nr).NE.5) THEN !not inextensible
              NBP=NBH(NH_LOC(NH_LOC(0,nx),nx),1,ne)   !basis fn for hyd pres
            ELSE IF(KTYP52(nr).EQ.5) THEN !incomp+inext
              NBP=NBH(NH_LOC(NH_LOC(0,nx)-1,nx),1,ne) !basis fn for hyd pres
            ENDIF
            IF(NST(NBP).GT.0) PRESS_NODAL=.TRUE.
          ENDDO !noelem
        ENDIF !NH_LOC
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          NHP(np,nr)=NHT50(KTYP52(nr),KTYP51(nr))
          IF(KTYP51(nr).EQ.1.AND.NJT.EQ.1) THEN
            NHP(np,nr)=1
C TVK 31/05/2000 altered for ktyp52=6
          ELSE IF(KTYP51(nr).EQ.3.AND.KTYP52(nr).GE.2.AND.
     '        .NOT.PRESS_NODAL.AND.KTYP52(nr).NE.6) THEN
            IF(KTYP52(nr).NE.5) THEN !not inextensible
              NHP(np,nr)=NHP(np,nr)-1
            ELSE IF(KTYP52(nr).EQ.5) THEN !incomp+inext
              NHP(np,nr)=NHP(np,nr)-2
            ENDIF
          ENDIF
        ENDDO !nonode (np)

      ELSE !all other problem types
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nc=1,NCT(nr,nx)
C JPP 1Apr2003 Assumption can cause problems if first nhx
C  number associated with field using auxiliary basis (i.e. no nodes)
C              nb=NBH(NH_LOC(1,nx),nc,ne)
            DO nhx=1,NH_LOC(0,nx)
              nb=NBH(NH_LOC(nhx,nx),nc,ne)
              DO nn=1,NNT(nb)
                np=NPNE(nn,nb,ne)
                IF(NHE(ne).GT.NHP(np,nr)) NHP(np,nr)=NHE(ne)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                  WRITE(OP_STRING,'('' np='',I4,'' ne='',I4,'
     '              //''' nn='',I2,'' NHP(np,nr)='',I1)')
     '              np,ne,nn,NHP(np,nr)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                ENDIF
              ENDDO !nn
            ENDDO !nhx
          ENDDO !nc
        ENDDO !noelem (ne)
      ENDIF

C     Find max NHP for zero'th nr location
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        IF(NHP(np,nr).GT.NHP(np,0)) NHP(np,0)=NHP(np,nr)
      ENDDO !nonode (np)

      CALL EXITS('CALC_NHP')
      RETURN
 9999 CALL ERRORS('CALC_NHP',ERROR)
      CALL EXITS('CALC_NHP')
      RETURN 1
      END


