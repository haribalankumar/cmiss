      SUBROUTINE CALC_NENP_VORO(NBJ,NEELEM,NENP,NPNE,NPNODE,nr,
     '  ERROR,*)

C#### Subroutine: CALC_NENP_VORO
C###  Description:
C###    Calculates NENP, ie given a node, find an element that
C###    belongs to it for region nr.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'mxch.inc'
!     Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ne,nn,noelem,nonode,np,nb

      CALL ENTERS('CALC_NENP_VORO',*9999)

      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        DO ne=0,NEPM
          NENP(np,ne,nr)=0
        ENDDO
      ENDDO

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        nb=NBJ(1,ne)
        DO nn=1,NNT(nb)
          np=NPNE(nn,nb,ne)
          NENP(np,0,nr)=NENP(np,0,nr)+1
          CALL ASSERT(NENP(np,0,nr).LE.NEPM,
     '      '>>Increase NEPM',ERROR,*9999)
          NENP(np,NENP(np,0,nr),nr)=ne
        ENDDO
      ENDDO

      IF(NRT.EQ.1) THEN
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          NENP(np,0,0)=NENP(np,0,nr)
          DO ne=1,NENP(np,0,nr)
            NENP(np,ne,0)=NENP(np,ne,nr)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('CALC_NENP_VORO')
      RETURN
 9999 CALL ERRORS('CALC_NENP_VORO',ERROR)
      CALL EXITS('CALC_NENP_VORO')
      RETURN 1
      END


