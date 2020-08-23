      SUBROUTINE CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr,ERROR,*)

C#### Subroutine: CALC_NENP_1D
C###  Description:
C###    CALC_NENP_1D calculates the list of elements surrounding a node
C###    for 1D elements.
C***    This routine should be much faster than the general routine
C***    CALC_NENP as it is specific to the 1D case.
C***  Created by Merryn Howatson Tawhai, 24 October 1997

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nb,NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NPNE(NNM,NBFM,NEM),nr
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne,noelem,np1,np2
      CHARACTER STRING*255

      CALL ENTERS('CALC_NENP_1D',*9999)

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        np1=NPNE(1,nb,ne)
        np2=NPNE(2,nb,ne)
        NENP(np1,0,nr)=0
        NENP(np2,0,nr)=0
      ENDDO !noelem (ne)
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        np1=NPNE(1,nb,ne)
        np2=NPNE(2,nb,ne)
        NENP(np1,0,nr)=NENP(np1,0,nr)+1
        NENP(np2,0,nr)=NENP(np2,0,nr)+1
        WRITE(STRING,'(''>>Increase NEPM, node'',I5,'' has'',I5)') np1,
     &    NENP(np1,0,nr)
        CALL ASSERT(NENP(np1,0,nr).LE.NEPM,STRING,ERROR,
     &    *9999)
        NENP(np1,NENP(np1,0,nr),nr)=ne
        NENP(np2,NENP(np2,0,nr),nr)=ne
      ENDDO !noelem (ne)

      CALL EXITS('CALC_NENP_1D')
      RETURN
 9999 CALL ERRORS('CALC_NENP_1D',ERROR)
      CALL EXITS('CALC_NENP_1D')
      RETURN 1
      END


