      SUBROUTINE NENXI_1D(nb,NEELEM,NENP,NPNE,nr,NXI,ERROR,*)

C#### Subroutine: NENXI_1D
C###  Description:
C###    NENXI_1D finds elements surrounding element ne for 1D elements.
C***    This routine is specific for 1D elements, with no more than two
C***    adjoining elements in the Xi(1) direction, and where adjoining
C***    elements in the Xi(1) direction must have a higher global
C***    element number than the subject element.
C***  Created by Merryn Howatson Tawhai, 24 October 1997

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nb,NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NPNE(NNM,NBFM,NEM),nr,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,ne,ne0,NE_SURROUND,noelem,noelem0,np0,np1,np2
      CHARACTER STRING*255

      CALL ENTERS('NENXI_1D',*9999)

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        DO i=-1,1
          DO j=0,NEIM
            NXI(i,j,ne)=0
          ENDDO !j
        ENDDO !i
      ENDDO !noelem (ne)
      DO noelem=1,NEELEM(0,nr)
         ne=NEELEM(noelem,nr)
         np1=NPNE(1,nb,ne)
         np2=NPNE(2,nb,ne)
         !Find elements in -Xi1 direction
         NE_SURROUND=NENP(np1,0,nr) !# of elements surrounding np1
         DO noelem0=1,NE_SURROUND !for each neighbouring element
            ne0=NENP(np1,noelem0,nr) !ne# of neighbouring element
            np0=NPNE(2,nb,ne0)  !end np# of neighbouring element
            !If np0 = np1, then will be a parent
            IF(np0.EQ.np1)THEN
               NXI(-1,0,ne)=NXI(-1,0,ne)+1 !count # of parents
               NXI(-1,NXI(-1,0,ne),ne)=ne0 !store parent
            ENDIF               !np0
         ENDDO                  !noelem0
         !Find elements in Xi1 direction
         NE_SURROUND=NENP(np2,0,nr)
         DO noelem0=1,NE_SURROUND
            ne0=NENP(np2,noelem0,nr)
            np0=NPNE(1,nb,ne0)
            IF(np0.EQ.np2)THEN
              NXI(1,0,ne)=NXI(1,0,ne)+1
              WRITE(STRING,'(''>>Increase NEIM. Try '',I4)') NXI(1,0,ne)
              CALL ASSERT(NXI(1,0,ne).LE.NEIM,STRING,ERROR,*9999)
              NXI(1,NXI(1,0,ne),ne)=ne0
            ENDIF               !np0
         ENDDO                  !noelem0
      ENDDO                     !noelem

      CALL EXITS('NENXI_1D')
      RETURN
 9999 CALL ERRORS('NENXI_1D',ERROR)
      CALL EXITS('NENXI_1D')
      RETURN 1
      END


