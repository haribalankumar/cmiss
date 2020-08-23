      SUBROUTINE DXRCDXI(NITB,nr,DXRCXI,XG,ERROR,*)

C#### Subroutine: DXRCDXI
C###  Description:
C###    DXRCDXI calculates DXRCX(njj,ni) = dXRC(njj)/dXI(ni) from XG,
C###    where XRC(njj) are rectangular Cartesian and XI(ni) are local element
C###    coordinates.  If NJT < 3 then the unused
C###    elements in column ni are initiallized to zero for FIBRE_REF_VECS.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NITB,nr
      REAL*8  DXRCXI(3,3),DXXI(3,3),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,nj,njj,njj2,NJTR,NU1(3)
      REAL*8 DXRCX(3,3),X(3)

      DATA NU1/2,4,7/

      CALL ENTERS('DXRCDXI',*9999)

      NJTR=NJ_LOC(NJL_GEOM,0,nr)

      IF(ITYP10(nr).EQ.1) THEN !ref coords are already rc
        DO njj=1,NJTR
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          DO ni=1,NITB
            DXRCXI(njj,ni)=XG(nj,NU1(ni))
          ENDDO !ni
        ENDDO !njj
        DO njj=NJTR+1,3
          DO ni=1,NITB
            DXRCXI(njj,ni)=0.0d0 !for FIBRE_REF_VECS
          ENDDO !ni
        ENDDO !njj
      ELSE !convert ref coords to rc coords
        DO njj=1,NJTR
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          X(njj)=XG(nj,1)
          DO ni=1,NITB
            DXXI(njj,ni)=XG(nj,NU1(ni))
          ENDDO !ni
        ENDDO !njj
        DO njj=NJTR+1,3
          X(njj)=0.0d0 !for DXRCDX
        ENDDO !njj
        CALL DXRCDX(ITYP10(nr),NJTR,DXRCX,X,ERROR,*9999)
        DO ni=1,NITB
          DO njj=1,3
            DXRCXI(njj,ni)=0.0d0
          ENDDO !njj
          DO njj=1,NJTR
            DO njj2=1,NJTR
              DXRCXI(njj2,ni)=
     '          DXRCXI(njj2,ni)+DXRCX(njj2,njj)*DXXI(njj,ni)
            ENDDO !njj2
          ENDDO !njj
        ENDDO !ni
      ENDIF !ITYP10

      CALL EXITS('DXRCDXI')
      RETURN
 9999 CALL ERRORS('DXRCDXI',ERROR)
      CALL EXITS('DXRCDXI')
      RETURN 1
      END


