      SUBROUTINE D2XRCDXI(NITB,nr,DXRCXI,D2XRCXI,XG,ERROR,*)

C#### Subroutine: D2XRCDXI
C###  Description:
C###    D2XRCDXI calculates DXRCX(njj,ni) = dXRC(njj)/dXI(ni) and
C###    D2XRCX(njj,mi,ni) = d2XRC(njj)/dXI(mi)dXI(ni) from XG, where
C###    XRC(njj) are rectangular Cartesian and XI(ni) are local element
C###    coordinates.  If NJT < 3 then the unused njj elements are
C###    initiallized to zero for FIBRE_REF_VECS and D_FIBRE_REF_VECS.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NITB,nr
      REAL*8  DXRCXI(3,3),D2XRCXI(3,3,3),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,ni,nj,njj,NJTR,NU1(3),NU2(3,3)
      REAL*8 DXXI(3,3),D2XXI(3,3,3),X(3)
C      INTEGER njj2
C      REAL*8 D2XRCX(3,3,3),DXRCX(3,3)

      DATA NU1/2,4,7/
      DATA NU2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('D2XRCDXI',*9999)

      NJTR=NJ_LOC(NJL_GEOM,0,nr)

      IF(ITYP10(nr).EQ.1) THEN !ref coords are already rc
        DO njj=1,NJTR
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          DO ni=1,NITB
            DXRCXI(njj,ni)=XG(nj,NU1(ni))
            DO mi=1,NITB
              D2XRCXI(njj,mi,ni)=XG(nj,NU2(mi,ni))
            ENDDO
          ENDDO
        ENDDO !njj
        DO njj=NJTR+1,3
          DO ni=1,NITB
            DXRCXI(njj,ni)=0.0d0 !for FIBRE_REF_VECS
            DO mi=1,NITB
              D2XRCXI(njj,mi,ni)=0.0d0 !for D_FIBRE_REF_VECS
            ENDDO
          ENDDO
        ENDDO
      ELSE !convert ref coords to rc coords
        DO njj=1,NJTR
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          X(njj)=XG(nj,1)
          DO ni=1,NITB
            DXXI(njj,ni)=XG(nj,NU1(ni))
            DO mi=1,NITB
              D2XXI(njj,mi,ni)=XG(nj,NU2(mi,ni))
            ENDDO
          ENDDO
        ENDDO !njj
        DO njj=NJTR+1,3
          X(njj)=0.0d0 !for DXRCDX
        ENDDO
C!!! KAT 28Sep98: This routine is not written yet
C        CALL D2XRCDX(ITYP10(nr),NJTR,DXRCX,D2XRCX,X,ERROR,*9999)
C        CALL DXRCDX(ITYP10(nr),NJTR,DXRCX,X,ERROR,*9999)
C        DO ni=1,NITB
C          DO njj=1,3
C            DXRCXI(njj,ni)=0.0d0
C          ENDDO
C          DO njj=1,NJTR
C            DO njj2=1,NJTR
C              DXRCXI(njj2,ni)=
C     '          DXRCXI(njj2,ni)+DXRCX(njj2,njj)*DXXI(njj,ni)
C            ENDDO
C          ENDDO
C        ENDDO
C!!!    Need some more code in here for D2XRCXI
        ERROR='>>Only implemented for rectangular coordinates'
        GO TO 9999

      ENDIF

      CALL EXITS('D2XRCDXI')
      RETURN
 9999 CALL ERRORS('D2XRCDXI',ERROR)
      CALL EXITS('D2XRCDXI')
      RETURN 1
      END


