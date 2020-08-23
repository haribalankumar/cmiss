      SUBROUTINE SGPLOT(INDEX,ISEG,ISPLOT,iw,nb,ne,nh,
     '  XE,ZE,CSEG,ERROR,*)

C#### Subroutine: SGPLOT
C###  Description:
C###    SGPLOT creates element 3D Phigs plots.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISPLOT,iw,nb,ne,nh
      REAL*8 XE(NSM,NJM),ZE(NSM,NHM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER INDEX_OLD,nj,nn,NPTS(1)
      REAL*8 PTS(3,16)

      CALL ENTERS('SGPLOT',*9999)
      CALL OPEN_SEGMENT(ISPLOT,ISEG,iw,'PLOT',INDEX,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

      NPTS(1)=4
      DO nj=1,2
        PTS(nj,1)=XE(1,nj)
        PTS(nj,2)=XE(2,nj)
        PTS(nj,3)=XE(4,nj)
        PTS(nj,4)=XE(3,nj)
      ENDDO
      PTS(3,1)=ZE(1,nh)
      PTS(3,2)=ZE(2,nh)
      PTS(3,3)=ZE(4,nh)
      PTS(3,4)=ZE(3,nh)
      IF(DOP) THEN
        DO nj=1,3
          WRITE(OP_STRING,'('' PTS('',I1,'',nn):'',9E10.3)')
     '      nj,(PTS(nj,nn),nn=1,NNT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF
      CALL FILL_AREA(1,iw,NPTS(1),PTS,ERROR,*9999)

      CALL CLOSE_SEGMENT(ISPLOT,iw,ERROR,*9999)
      CALL EXITS('SGPLOT')
      RETURN
 9999 CALL ERRORS('SGPLOT',ERROR)
      CALL EXITS('SGPLOT')
      RETURN 1
      END


