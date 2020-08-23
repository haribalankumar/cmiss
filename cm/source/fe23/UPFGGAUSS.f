      SUBROUTINE UPFGGAUSS(FROM,TO,FG,FGNK,FGNV,CONST,INDICES,
     '      NBH,NEELEM,NUMVALUES,POS,SCALE,YG,ERROR,*)

C#### Subroutine: UPFGGAUSS
C###  Description:
C###    This subroutine copies values and derivatives at node points
C###    between FG and XG/YG/CG


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NUMVALUES,FGNK(2,0:NUMVALUES),FGNV(2,0:NUMVALUES),
     '  INDICES(10,2),NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),POS
      REAL*8  CONST,FG(NKM,NVM,NUMVALUES),SCALE,YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*),FROM*(*),TO*(*)
!      LOGICAL
!     Local Variables
      INTEGER count,ne,ng,nh,niy,noelem,nr
!      INTEGER nj,nm,nx
!      REAL*8
!      LOGICAL
!      CHARACTER

      CALL ENTERS('UPFGGAUSS',*9999)

C      nm=INDICES(1,POS)
      niy=INDICES(2,POS)
      nh=INDICES(3,POS)
C      nj=INDICES(4,POS)
C      nx=INDICES(5,POS)
      nr=INDICES(6,POS)
      IF(TO(1:2).EQ.'FG')THEN
        IF(FROM(1:8).EQ.'CONSTANT')THEN
          count=0
          nh=INDICES(3,1)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ng=1,NGT(NBH(nh,1,ne))
              count=count+1
              FGNK(POS,count)=1
              FGNV(POS,count)=1
              FG(1,1,count)=SCALE*CONST
            ENDDO
          ENDDO
          FGNK(POS,0)=count
        ELSEIF(FROM(1:2).EQ.'YG')THEN
          count=0
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ng=1,NGT(NBH(nh,1,ne))
              count=count+1
              FGNK(POS,count)=1
              FGNV(POS,count)=1
              FG(1,1,count)=SCALE*YG(niy,ng,ne)
            ENDDO
          ENDDO
          FGNK(POS,0)=count
        ELSE
          ERROR='>> Not Implemented'
          GOTO 9999
        ENDIF
      ELSEIF(FROM(1:2).EQ.'FG')THEN
        IF(TO(1:2).EQ.'YG')THEN
          count=0
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ng=1,NGT(NBH(nh,1,ne))
              count=count+1
              YG(niy,ng,ne)=FG(1,1,count)
            ENDDO
          ENDDO
        ELSE
          ERROR='>> Not Implemented'
          GOTO 9999
        ENDIF
      ENDIF

      CALL EXITS('UPFGGAUSS')
      RETURN
 9999 CALL ERRORS('UPFGGAUSS',ERROR)
      CALL EXITS('UPFGGAUSS')
      RETURN 1
      END


