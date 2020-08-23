      SUBROUTINE UPFGELEMS(FROM,TO,FG,FGNK,FGNV,CONST,INDICES,
     '  NELISTL,NUMVALUES,POS,SCALE,XAB,ERROR,*)

C#### Subroutine: UPFGELEMS
C###  Description:
C###    This subroutine copies values and derivatives at node points
C###    between FG and XP/YP/CP


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NUMVALUES,FGNK(2,0:NUMVALUES),FGNV(2,0:NUMVALUES),
     '  INDICES(10,2),NELISTL(0:NEM),POS
      REAL*8  CONST,FG(NKM,NVM,NUMVALUES),SCALE,
     '  XAB(NORM,NEM)
      CHARACTER FROM*(*),TO*(*),ERROR*(*)
!      LOGICAL
!     Local Variables
      INTEGER nj,nk,nk1,nk2,nkc,noelem,ne,nv
C currently not used      
c      INTEGER nh,niy,nm,nr,nx
!      REAL*8
!      LOGICAL
!      CHARACTER

      CALL ENTERS('UPFGELEMS',*9999)
c      nm=INDICES(1,POS)
c      niy=INDICES(2,POS)
c      nh=INDICES(3,POS)
      nj=INDICES(4,POS)
c      nx=INDICES(5,POS)
c      nr=INDICES(6,POS)
      
      IF(TO(1:2).EQ.'FG')THEN
        IF(FROM(1:8).EQ.'CONSTANT')THEN
          FGNK(POS,0)=NELISTL(0)
          DO noelem=1,NELISTL(0)
            ne=NELISTL(noelem)
            FGNK(POS,noelem)=NKM
            FGNV(POS,noelem)=FGNV(1,noelem)
            DO nv=1,FGNV(POS,noelem)
              FG(1,nv,noelem)=SCALE*CONST
              nkc=0
              DO nk=2,FGNK(2,noelem)
                nkc=nkc+1
                FG(nk,nv,noelem)=0.0D0
              ENDDO
            ENDDO
          ENDDO
         ELSEIF(FROM(1:3).EQ.'XAB')THEN
          FGNK(POS,0)=NELISTL(0)
          DO noelem=1,NELISTL(0)
            ne=NELISTL(noelem)
            IF(INDICES(7,POS).EQ.0)THEN
              nk1=1
              nk2=1
              FGNK(POS,noelem)=1
            ELSE
              nk1=INDICES(7,POS)
              nk2=INDICES(7,POS)
              FGNK(POS,noelem)=1
            ENDIF
            FGNV(POS,noelem)=1
            DO nv=1,1
              nkc=0
              DO nk=nk1,nk2
                nkc=nkc+1
                FG(nkc,nv,noelem)=SCALE*XAB(nj,ne)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSEIF(FROM(1:2).EQ.'FG')THEN
        IF(TO(1:3).EQ.'XAB')THEN
          DO noelem=1,NELISTL(0)
            ne=NELISTL(noelem)
            IF(INDICES(7,POS).EQ.0)THEN
              nk1=1
              nk2=1
              FGNK(POS,noelem)=1
            ELSE
              nk1=INDICES(7,POS)
              nk2=INDICES(7,POS)
              FGNK(POS,noelem)=1
            ENDIF
            nkc=0
            DO nk=nk1,nk2
              nkc=nkc+1
              DO nv=1,1
                XAB(nj,ne)=FG(nkc,nv,noelem)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('UPFGELEMS')
      RETURN
 9999 CALL ERRORS('UPFGELEMS',ERROR)
      CALL EXITS('UPFGELEMS')
      RETURN 1
      END


