      SUBROUTINE UPFGDATA(FROM,TO,FG,FGNK,FGNV,CONST,INDICES,NUMVALUES,
     &  POS,SCALE,ZD,ERROR,*)

C#### Subroutine: UPFGDATA
C###  Description:
C###    This subroutine copies values and derivatives at data points
C###    between FG and ZD


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NUMVALUES,FGNK(2,0:NUMVALUES),FGNV(2,0:NUMVALUES),
     '  INDICES(10,2),POS
      REAL*8  CONST,FG(NKM,NVM,NUMVALUES),SCALE,ZD(NJM,NDM)
      CHARACTER ERROR*(*),FROM*(*),TO*(*)
!     Local Variables
      INTEGER nd,nj

      CALL ENTERS('UPFGDATA',*9999)

      nj=INDICES(4,POS)
c      nr=INDICES(6,POS)
      IF(TO(1:2).EQ.'FG')THEN
        IF(FROM(1:8).EQ.'CONSTANT')THEN
          DO nd=1,NDT
            FGNK(POS,nd)=1
            FGNV(POS,nd)=1
            FG(1,1,nd)=SCALE*CONST
          ENDDO
          FGNK(POS,0)=NDT
        ELSEIF(FROM(1:2).EQ.'ZD')THEN
          DO nd=1,NDT
            FGNK(POS,nd)=1
            FGNV(POS,nd)=1
            FG(1,1,nd)=SCALE*ZD(nj,nd)
          ENDDO
          FGNK(POS,0)=NDT
        ELSE
          ERROR='>> Not Implemented'
          GOTO 9999
        ENDIF
      ELSEIF(FROM(1:2).EQ.'FG')THEN
        IF(TO(1:2).EQ.'ZD')THEN
          DO nd=1,NDT
            ZD(nj,nd)=FG(1,1,nd)
          ENDDO
        ELSE
          ERROR='>> Not Implemented'
          GOTO 9999
        ENDIF
      ENDIF

      CALL EXITS('UPFGDATA')
      RETURN
 9999 CALL ERRORS('UPFGDATA',ERROR)
      CALL EXITS('UPFGDATA')
      RETURN 1
      END


