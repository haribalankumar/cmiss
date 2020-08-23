      SUBROUTINE SGOBJE(INDEX,ISEG,ISOBJE,iw,ND1,ND2,noobje,nopart,
     '  CSEG,ZDD,ERROR,*)

C#### Subroutine: SGOBJE
C###  Description:
C###    Create object part segments ISOBJE(iw,noobje,nopart).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISOBJE,iw,ND1,ND2,noobje,nopart
      REAL*8 ZDD(NJM,NDM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,INDEX_OLD
c      REAL*8 X(3)
      CHARACTER CHAR2*2

      CALL ENTERS('SGOBJE',*9999)
      WRITE(CHAR2,'(I2)') NOOBJE
      CALL OPEN_SEGMENT(ISOBJE,ISEG,iw,'OBJECT '//CHAR2(1:2),INDEX,
     '  INDEX_OLD,NOPART,1,CSEG,ERROR,*9999)

      WRITE(CHAR2,'(I2)') NOPART
      CALL STRING_TRIM(CHAR2,IBEG,IEND)
c      X(1)=XWC
c      X(2)=YWC
c      X(3)=ZWC
C      CALL TEXT(INDEX,iw,CHAR2(IBEG:IEND),X,ERROR,*9999)
      CALL POLYMARKER(INDEX,iw,ND2-ND1+1,ZDD(1,ND1),ERROR,*9999)

      CALL CLOSE_SEGMENT(ISOBJE,iw,ERROR,*9999)

      CALL EXITS('SGOBJE')
      RETURN
 9999 CALL ERRORS('SGOBJE',ERROR)
      CALL EXITS('SGOBJE')
      RETURN 1
      END


