      SUBROUTINE DRALIG(ISALIG,ISEG,CSEG,STRING,ERROR,*)

C#### Subroutine: DRALIG
C###  Description:
C###    DRALIG draws alignment segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'alig00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISALIG(NWM),ISEG(*)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYMARKER,iw,IWK(6),N3CO,noiw,NTIW
      REAL*8 RFROMC
      LOGICAL CBBREV

      CALL ENTERS('DRALIG',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)


C---------------------------------------------------------------------

C#### Command: FEM draw alignment
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation (GX window) to draw the
C###    points on.
C###  Parameter:    <x=X_SPACING#[0.1]>
C###    Specify the spacing of the points in the x direction.
C###  Parameter:    <y=Y_SPACING#[0.1]>
C###    Specify the spacing of the points in the y direction.
C###  Parameter:    <z=Z_SPACING#[0.1]>
C###    Specify the spacing of the points in the z direction.
C###  Description:
C###    Draw an array of reference points on the specified workstation.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<x=X_SPACING#[0.1]>'
        OP_STRING(4)=BLANK(1:15)//'<y=Y_SPACING#[0.1]>'
        OP_STRING(5)=BLANK(1:15)//'<z=Z_SPACING#[0.1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRALIG',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'X',1,noco+1,NTCO,N3CO)) THEN
          DX_ALIG=RFROMC(CO(N3CO+1))
        ELSE
          DX_ALIG=0.1D0
        ENDIF
        IF(CBBREV(CO,'Y',1,noco+1,NTCO,N3CO)) THEN
          DY_ALIG=RFROMC(CO(N3CO+1))
        ELSE
          DY_ALIG=0.1D0
        ENDIF
C AJP 7/10/95 DZ_ALIGN never used in CMISS
c        IF(CBBREV(CO,'Z',1,noco+1,NTCO,N3CO)) THEN
c          DZ_ALIG=RFROMC(CO(N3CO+1))
c        ELSE
c          DZ_ALIG=0.1D0
c        ENDIF

        INDEX=INDEX_POLYMARKER(0,'POINT','SIZE1','BLACK')

        ALIGNMENT_ON=.TRUE.
        DO noiw=1,NTIW
          IW=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          CALL SGALIG(INDEX,ISALIG(iw),ISEG,iw,CSEG,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO

      ENDIF

      CALL EXITS('DRALIG')
      RETURN
 9999 CALL ERRORS('DRALIG',ERROR)
      CALL EXITS('DRALIG')
      RETURN 1
      END


