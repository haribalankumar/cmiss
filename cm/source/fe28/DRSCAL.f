      SUBROUTINE DRSCAL(ISEG,ISSCAL,CSEG,STRING,ERROR,*)

C#### Subroutine: DRSCAL
C###  Description:
C###    DRSCAL draws scale segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gks001.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSCAL(NWM,NGRSEGM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IBEG3,IEND,IEND1,IEND2,IEND3,INDEX,
     '  INDEX_TEXT,iw,IWK(6),N3CO,noiw,NTIW
      REAL*8 RFROMC
      CHARACTER CHAR1*11,CHAR2*11,CHAR3*9
      LOGICAL CBBREV,COLOUR

      CALL ENTERS('DRSCAL',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C        CHAR1=CFROMR(ZMINI,'(E11.3)')
        WRITE(CHAR1,'(E11.3)') ZMINI
        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
C        CHAR2=CFROMR(ZMAXI,'(E11.3)')
        WRITE(CHAR2,'(E11.3)') ZMAXI
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        IF(NINDICES.GE.150) THEN
          CHAR3='colour'
        ELSE
          CHAR3='greyscale'
        ENDIF
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)

C---------------------------------------------------------------------

C#### Command: FEM draw scale
C###  Parameter:    <(greyscale/colour)[greyscale]>
C###    Specify whether the scale is to be colour or greyscale.
C###  Parameter:    <z=ZMIN#[0.0] to ZMAX#[1.0]>
C###    Specify the upper and lower bounds of the scale.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Description:
C###    Draws a scale segment, in either colour or greyscale. The scale
C###    starts at ZMIN and goes up to ZMAX. The scale is drawn on the
C###    specified workstation.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)
     '    //'<(greyscale/colour)['//CHAR3(IBEG3:IEND3)//']>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<z=ZMIN#['//CHAR1(IBEG1:IEND1)//' to ZMAX#['
     '    //CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(5)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRSCAL',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_TEXT(0,'WIDTH1','FONT1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')
        ENDIF
        IF(CBBREV(CO,'COLOUR',1,noco+1,NTCO,N3CO)) THEN
          COLOUR=.TRUE.
        ELSE
          IF(NINDICES.GE.150) THEN
            COLOUR=.TRUE.
          ELSE
            COLOUR=.FALSE.
          ENDIF
        ENDIF
        IF(CBBREV(CO,'Z',1,noco+1,NTCO,N3CO)) THEN
          ZMINI=RFROMC(CO(N3CO+1))
        ENDIF
        IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
          ZMAXI=RFROMC(CO(N3CO+1))
        ENDIF

        IF(ADD) THEN
          NTSCAL=NTSCAL+1
        ELSE IF(NTSCAL.EQ.0) THEN
          NTSCAL=1
        ENDIF
        CALL ASSERT(NTSCAL.LE.NRM,'>>NRM too small',ERROR,*9999)
        ZDIFF=ZMAXI-ZMINI
        DO noiw=1,NTIW
          IW=IWK(noiw)
C CPB 10/9/92 Changed workstation update mode
          CALL ACWK(iw,1,ERROR,*9999)
          CALL SGSCAL(INDEX,ISEG,ISSCAL(iw,NTSCAL),iw,NTSCAL,COLOUR,
     '      CSEG,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('DRSCAL')
      RETURN
 9999 CALL ERRORS('DRSCAL',ERROR)
      CALL EXITS('DRSCAL')
      RETURN 1
      END


