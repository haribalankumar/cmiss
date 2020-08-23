      SUBROUTINE CHECK_GX_OPEN(ERROR,*)

C#### Subroutine: CHECK_GX_OPEN
C###  Description:
C###    CHECK_GX_OPEN opens GX if not already open.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
      INCLUDE 'disp00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'gks001.cmn'
      INCLUDE 'gx.inc'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ERR

      CALL ENTERS('CHECK_GX_OPEN',*9999)

      IF(.NOT.GKS) THEN
        IF(DOP) THEN
          WRITE(OP_STRING,*)' Setup opening GX'
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(WINDOW_TYPE(1:5).EQ.'MOTIF') THEN
          CALL OPGRPH(ERR)
          IF(err.ne.0) THEN
            ERROR='>>Error from GX'
            GOTO 9999
          ENDIF
          COLOUR_WS=.TRUE.
          MAXCOLOURS=240
C CPB 8/3/94 Not sure about these
          NINDICES = gxNSPECT-gxSPOFF
        ENDIF
        GKS=.TRUE.
        CALL SET_COLOUR_LUT(COLOUR_LUT,ERROR,*9999)
      ENDIF

      CALL EXITS('CHECK_GX_OPEN')
      RETURN
 9999 CALL ERRORS('CHECK_GX_OPEN',ERROR)
      CALL EXITS('CHECK_GX_OPEN')
      RETURN 1
      END


