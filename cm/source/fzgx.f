C**** GX Graphics Library dummy routines

      SUBROUTINE gxwait(A,ERR)
      INTEGER ERR
      DIMENSION A(*)
      CHARACTER ERROR*10

C     Do nothing here
C     Set error code so that this will not crash if there
C       is no window present
      ERR=0

      RETURN
      END

      SUBROUTINE CHECK_GX_OPEN(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub CHECK_GX_OPEN')
      RETURN 1
      END

      SUBROUTINE DETECT_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub DETECT_GX')
      RETURN 1
      END

      SUBROUTINE FILL_AREA_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub FILL_AREA_GX')
      RETURN 1
      END

      SUBROUTINE FILL_LINE_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub FILL_LINE_GX')
      RETURN 1
      END

      SUBROUTINE LOCATOR_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub LOCATOR_GX')
      RETURN 1
      END

      SUBROUTINE OPEN_SEGMENT_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub OPEN_SEGMENT_GX')
      RETURN 1
      END

      SUBROUTINE POLYLINE_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub POLYLINE_GX')
      RETURN 1
      END

      SUBROUTINE POLYMARKER_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub POLYMARKER_GX')
      RETURN 1
      END

      SUBROUTINE TEXT_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub TEXT_GX')
      RETURN 1
      END

      SUBROUTINE VECTOR_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub VECTOR_GX')
      RETURN 1
      END

      SUBROUTINE VISIB_GX(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module fegx.f: '//
     '  'need sub VISIB_GX')
      RETURN 1
      END

      SUBROUTINE gxSAVE_WINDOW_PS(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub gxSAVE_WINDOW_PS')
      RETURN
      END

      SUBROUTINE setwin(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub setwin')
      RETURN
      END

      SUBROUTINE clobjt(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub clobjt')
      RETURN
      END

      SUBROUTINE clwin(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub clwin')
      RETURN
      END

      SUBROUTINE unswin(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub unswin')
      RETURN
      END

      SUBROUTINE dlobjt(A)
      DIMENSION A(*)
      CHARACTER ERROR*10
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub dlobjt')
      RETURN
      END

      SUBROUTINE fpick(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub fpick')
      RETURN
      END

      SUBROUTINE clgrph(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub clgrph')
      RETURN
      END

      SUBROUTINE opwin(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub opwin')
      RETURN
      END

      SUBROUTINE setview(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub setview')
      RETURN
      END

      SUBROUTINE SAVE_WINDOW_PPM(A)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with GX Graphics library: '//
     '  'need sub SAVE_WINDOW_PPM')
      RETURN
      END
