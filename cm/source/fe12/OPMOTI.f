      SUBROUTINE OPMOTI(IBT,NBH,NEELEM,NHP,NKH,NPNODE,nr,FIX,YP,
     '  ERROR,*)

C#### Subroutine: OPMOTI
C###  Description:
C###    OPMOTI outputs motion parameters for region nr.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'four00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'moti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHP(NPM),NKH(NHM,NPM,NCM),NPNODE(0:NP_R_M,0:NRM),
     '  nr
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER IBEG1,IEND1,nb,nc,nh,nhh2,nhx,nk,NO_COEFFS,nonode,np,
     '  NT_COEFFS,nx,NYTOT
      CHARACTER CHAR1*5,TITLE1(2)*20

      DATA TITLE1/'Fourier coefficients',
     '            'Spreadsheet column'/

      CALL ENTERS('OPMOTI',*9999)

      nx=1 !temporary
      IF(KTYP58(nr).LE.2) THEN !Motion type is Geom coords or Displ's
        WRITE(OP_STRING,'('' Motion type is '',A)')
     '    TITLE1(MOTION_TYPE)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(MOTION_TYPE.EQ.1) THEN !Fourier coefficients
          WRITE(OP_STRING,'('' Ang. freq of fundamental = '',E11.4)')
     '      OMEGA
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nhh2=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NJ_LOC(NJL_GEOM,nhh2,nr)
            nc=1 !Temporary AJP 17-12-91
            nb=NBH(nh,nc,NEELEM(1,1))
            WRITE(OP_STRING,'('' For geometric variable'',I2)') nh
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''   Fourier basis number    = '',I2)')
     '        nb
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''   Number of harmonics     = '',I2)')
     '        (IBT(2,NIT(nb),nb)-1)/2
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

        NYTOT=0
        nc=1 !Temporary AJP 17-12-91
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          WRITE(CHAR1,'(I5)') np
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            nc=1 !Temporary AJP 17-12-91
            nb=NBH(nh,nc,NEELEM(1,nr))
            NT_COEFFS=IBT(2,NIT(nb),nb)
            DO nk=1,NKH(nh,np,nc)
              CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
              IF(FIX(NYTOT+1,5).OR.KTYP8.EQ.5) THEN
                FORMAT='(/'' Node '//CHAR1(IBEG1:IEND1)//' Variable '','
     '            //'I1,'//''' Derivative '',I1,'' Fourier coeffs: '')'
                WRITE(OP_STRING,FORMAT) nh,nk
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(11D12.4)')
     '            (YP(NYTOT+(NO_COEFFS-1)*NKH(nh,np,nc)+nk,MOTION_IY),
     '             NO_COEFFS=1,NT_COEFFS)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            NYTOT=NYTOT+NKH(nh,np,nc)
          ENDDO
        ENDDO

      ELSE IF(KTYP58(nr).EQ.3) THEN !Motion type is Flow
        WRITE(OP_STRING,'('' Motion type is '',A,'' for flow at '','
     '    //'''node '',I4)') TITLE1(MOTION_TYPE),NP_MOTION
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(MOTION_TYPE.EQ.1) THEN
          nb=NB_MOTION
          NT_COEFFS=IBT(2,NIT(nb),nb)
          WRITE(OP_STRING,'('' Fourier basis number     = '',I2)') nb
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Number of harmonics      = '',I2)')
     '      (IBT(2,NIT(nb),nb)-1)/2
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Ang. freq of fundamental = '',E11.4)')
     '      OMEGA
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Fourier coeffs: '',9E12.3)')
     '      (FLOW_COEFFS(NO_COEFFS),NO_COEFFS=1,NT_COEFFS)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

      ELSE IF(KTYP58(nr).EQ.4) THEN !Motion type is Lung gas flow
        WRITE(OP_STRING,
     '    '('' Motion type is lung gas flow at node '',I4)') NP_MOTION
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Number of breaths       = '',I2)') NT_CYCLES
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Lung tidal volume       = '',E11.4)') FLOW_COEFFS(1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Duration of inspiration = '',E11.4)') FLOW_COEFFS(2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Duration of breath-hold = '',E11.4)') FLOW_COEFFS(3)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '('' Duration of expiration  = '',E11.4)') FLOW_COEFFS(4)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ENDIF

      CALL EXITS('OPMOTI')
      RETURN
 9999 CALL ERRORS('OPMOTI',ERROR)
      CALL EXITS('OPMOTI')
      RETURN 1
      END


