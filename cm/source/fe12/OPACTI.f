      SUBROUTINE OPACTI(NBJ,NEELEM,nr,FULL_OUTPUT,FEXT,ERROR,*)

C#### Subroutine: OPACTI
C###  Description:
C###    OPACTI outputs active muscle model parameters.

      IMPLICIT NONE
      INCLUDE 'acti00.cmn'
      INCLUDE 'acti01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp50.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),nr
      REAL*8 FEXT(NIFEXTM,NGM,NEM)
      CHARACTER ERROR*(*)
      LOGICAL FULL_OUTPUT
!     Local Variables
      INTEGER i,n,ne,ng,noelem

      CALL ENTERS('OPACTI',*9999)

      WRITE(OP_STRING,'(/'' Region '',I2)') nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(KTYP59(nr).EQ.1) THEN    !SS tension-length-Ca relation  (set Cai)
        WRITE(OP_STRING,'('' SS tension-length-Ca relation'','
     '    //'/''  Max isometric tension at ext.ratio=1 (Tref) = '','
     '    //'D11.4,'' kPa'','
     '    //'/''  Non-dimensional slope parameter (beta)      = '','
     '    //'D11.4,'
     '    //'/''  c50 for [Ca]i saturation curve (0<c<1)      = '','
     '    //'D11.4,'
     '    //'/''  Hill coeff. for [Ca]i saturation curve (h)  = '','
     '    //'D11.4,'
     '    //'/''  Maximum [Ca]i (Ca_max)                      = '','
     '    //'D11.4,'
     '    //'/''  Current [Ca]i at Gauss pt 1 in element 1    = '','
     '    //'D11.4)')
     '      Tref,T0_beta,Ca_c50,Ca_h,Ca_max,FEXT(4,1,1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP59(nr).EQ.2) THEN    !Dynamic HMT (from grid problem)

C       No active material properties are defined here as they are
C       set up in IPCELL for HMT.
        WRITE(OP_STRING,'('' Active material properties are stored'
     '    //' in HMT cell/grid problem: use "fem list cell ..."'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(KTYP59(nr).EQ.3) THEN  !unused

      ELSE IF(KTYP59(nr).EQ.4) THEN  !SS TCa-length relation (set TCa)

        WRITE(OP_STRING,'('' SS tension-length-Ca relation'','
     '    //'/''  Non-dimensional slope parameter (beta)      = '','
     '    //'D11.4,'
     '    //'/''  Current TCa at Gauss pt 1 in element 1      = '','
     '    //'D11.4)')
     '      T0_beta,FEXT(4,1,1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C old MPN 6Aug2014: hair cell code not used
C      ELSE IF(KTYP59(nr).EQ.4) THEN  !Inner hair cell
C        FORMAT='(''   Outer hair cell active stiffness:'''
C     '    //'/'' Resonant frequency  = '',E11.4,'' rad'','
c     '    //'/'' Damping factor      = '',E11.4,'
C     '    //'/'' Phase shift         = '',E11.4,'' rad'','
C     '    //'/'' Amplitude in x dir. = '',E11.4,'' kPa/m'','
C     '    //'/'' Amplitude in y dir. = '',E11.4,'' kPa/m'','
c     '    //'/'' Element numbers     = '',10I5)'
C        WRITE(OP_STRING,FORMAT) ACTIVE_FREQUENCY,ACTIVE_DAMPING,
C     '    ACTIVE_PHASE,ACTIVE_AMPLITUDE_1,ACTIVE_AMPLITUDE_2,
C     '    (NE_ACTI(n),n=1,NE_ACTI(0))
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ENDIF

      IF(FULL_OUTPUT) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO ng=1,NGT(NBJ(1,ne))
            FORMAT='('' FEXT(i,ng='',I3,'',ne='',I4,''): '',8E11.4)'
            WRITE(OP_STRING,FORMAT) ng,ne,(FEXT(i,ng,ne),i=1,8)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !ng
        ENDDO !noelem
      ENDIF !full_output

      CALL EXITS('OPACTI')
      RETURN
 9999 CALL ERRORS('OPACTI',ERROR)
      CALL EXITS('OPACTI')
      RETURN 1
      END


