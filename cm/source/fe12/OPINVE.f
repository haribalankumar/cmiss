      SUBROUTINE OPINVE(ERROR,*)

C#### Subroutine: OPINVE
C###  Description:
C###    OPINVE outputs regularisation options for inverting
C###    the transfer matrix.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      CHARACTER ERROR*(*)

C*** adding output for inverse approaches

      CALL ENTERS('OPINVE',*9999)

C*** Potential Inverse
      IF(INV_APPROACH.EQ.1) THEN

        WRITE(OP_STRING,'('' Using a potential inverse approach '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

CC Modified JMB 14-FEB-2000
        IF(ICALC_TRANSFER.EQ.1) THEN
          WRITE(OP_STRING,'('' The inverse transfer matrix is '','
     '      //'''calculated explicitly'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(IREGULARISE.EQ.1) THEN
            WRITE(OP_STRING,'('' No stabilisation scheme is used'','
     '        //''' in the inversion '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(IREGULARISE.EQ.2) THEN
            WRITE(OP_STRING,'('' SVD is used to stabilise the '','
     '        //'''inverse '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' SVD cutoff is '',D12.4)')RMIN_SVD
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(IREGULARISE.EQ.3) THEN
            WRITE(OP_STRING,'('' Twomey regularisation is used '','
     '        //'''to stabilise the inverse'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(IREGULARISE.EQ.4) THEN
            WRITE(OP_STRING,'('' Zero-order Tikhonov '','
     '        //'''regularisation is used to stabilise the '','
     '        //'''inverse'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(IREGULARISE.EQ.5) THEN
            WRITE(OP_STRING,'('' First-order Tikhonov '','
     '        //'''regularisation is used to stabilise the '','
     '        //'''inverse'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(IREGULARISE.EQ.6) THEN
            WRITE(OP_STRING,'('' Second-order Tikhonov  '','
     '        //'''regularisation is used to stabilise the '','
     '        //'''inverse'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(IREGULARISE.EQ.7) THEN
            WRITE(OP_STRING,'('' Local regularisation '','
     '        //'''is used to stabilise the inverse'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(IREGULARISE.GE.4.AND.IREGULARISE.LE.6) THEN
            WRITE(OP_STRING,'('' The Tikhonov parameter is '',D12.4)')
     '        TIKH_VALUE
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

CC New JMB 14-FEB-2000
        ELSEIF(ICALC_TRANSFER.EQ.2) THEN
          WRITE(OP_STRING,'('' The inverse solution is calculated'','
     '      //''' from matrix equations'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ISTABILISE.EQ.1) THEN
            WRITE(OP_STRING,'('' No stabilisation scheme is used'','
     '        //''' in the inversion '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(ISTABILISE.EQ.2) THEN
            WRITE(OP_STRING,'('' Truncated SVD scheme is used to'','
     '        //''' stabilise the inverse '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(ISTABILISE.EQ.3) THEN
            WRITE(OP_STRING,'('' Tikhonov regularisation is used '','
     '        //'''to stabilise the inverse '')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            ERROR='Unknown stabilising scheme'
            GOTO 9999
          ENDIF
          IF(ISTABILISE.NE.1) THEN
            IF(ICONSTRAINT.EQ.1) THEN ! Additional constaints
              WRITE(OP_STRING,'('' Additional contraint is       :'','
     '          //''' none'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(ICONSTRAINT.EQ.2) THEN
              WRITE(OP_STRING,'('' Additional contraint is       :'','
     '          //''' a surface gradient'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(ICONSTRAINT.EQ.3) THEN
              WRITE(OP_STRING,'('' Additional contraint is       :'','
     '          //''' a surface Laplacian'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE
              ERROR='Unknown additional constraints'
              GOTO 9999
            ENDIF
            IF(ICOUPLING.EQ.1) THEN ! Additional coupling
              WRITE(OP_STRING,'('' Additional coupling is        :'','
     '          //''' none'')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(ICOUPLING.EQ.2) THEN
              WRITE(OP_STRING,'('' Additional coupling is        :'','
     '          //''' Greensite'')')
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE
              ERROR='Unknown additional coupling'
              GOTO 9999
            ENDIF
            IF(IREGULARISE.EQ.1) THEN
              WRITE(OP_STRING,'('' GCV criterion used for'','
     '          //''' determining regularisation parameters'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(IREGULARISE.EQ.2) THEN
              WRITE(OP_STRING,'('' L-curve criterion used for'','
     '          //''' determining regularisation parameters'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(IREGULARISE.EQ.3) THEN
              WRITE(OP_STRING,'('' Picard criterion used'','
     '          //''' for determining regularisation parameters'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(IREGULARISE.EQ.4) THEN
              WRITE(OP_STRING,'('' Quasi-optimality criterion used'','
     '          //''' for determining regularisation parameters'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(IREGULARISE.EQ.5) THEN
              WRITE(OP_STRING,'('' Optimal criterion used for'','
     '          //''' determining regularisation parameters'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(IREGULARISE.EQ.6) THEN
              WRITE(OP_STRING,'('' CRESO criterion used for'','
     '          //''' determining regularisation parameters'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSEIF(IREGULARISE.EQ.7) THEN
              WRITE(OP_STRING,'('' Zero-crossing criterion used '','
     '          //'''for determining regularisation parameters'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDIF
C*** Activation inverse
      ELSEIF(INV_APPROACH.EQ.2) THEN

        WRITE(OP_STRING,'('' Using an activation inverse '','
     '    //'''approach '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(ICALC_TRANSFER.EQ.1) THEN ! Imaging technique
          WRITE(OP_STRING,'('' Imaging technique              :'','
     '      //''' zero-crossing'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          ERROR='Unknown imaging approach'
          GOTO 9999
        ENDIF

      ELSE
          ERROR='Unknown Inverse Approach'
          GOTO 9999
      ENDIF !INV_APPROACH

      CALL EXITS('OPINVE')
      RETURN
 9999 CALL ERRORS('OPINVE',ERROR)
      CALL EXITS('OPINVE')
      RETURN 1
      END


