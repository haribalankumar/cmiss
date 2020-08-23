      SUBROUTINE MAKESPLINE(DIRECTION,nknots,layer,ZD,
     '  DIAG,UPPERDIAG,P,ERROR,*)

C#### Subroutine: MAKESPLINE
C###  Description:
C###    Make a 1D spline given information about the geometry and
C###    its associated dependent variable.
C###

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'

      INTEGER NKNOTSM
      PARAMETER (NKNOTSM=40)
! Parameters
      REAL*8 ZD(NJM,NDM),
     '  P(NKNOTSM), DIAG(NKNOTSM), UPPERDIAG(NKNOTSM)
      INTEGER nknots,layer
      CHARACTER ERROR*(*), DIRECTION*(*)
! Local Variable
      REAL*8 W(NKNOTSM)
      INTEGER nd,offset

! Function
      REAL*8 DATA_DIST


      CALL ENTERS('MAKESPLINE',*9999)

      IF(DIRECTION(1:10).EQ.'HORIZONTAL') THEN
        offset=(layer-1)*nknots
      ELSEIF (DIRECTION(1:8).EQ.'VERTICAL') THEN
        offset=NDT
      ELSE
        ERROR='>>Unknown Parameter, Update Code '
        GOTO 9999
      ENDIF

      UPPERDIAG(1)=DATA_DIST(offset+1,offset+2,ZD) !first 2 points
      DO nd=2,nknots-1
        DIAG(nd)=2.D0*DATA_DIST(offset+nd+1,offset+nd-1,ZD)
        UPPERDIAG(nd)=DATA_DIST(offset+nd+1,offset+nd,ZD)
      ENDDO

C*** Wrap around to beginning
      IF(DIRECTION(1:10).EQ.'HORIZONTAL') THEN
        DIAG(nknots)=2.D0*DATA_DIST(offset+1,offset+nknots-1,ZD)
        UPPERDIAG(nknots)=DATA_DIST(offset+1,offset+nknots,ZD)
      ENDIF

      DO nd=2,nknots-1
        W(nd)=6.D0*(
     '      (ZD(NJT+1,offset+nd+1)-ZD(NJT+1,offset+nd))/UPPERDIAG(nd)
     '    - (ZD(NJT+1,offset+nd)-ZD(NJT+1,offset+nd-1))/UPPERDIAG(nd-1)
     '    )
      ENDDO

      IF(DIRECTION(1:10).EQ.'HORIZONTAL') THEN
        W(nknots)=6.D0*(
     '      (ZD(NJT+1,offset+1)-ZD(NJT+1,offset+nknots))
     '    /UPPERDIAG(nknots)
     '    - (ZD(NJT+1,offset+nknots)-ZD(NJT+1,offset+nknots-1))
     '    /UPPERDIAG(nknots-1)
     '    )
      ENDIF

      IF(DOP) THEN
        IF(nknots.GT.1) THEN
          IF(DIRECTION(1:10).EQ.'HORIZONTAL') THEN
            WRITE(OP_STRING,'(''  W '',10(E11.3))')
     '        (W(nd),nd=2,nknots)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  D '',10(E11.3))')
     '        (DIAG(nd),nd=2,nknots)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'(''  W '',10(E11.3))')
     '        (W(nd),nd=2,nknots-1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  D '',10(E11.3))')
     '        (DIAG(nd),nd=2,nknots-1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

C*** set the values of the second derivative
      P(1)=0.D0
      P(nknots)=0.D0

C*** Here follows the tri-diag SOLVER ... fixed for a linear cubic
C    NOTE This destroys DIAG !!!
      DO nd=2,nknots-2
        W(nd+1)=W(nd+1)-W(nd)*UPPERDIAG(nd)/DIAG(nd)
        DIAG(nd+1)=DIAG(nd+1)-UPPERDIAG(nd)*UPPERDIAG(nd)/DIAG(nd)
      ENDDO

C LKC 27-FEB-1999 - Changing the loop bounds for the horizontal section
C from nknots to  (nknots-1) - makes it consistent for hori and vert.

      IF(DIRECTION(1:10).EQ.'HORIZONTAL') THEN
        nd=nknots-1
        W(nknots)=W(nknots)-W(nd)*UPPERDIAG(nd)/DIAG(nd)
        DIAG(nknots)=DIAG(nknots)-UPPERDIAG(nd)*UPPERDIAG(nd)/DIAG(nd)
        DO nd=nknots-1,2,-1
          P(nd)=(W(nd)-UPPERDIAG(nd)*P(nd+1))/DIAG(nd)
        ENDDO
      ELSE
        DO nd=nknots-1,2,-1
          P(nd)=(W(nd)-UPPERDIAG(nd)*P(nd+1))/DIAG(nd)
        ENDDO
      ENDIF


      CALL EXITS('MAKESPLINE')
      RETURN
 9999 CALL ERRORS('MAKESPLINE',ERROR)
      CALL EXITS('MAKESPLINE')
      RETURN 1
      END



