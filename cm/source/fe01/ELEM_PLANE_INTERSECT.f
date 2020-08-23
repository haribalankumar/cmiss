      SUBROUTINE ELEM_PLANE_INTERSECT(IBT,IDO,INP,NBJ,ndep,ne,      
     '  NMAX_ITERATIONS,TOLERANCE,XE,XI,XNORM,XPOINT,ERROR,      
     '  SOLVED,*)

C#### Subroutine: ELEM_PLANE_INTERSECT
C###  Description:
C###    Calculates the intersection point in terms of xi coords of a
C###    plane and element given one xi coordinate as fixed and the
C###    other xi coordinate a dependant variable.      
C###    The initial value of xi coord being solved for is used as a
C###    seed value for the numerical method.  The Newton-Rhapson
C###    technique is used to solve for the dependant xi coordinate
C###    The plane is specified in point normal form
C**** Created by Peter Bier, April 2003
      
      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      
!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),      
     '  NBJ(NJM,NEM),ndep,ne,NMAX_ITERATIONS
      REAL*8 TOLERANCE,XE(NSM,NJM),XNORM(3),XPOINT(3)
      CHARACTER ERROR
      LOGICAL SOLVED
!     local variables
      INTEGER nb,nblah,nblahmax,ni,nj,nu
      REAL*8 D,F,FDERIV,FDX(3),FX(3),XI(3),XI_NEXT(3)
      LOGICAL FAILED
!     functions
      REAL*8 PXI
      
      CALL ENTERS('ELEM_PLANE_INTERSECT',*9999)
      NMAX_ITERATIONS = 20
      
C     The dependant variable is XI(ndep) where ndep = 1 or 2.
C     If ndep = 1 then we are solving for xi1, if ndep = 2 we are solving
C     for xi2.
C     The fixed XI coord will remain unchanged.

C     SUMMARY OF NEWTON-RHAPSON for this probelm:
C     We want to find a point (x,y,z) on the plane
C     ax + by + cz - d = 0
C     where x,y,z are on the element so
C     x=f(xi1,xi2), y=g(xi1,xi2), z=h(xi1,xi2)
C     where one of the xis is fixed
C     This means we want to solve F(xi1) = 0 where
C     F(xi1) = a*f(xi1,const) + b*g(xi1,const) + c*h(xi1,const) - d
C     F'(xi1) = a*f'(xi1,const) + b*g'(xi1,const) + c*h'(xi1,const)
C     
C     OR we want to solve F(xi2) = 0 where
C     F(xi2) = a*f(const,xi2) + b*g(const,xi2) + c*h(const,xi2) - d
C     F'(xi2) = a*f'(const,xi2) + b*g'(const,xi2) + c*h'(const,xi2)

C     To solve for dependant xi variable (xi1 or xi2) use Newton-Rhapson
C     xi_nplus1 = xi_n - F(xi_n) / F'(xi_n)
C     This scheme fails if a solution within the tolerance limit is not
C     found before the maximum number of iterations is reached or if the
C     value of F' is zero

C     D is the dot product of the normal vector and a point on the plane

C     make normal vector unit length
      CALL NORMALISE(3,XNORM,ERROR,*9999)
      
      D= 0.0d0
      DO nj = 1,3
        D = D + XNORM(nj) * XPOINT(nj)
      ENDDO

      IF (ndep .EQ. 1) THEN
        nu = 2 !derivative with respect to xi1
      ELSE
        nu = 4 !derivative with respect to xi2
      ENDIF

      nb = NBJ(1,ne)
      
      FAILED = .FALSE.
      SOLVED = .FALSE.
      ni = 1

      IF( .FALSE. ) THEN
C       Added to see form of objective funtion
        WRITE(OP_STRING,*) 'element ',ne, ' ndep ', ndep
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) 'plane normal ', XNORM
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) 'plane point ', XPOINT
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
        WRITE(*,*) '***********************'

        nblahmax = 10
        DO nblah = 1,nblahmax
C         calculate value of F for current xi
          XI(ndep) = 0.0d0 + nblah/(1.0d0*nblahmax)
          F = 0.0d0
          DO nj = 1,3
            FX(nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,1,XI,XE(1,nj))
            F = F + XNORM(nj) * FX(nj)
          ENDDO
          F = F - D
          
          FDERIV = 0.0d0
          DO nj = 1,3
            FDX(nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,nu,XI,XE(1,nj))
            FDERIV = FDERIV + XNORM(nj) * FDX(nj)
          ENDDO
          
          WRITE(OP_STRING,'(''X '',3F10.4,'' F '',F12.4,'' D '',F8.4)')        
     '      FX,F,D
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''DX '',3F10.4,'' DF '',F8.4)') FDX,FDERIV
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        
        ENDDO
        
        WRITE(*,*) '***********************'
        XI(ndep) = 0.5d0
      ENDIF
C*******************************************
      
      DO WHILE( ni .LE. NMAX_ITERATIONS .AND. .NOT.FAILED
     '  .AND. .NOT.SOLVED )

C       calculate value of FDERIV for current xi
        FDERIV = 0.0d0
        DO nj = 1,3
          FDX(nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '      INP(1,1,nb),nb,nu,XI,XE(1,nj))
          FDERIV = FDERIV + XNORM(nj) * FDX(nj)
        ENDDO

        IF( FDERIV .EQ. 0.0d0 ) THEN ! routine does not work if FDERIV = 0
          FAILED = .TRUE.
        ELSE
C         calculate value of F for current xi
          F = 0.0d0
          DO nj = 1,3
            FX(nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,1,XI,XE(1,nj))
            F = F + XNORM(nj) * FX(nj)
          ENDDO
          F = F - D

C         calculate next xi value
          XI_NEXT(ndep) = XI(ndep) - F / FDERIV

C          WRITE(OP_STRING,*) 'xi_next', XI_NEXT(ndep), F, FDERIV,          
C     '      XI(ndep)
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C         allow numbers outside 0 to 1 range during working
C         otherwise we miss some points whos initial xi_next values are
C         a little too big or a little too small 
          IF( DABS(XI_NEXT(ndep) - XI(ndep)) .LE. TOLERANCE ) THEN
            IF( (XI_NEXT(ndep) .GT. 1.0d0) .OR.
     '        (XI_NEXT(ndep) .LT. 0.0d0)) THEN
              FAILED = .TRUE. ! xi solution must be between zero and one
            ELSE
              SOLVED = .TRUE. ! xi solved within tolerance
C              WRITE(OP_STRING,*) 'SOLVED!', FX          
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              
            ENDIF
          ENDIF
          
C         store latest value in current xi
          XI(ndep) = XI_NEXT(ndep)
          
        ENDIF
        ni = ni + 1
      ENDDO

C     after exit, if FOUND is true xi has the solution value within
C     range 0 to 1
      IF(DOP .AND. .NOT.SOLVED) THEN
        WRITE(OP_STRING,'(''Performed '',I2,'' iterations of NR '
     '    // 'method'')') NMAX_ITERATIONS
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''Last xi value was '',F16.8)') XI(ndep)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        IF( FAILED ) THEN
          WRITE(OP_STRING,*) 'Failed to find an xi value ' //
     '      ' within the range 0 to 1'
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,*) 'NR did not converge'
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF
          

      CALL EXITS('ELEM_PLANE_INTERSECT')
      RETURN
 9999 CALL ERRORS('ELEM_PLANE_INTERSECT',ERROR)
      CALL EXITS('ELEM_PLANE_INTERSECT')
      RETURN 1
      END


