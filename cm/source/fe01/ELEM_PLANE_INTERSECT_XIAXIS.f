      SUBROUTINE ELEM_PLANE_INTERSECT_XIAXIS(IBT,IDO,INP,NBJ,ndep,
     '  ne,NMAX_ITERATIONS,NSEED_POINTS,NXI_AXIS_INT,RADIUS,TOLERANCE,      
     '  XE,XI,XI_AXIS_INTERSECT,XNORM,XPOINT,ERROR,*)

C#### Subroutine: ELEM_PLANE_INTERSECT_XIAXIS
C###  Description:
C###    Calculates intersection points in terms of xi coords of a
C###    plane with and a line described by a fixed xi axis value
C###    Any point not within RADIUS of the point on the plane is
C###    discarded.  The number of seed points to use is specified.
C###    In a small number of special cases a different seed point
C###    may lead to a second intersection along the axis if the
C###    element is curved enough.
C**** Created by Peter Bier, April 2003      
      
      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'

!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),ndep,ne,NMAX_ITERATIONS,      
     '  NSEED_POINTS,NXI_AXIS_INT
      REAL*8 RADIUS,TOLERANCE,XE(NSM,NJM),XI(3),
     '  XI_AXIS_INTERSECT(6,3),XNORM(3),XPOINT(3)
      CHARACTER ERROR
!     local variables
      INTEGER nb,nin,nj,np,nseed
      REAL*8 D,DP,INTERCEPT_POINT(6,3),LIST_POINT(3),XI_SEED,      
     '  X_INTERCEPT(3),SEED_INC
      LOGICAL INLIST,SOLVED
!     functions
      REAL*8 DISTANCE,PXI
      
      CALL ENTERS('ELEM_PLANE_INTERSECT_XIAXIS',*9999)


C     The dependant variable is XI(ndep) where ndep = 1 or 2.
C     If ndep = 1 then we are solving for xi1, if ndep = 2 we are solving
C     for xi2.
C     The fixed XI coord will remain unchanged.  A number of different
C     seed values are tried for the dependant variable

      nin = 0
      IF( NSEED_POINTS .LE. 1 ) THEN
        NSEED_POINTS = 1
        XI_SEED = 0.5d0 ! seed point (only one used)
        SEED_INC = 0.0d0
      ELSE
        XI_SEED = 0.0d0 ! seed point to start with
        SEED_INC = 1.0d0 /(NSEED_POINTS - 1)
      ENDIF
      
      XI(ndep) = XI_SEED ! initial seed point

      DO nseed = 1,NSEED_POINTS

        CALL ELEM_PLANE_INTERSECT(IBT,IDO,INP,NBJ,ndep,ne,      
     '    NMAX_ITERATIONS,TOLERANCE,XE,XI,XNORM,XPOINT,ERROR,      
     '    SOLVED,*9999)
        IF( SOLVED ) THEN
C         find distance from intercept point to plane point
          nb = NBJ(1,ne)
          DO nj = 1,3
            X_INTERCEPT(nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,1,XI,XE(1,nj))
          ENDDO ! nj = 1,3          

          D = DISTANCE(3,XPOINT,X_INTERCEPT)
C         only consider adding point to the list if it is within the
C         RADIUS of interest                    
          IF( (D .LT. RADIUS)) THEN ! D .LT. RADIUS

            IF(nin .EQ. 0) THEN !first point so just add to list
              nin = nin + 1
              DO nj = 1,3
                INTERCEPT_POINT(nin,nj) = X_INTERCEPT(nj)
                XI_AXIS_INTERSECT(nin,nj) = XI(nj)
              ENDDO
            ELSE ! check to see if point already in list before adding
              INLIST = .FALSE.
              DO np = 1,nin
                DO nj = 1,3
                  LIST_POINT(nj) = INTERCEPT_POINT(np,nj)
                ENDDO            
                DP = DISTANCE(3,X_INTERCEPT,LIST_POINT)
                IF( DP .LT. TOLERANCE ) THEN
                  INLIST = .TRUE.
                ENDIF
              ENDDO ! np = 1,nin 

              IF( .NOT.INLIST ) THEN ! add point as not in list
                nin = nin + 1
                DO nj = 1,3
                  INTERCEPT_POINT(nin,nj) = X_INTERCEPT(nj)
                  XI_AXIS_INTERSECT(nin,nj) = XI(nj)
                ENDDO
              ENDIF ! .NOT.INLIST              
            ENDIF ! nin .EQ. 0
          ENDIF ! D .LT. RADIUS
        ENDIF ! solved

        XI(ndep) = XI(ndep) + SEED_INC ! increment to next seed point to check

      ENDDO ! nseed = 1,NSEED_POINTS
      NXI_AXIS_INT = nin
      
      CALL EXITS('ELEM_PLANE_INTERSECT_XIAXIS')
      RETURN
 9999 CALL ERRORS('ELEM_PLANE_INTERSECT_XIAXIS',ERROR)
      CALL EXITS('ELEM_PLANE_INTERSECT_XIAXIS')
      RETURN 1

      END
      

          
