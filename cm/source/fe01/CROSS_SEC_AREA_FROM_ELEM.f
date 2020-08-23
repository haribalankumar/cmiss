      SUBROUTINE CROSS_SEC_AREA_FROM_ELEM(IBT,IDO,INP,NBJ,ne,
     '  AREA_FROM_ELEM,RADIUS,TOLERANCE,XE,XNORM,XPOINT,ERROR,*)

C#### Subroutine: CROSS_SEC_AREA_FROM_ELEM
C###  Description: Calculates the cross section area contribution
C###    for a particular element from the intersection of a plane in point
C###    normal from with the surface mesh element
C###    To calculate the contribution first locate the intersection points of
C###    the plane and the edges of the surface element.  Interior points are 
C###    then generated between the edge intersections and used with the centre
C###    point of the plane to calculate a triangular approximation to the area.
C###    The number of triangles used to approximate the area is doubled
C###    until the area contribution from the element converges to within
C###    the specified tolearnce level
C**** Created by Peter Bier, April 2003      
      
      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
          
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),      
     '  NBJ(NJM,NEM),ne
      REAL*8 AREA_FROM_ELEM,RADIUS,TOLERANCE,XE(NSM,NJM),XNORM(3),      
     '  XPOINT(3)
      CHARACTER ERROR*(*)
!     Parameter list

!     local variables
      INTEGER nb,ndep,nfixed,nin,ninaxis,nj,NINTOTAL,      
     '  NMAX_ITERATIONS,np,NPOINTS,NSEED_POINTS,nxi,NXI_AXIS_INT
      REAL*8 ELEM_XI_INTERSECT(4,3),END_XI(3),NEW_SECTOR_AREA,
     '  NEW_XI(3),OLD_SECTOR_AREA,P1(3),P2(3),START_XI(3),
     '  TRI_AREA,TRIANGLE_POINTS(257,3),
     '  TRIANGLE_XI(257,3),XI(3),XI_AXIS_INTERSECT(6,3),XI_INCREMENT
      LOGICAL CONVERGED,FOUND
!     functions
      REAL*8 PXI,TRIANGLE_AREA
      
      CALL ENTERS('CROSS_SEC_AREA_FROM_ELEM',*9999)

      NMAX_ITERATIONS = 10
      XI(3) = 0.0d0 !not used in calculations but must be set

      NSEED_POINTS = 1 ! only use one seed point

C     first locate any pairs of intersection points on the edges
C     Each pair is the start and end of a line.  
      nb = NBJ(1,ne)
      nin = 0 ! intersection point loop variable
      NINTOTAL = 0
      
      DO ndep = 1,2 ! dependent xi direction is first xi1 then xi2
        nfixed = MOD(ndep,2) + 1 ! nfixed is other xi coordinate
        DO nxi = 1,2
          XI(nfixed) = 1.0d0 * (nxi - 1) ! fixed xi is either 0 or 1
          XI(ndep) = 0.5d0 ! seed with a half
C         find all intersection points along a fixed xi axis
          CALL ELEM_PLANE_INTERSECT_XIAXIS(IBT,IDO,INP,NBJ,ndep,ne,
     '      NMAX_ITERATIONS,NSEED_POINTS,NXI_AXIS_INT,RADIUS,TOLERANCE,          
     '      XE,XI,XI_AXIS_INTERSECT,XNORM,XPOINT,ERROR,*9999)

C         add intersection point(s) to element list          
          DO ninaxis = 1,NXI_AXIS_INT ! 1 to total axis intercepts

            nin = nin + 1
            NINTOTAL = nin
            DO nj = 1,3                            
              ELEM_XI_INTERSECT(nin,nj) =
     '          XI_AXIS_INTERSECT(ninaxis,nj)
            ENDDO !nj = 1,3

          ENDDO !nin = 1,XI_AXIS_INTERCEPT
        ENDDO
      ENDDO

      IF( DOP ) THEN
        WRITE(OP_STRING,'(I3,'' pairs of intersection points'')') 
     '    NINTOTAL/2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C     ELEM_XI_INTERSECT should now contain all intersection points of
C     plane with the element.  Expect either 0 or 2 intersection points
C     for most cases.
            
      IF( MOD(NINTOTAL,2) .NE. 0 ) THEN ! odd number intersection points
        WRITE(OP_STRING,*) '>>Error: could not find even number of '//
     '    'intersection points for elem ', ne
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*) 'Returning zero for area contribution for '//
     '    'element'
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        NINTOTAL=0
      ELSEIF( NINTOTAL .GT. 2) THEN
        WRITE(OP_STRING,'(''>>WARNING: found'',I3,''intersections'')') 
     '    NINTOTAL
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''Code written but not tested '//
     '    'for >2points'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)        
      ENDIF
      
C     For each arc (specified in terms of end points) calculate the
C     surface area contribution by using triangles, double the number of
C     triangles used each time until the area approximation converges
     
      nin = 1
      AREA_FROM_ELEM = 0.0d0

      DO WHILE (nin .LT. NINTOTAL)
        
        DO nj = 1,3
          START_XI(nj) = ELEM_XI_INTERSECT(nin,nj)
          END_XI(nj) = ELEM_XI_INTERSECT(nin+1,nj)
        ENDDO
        
C       work out which xi direction to iterate along
C       along, choose the one with most variation
        IF( DABS(START_XI(1) - END_XI(1)) .LT.
     '    DABS(START_XI(2) - END_XI(2)) ) THEN
          ndep = 1 !iterate along xi1
          nfixed = 2
        ELSE
          ndep = 2 !iterate along xi2
          nfixed = 1
        ENDIF

        XI_INCREMENT = END_XI(nfixed) - START_XI(nfixed) ! may be +ve or -ve
        DO nj = 1,3
          TRIANGLE_XI(1,nj) = START_XI(nj)
          TRIANGLE_XI(2,nj) = END_XI(nj)
          TRIANGLE_POINTS(1,nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,1,START_XI,XE(1,nj)) ! evaluate start point coords
          P1(nj) = TRIANGLE_POINTS(1,nj)
          TRIANGLE_POINTS(2,nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '        INP(1,1,nb),nb,1,END_XI,XE(1,nj)) ! evaluate end point coords
          P2(nj) = TRIANGLE_POINTS(2,nj)
        ENDDO
        NPOINTS = 2
C       calculate first triangle area
        OLD_SECTOR_AREA = TRIANGLE_AREA(XPOINT,P1,P2)
C        WRITE(*,*) 'INITIAL AREA is ',OLD_SECTOR_AREA
        CONVERGED = .FALSE.
        
        DO WHILE( .NOT.CONVERGED .AND. NPOINTS.LE.129)
          
C         copy old points to odd array indexes to leave space for new
C         points
          DO np = NPOINTS,1,-1 ! start at end so no data overwritten
            DO nj = 1,3
              TRIANGLE_POINTS(2 * np -1,nj) = TRIANGLE_POINTS(np,nj)
              TRIANGLE_XI(2 * np - 1,nj) = TRIANGLE_XI(np,nj)
            ENDDO
          ENDDO
          
          XI_INCREMENT = XI_INCREMENT / 2.0d0
C         calaculate new point xi coords
          DO np = 1, NPOINTS - 1
            NEW_XI(nfixed) = TRIANGLE_XI(2 * np-1,nfixed)
     '        +  XI_INCREMENT !next point is last one plus xi increment
            NEW_XI(ndep) = TRIANGLE_XI(2 * np-1,ndep) !good initial guess,
C           initial guess is the other xi direction 
            NEW_XI(3) = 0.0d0
            
C         Calculate the xi coords of the xi dependant direction
            CALL ELEM_PLANE_INTERSECT(IBT,IDO,INP,NBJ,ndep,ne,
     '        NMAX_ITERATIONS,TOLERANCE,XE,NEW_XI,XNORM,XPOINT,ERROR,            
     '        FOUND,*9999)

C           when very close to an element edge it is possible for
C           intersections for the mid points to be outside the element
C           if this is the case then use points on the element edge instead
            IF(.NOT.FOUND) THEN
              WRITE(OP_STRING,'(1X,A,3F10.5)') 'Warning: Could not ' //
     '          'find plane intersection within elem for XI ', NEW_XI
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

              IF( NEW_XI(ndep) .GT. 1.0d0) THEN
                NEW_XI(ndep) = 1.0d0
              ELSEIF( NEW_XI(ndep) .LT. 0.0d0) THEN
                NEW_XI(ndep) = 0.0d0
              ENDIF
              WRITE(OP_STRING,'(1X,A,3F10.5)') 'Using XI ', NEW_XI
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                            
            ENDIF

C           store coords of new point
            TRIANGLE_XI(2 * np,nfixed) = NEW_XI(nfixed)
            TRIANGLE_XI(2 * np,ndep) = NEW_XI(ndep)
            TRIANGLE_XI(2 * np,3) = NEW_XI(3)

            DO nj = 1,3
              TRIANGLE_POINTS(2*np,nj) = PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '          INP(1,1,nb),nb,1,NEW_XI,XE(1,nj)) ! evaluate postion at XI axis intersect
            ENDDO

          ENDDO
          
          NPOINTS = 2 * NPOINTS - 1 ! store new number of points in arrays

C         calculate new sector area
          np = 1
          NEW_SECTOR_AREA = 0.0d0
          DO WHILE(np .LT. NPOINTS)
            DO nj = 1,3
              P1(nj) = TRIANGLE_POINTS(np,nj)
              P2(nj) = TRIANGLE_POINTS(np+1,nj)
            ENDDO        
            TRI_AREA = TRIANGLE_AREA(XPOINT,P1,P2)
            NEW_SECTOR_AREA = NEW_SECTOR_AREA + TRI_AREA
            np = np + 1 
          ENDDO      

C          WRITE(*,*) 'NEW SECTOR AREA is',NEW_SECTOR_AREA
          IF( DABS(NEW_SECTOR_AREA - OLD_SECTOR_AREA) .LT. TOLERANCE )          
     '      THEN
            CONVERGED = .TRUE.
            AREA_FROM_ELEM = NEW_SECTOR_AREA
          ELSE
            OLD_SECTOR_AREA = NEW_SECTOR_AREA
          ENDIF
        ENDDO  ! .NOT.CONVERGED
        nin = nin + 2 ! advance to next pair
      ENDDO
      
      
      CALL EXITS('CROSS_SEC_AREA_FROM_ELEM')
      RETURN
 9999 CALL ERRORS('CROSS_SEC_AREA_FROM_ELEM',ERROR)
      CALL EXITS('CROSS_SEC_AREA_FROM_ELEM')
      RETURN
      END


