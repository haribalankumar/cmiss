      SUBROUTINE NORM_LATTICE(NDIM,nq,NWQ,DXDXIQ,DXDXIQ2,XNLOCAL,
     '  ERROR,*)

C#### Subroutine: NORM_LATTICE
C###  Description:
C###    NORM_LATTICE finds the normal vector XNLOCAL(nj) in
C###    the cartesian reference frame. This routine is based
C###    upon NORM31 but uses a modified NWQ array to find
C###    what direction the boundary lies in instead of NXQ.      
C###  See-Also: NORM31
C###  See-Also: NWQ      

      IMPLICIT NONE
      
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
      
!     Parameter List
      INTEGER NDIM,nq,NWQ(8)
      REAL*8 DXDXIQ(3,3),DXDXIQ2(3,3),XNLOCAL(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER METHOD,nj
      REAL*8 VECTOR1(3),VECTOR2(3),VECTOR3(3),VLENGTH

      CALL ENTERS('NORM_LATTICE',*9999)

      CALL ASSERT(USE_LAT.EQ.1,'>>NORM_LATTICE only for lattice method',
     '  ERROR,*9999)

      DO nj=1,3
        XNLOCAL(nj)=0.0d0
        VECTOR1(nj)=0.0d0
        VECTOR2(nj)=0.0d0
        VECTOR3(nj)=0.0d0
      ENDDO

      IF(NDIM.EQ.1) THEN
        IF(NWQ(1).EQ.1) THEN
          XNLOCAL(1)=-1.0d0
        ELSEIF(NWQ(1).EQ.2) THEN
          XNLOCAL(1)=1.0d0
        ELSEIF(NWQ(1).EQ.0) THEN
          ERROR='Grid point is not a boundary point'
          GOTO 9999
        ENDIF
      ELSEIF(NDIM.EQ.2) THEN
        METHOD=3
        IF(METHOD.EQ.3) THEN
          IF(NWQ(1).EQ.1.OR.NWQ(1).EQ.7.OR.NWQ(1).EQ.8.OR.
     '      NWQ(1).EQ.11.OR.NWQ(1).EQ.12.OR.NWQ(1).EQ.19.OR.
     '      NWQ(1).EQ.20.OR.NWQ(1).EQ.21.OR.NWQ(1).EQ.22) THEN
            XNLOCAL(1)=-DXDXIQ2(2,2)
            XNLOCAL(2)=DXDXIQ2(1,2)
          ELSEIF(NWQ(1).EQ.2.OR.NWQ(1).EQ.9.OR.NWQ(1).EQ.10.OR.
     '      NWQ(1).EQ.13.OR.NWQ(1).EQ.14.OR.NWQ(1).EQ.23.OR.
     '      NWQ(1).EQ.24.OR.NWQ(1).EQ.25.OR.NWQ(1).EQ.26) THEN
            XNLOCAL(1)=DXDXIQ2(2,2)
            XNLOCAL(2)=-DXDXIQ2(1,2)
          ELSEIF(NWQ(1).EQ.3.OR.NWQ(1).EQ.7.OR.NWQ(1).EQ.9.OR.
     '      NWQ(1).EQ.15.OR.NWQ(1).EQ.16.OR.NWQ(1).EQ.19.OR.
     '      NWQ(1).EQ.20.OR.NWQ(1).EQ.23.OR.NWQ(1).EQ.24) THEN
            XNLOCAL(1)=-DXDXIQ2(2,1)
            XNLOCAL(2)=DXDXIQ2(1,1)
          ELSEIF(NWQ(1).EQ.4.OR.NWQ(1).EQ.8.OR.NWQ(1).EQ.10.OR.
     '      NWQ(1).EQ.17.OR.NWQ(1).EQ.18.OR.NWQ(1).EQ.21.OR.
     '      NWQ(1).EQ.22.OR.NWQ(1).EQ.25.OR.NWQ(1).EQ.26) THEN
            XNLOCAL(1)=DXDXIQ2(2,1)
            XNLOCAL(2)=-DXDXIQ2(1,1)
          ELSE
            !Not quite this
            XNLOCAL(1)=-1.0d0
            XNLOCAL(2)=0.0d0
          ENDIF
          VLENGTH=DSQRT((XNLOCAL(1)**2)+(XNLOCAL(2)**2))
          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(1)=XNLOCAL(1)/VLENGTH
          IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(2)=XNLOCAL(2)/VLENGTH
        ENDIF
        IF(XNLOCAL(1).EQ.0.0d0.AND.XNLOCAL(2).EQ.0.0d0) THEN
          WRITE(*,*) "Zero normal at ",nq
          XNLOCAL(2)=1.0d0
        ENDIF        
      ELSEIF(NDIM.EQ.3) THEN
CC       Build the jacobian
C        DO ni=1,3
C          DO nj=1,3
C            JAC_TMP(nj,ni)=DXDXIQ(nj,ni)
C          ENDDO
C        ENDDO

        IF(NWQ(1).EQ.1) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)        
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
        ELSEIF(NWQ(1).EQ.7) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF        
        ELSEIF(NWQ(1).EQ.8) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF        
        ELSEIF(NWQ(1).EQ.11) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF        
        ELSEIF(NWQ(1).EQ.12) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF        
        ELSEIF(NWQ(1).EQ.19) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR2)
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR3)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF
          IF(VECTOR3(1).EQ.0.0d0.AND.VECTOR3(2).EQ.0.0d0.AND.
     '      VECTOR3(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR3"
          ENDIF        
        ELSEIF(NWQ(1).EQ.20) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR2)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR3)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF
          IF(VECTOR3(1).EQ.0.0d0.AND.VECTOR3(2).EQ.0.0d0.AND.
     '      VECTOR3(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR3"
          ENDIF        
        ELSEIF(NWQ(1).EQ.21) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR2)
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR3)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF
          IF(VECTOR3(1).EQ.0.0d0.AND.VECTOR3(2).EQ.0.0d0.AND.
     '      VECTOR3(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR3"
          ENDIF        
        ELSEIF(NWQ(1).EQ.22) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,2),VECTOR1)
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR2)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR3)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF
          IF(VECTOR3(1).EQ.0.0d0.AND.VECTOR3(2).EQ.0.0d0.AND.
     '      VECTOR3(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR3"
          ENDIF        
        ELSEIF(NWQ(1).EQ.2) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF        
        ELSEIF(NWQ(1).EQ.9) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF        
        ELSEIF(NWQ(1).EQ.10) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF          
        ELSEIF(NWQ(1).EQ.13) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF          
        ELSEIF(NWQ(1).EQ.14) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF          
        ELSEIF(NWQ(1).EQ.23) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR2)
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR3)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF
          IF(VECTOR3(1).EQ.0.0d0.AND.VECTOR3(2).EQ.0.0d0.AND.
     '      VECTOR3(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR3"
          ENDIF        
        ELSEIF(NWQ(1).EQ.24) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR2)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR3)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF
          IF(VECTOR3(1).EQ.0.0d0.AND.VECTOR3(2).EQ.0.0d0.AND.
     '      VECTOR3(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR3"
          ENDIF        
        ELSEIF(NWQ(1).EQ.25) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR2)
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR3)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF
          IF(VECTOR3(1).EQ.0.0d0.AND.VECTOR3(2).EQ.0.0d0.AND.
     '      VECTOR3(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR3"
          ENDIF        
        ELSEIF(NWQ(1).EQ.26) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR2)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR3)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF
          IF(VECTOR3(1).EQ.0.0d0.AND.VECTOR3(2).EQ.0.0d0.AND.
     '      VECTOR3(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR3"
          ENDIF        
        ELSEIF(NWQ(1).EQ.3) THEN
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR1)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF        
        ELSEIF(NWQ(1).EQ.15) THEN
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF        
        ELSEIF(NWQ(1).EQ.16) THEN
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,3),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF         
        ELSEIF(NWQ(1).EQ.4) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR1)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF        
        ELSEIF(NWQ(1).EQ.17) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR1)
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF         
        ELSEIF(NWQ(1).EQ.18) THEN
          CALL CROSS(DXDXIQ(1,3),DXDXIQ(1,1),VECTOR1)
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR2)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF
          IF(VECTOR2(1).EQ.0.0d0.AND.VECTOR2(2).EQ.0.0d0.AND.
     '      VECTOR2(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR2"
          ENDIF         
        ELSEIF(NWQ(1).EQ.5) THEN
          CALL CROSS(DXDXIQ(1,2),DXDXIQ(1,1),VECTOR1)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF        
        ELSE
          CALL CROSS(DXDXIQ(1,1),DXDXIQ(1,2),VECTOR1)
          IF(VECTOR1(1).EQ.0.0d0.AND.VECTOR1(2).EQ.0.0d0.AND.
     '      VECTOR1(3).EQ.0.0d0) THEN
            WRITE(*,*) "DEGENERATE VECTOR1"
          ENDIF        
        ENDIF

        XNLOCAL(1)=VECTOR1(1)+VECTOR2(1)+VECTOR3(1)
        XNLOCAL(2)=VECTOR1(2)+VECTOR2(2)+VECTOR3(2)
        XNLOCAL(3)=VECTOR1(3)+VECTOR2(3)+VECTOR3(3)       
      
        VLENGTH=DSQRT((XNLOCAL(1)**2)+(XNLOCAL(2)**2)+(XNLOCAL(3)**2))
        IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(1)=XNLOCAL(1)/VLENGTH
        IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(2)=XNLOCAL(2)/VLENGTH
        IF(VLENGTH.GE.ZERO_TOL) XNLOCAL(3)=XNLOCAL(3)/VLENGTH
        
        IF(XNLOCAL(1).EQ.0.0d0.AND.XNLOCAL(2).EQ.0.0d0.AND.
     '    XNLOCAL(3).EQ.0.0d0) THEN
          WRITE(*,*) "Zero normal at ",nq
          XNLOCAL(2)=1.0d0
        ENDIF
      ENDIF  

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Normal at '',I6, ''is '',3F8.4)') nq,
     '    XNLOCAL(1),XNLOCAL(2),XNLOCAL(3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Vector length '',F12.8)') VLENGTH
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DOP=.FALSE.

      CALL EXITS('NORM_LATTICE')
      RETURN
 9999 CALL ERRORS('NORM_LATTICE',ERROR)
      CALL EXITS('NORM_LATTICE')
      RETURN 1      
      END


