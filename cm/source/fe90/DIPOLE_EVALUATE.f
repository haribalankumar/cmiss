      SUBROUTINE DIPOLE_EVALUATE(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IVAL,
     '  NDIPOLES,np,NP_INTERFACE,nr,nx,CE,DIPOLE_CEN,DIPOLE_DIR,XP,W,
     '  X,Y,Z,ERROR,*)

C#### Subroutine: DIPOLE_EVALUATE
C###  Description:
C###    DIPOLE_EVALUATE evaluates the solution W at the given point
C###    x,y,z of a dipole inside concentric circles/spheres.

C**** If IVAL = 1 W contains the solution
C****         = 2 W contains the s1 derivative
C****         = 3 W contains the s2 derivative
C****         = 4 W contains the normal derivative
C****         = 5 W contains s1 derivative of the normal derivative
C****         = 6 W contains s2 derivative of the normal derivative
C****         = 7 W contains s1&s2 (cross) derivative
C****         = 8 W contains s1&s2 derivative of the normal derivative
C**** The actual solution used is based on ANAL_CHOICE (in anal00.cmn).
C**** For 2D:
C****   If ANAL_CHOICE(1) = 5 then soln is for a single circle
C****                     = 6 then soln is for multiple circles
C**** For 3D:
C****   If ANAL_CHOICE(1) = 5 then soln is for a dipole on the x axis
C****                     = 6                "                 y axis
C****                     = 7                "                 z axis
C****                     = 8 then the soln is the concentric spheres
C****                         soln where the source is a single Z dipole
C****                         with no flux boundary conditions on the
C****                         outer sphere.
C****                     = 14 then the dipole is eccentric, at any
C****                         orientation
C**** NOTE : Assumes that the geometry is spherical for the calculation
C****        of the derivatives (particularly the normal derivative).

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'mesh03.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter list
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),IVAL,NDIPOLES(NRM),np,
     '  NP_INTERFACE(0:NPM,0:3),nr,nx
      REAL*8 CE(NMM),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),XP(NKM,NVM,NJM),W,X,Y,Z
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER a_coeff,b_coeff,c_coeff,circle,index,ndip,nj,NSERIES,
     '  nv,MSERIES,sphere
      REAL*8 COND,DATAN_MOD,
     '  first_compo,first_compo_s1,first_compo_s2,second_compo,
     '  GRAD_POTENTIAL(3),M_PLGNDR_SIN,PHI,PLGNDR,
     '  POSITION(3),POTENTIAL,R,s1_coeff,s2_coeff,s1_dir(3),
     '  s2_dir(3),
     '  THETA,U,V,W_TEMP,XX(3),ZZ(3)

      CALL ENTERS('DIPOLE_EVALUATE',*9999)

      nv=1 !Temporary MPN 12-Nov-94
      W=0.0d0 !Eccentric multiple spheres is incremental
      IF(ANAL_CHOICE(nr).NE.14) THEN
C CPB 26/1/96 adding 2D dipole solutions
        IF(NJT.EQ.2) THEN
          R=DSQRT(X*X+Y*Y)
          THETA=DATAN_MOD(X,Y)
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(DIPOLE_EVALUATE_1)
            WRITE(OP_STRING,'('' Nodal coordinates:'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' X='',D12.4,'', Y='',D12.4)') X,Y
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' R='',D12.4,'', THETA='',D12.4)')
     '        R,THETA
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(DIPOLE_EVALUATE_1)
          ENDIF
          IF(ANAL_CHOICE(nr).EQ.5) THEN !Single circle
            IF(IVAL.EQ.1) THEN
              W=-(ANAL_A*DCOS(THETA)+ANAL_B*DSIN(THETA))/
     '          (PI*SIGMA(1)*R)
            ELSE IF(IVAL.EQ.2) THEN
              W=-(ANAL_A*DSIN(THETA)-ANAL_B*DCOS(THETA))/
     '          (PI*SIGMA(1)*R**2)
            ELSE IF(IVAL.EQ.4) THEN
              W=0.0d0 !no flux
            ELSE IF(IVAL.EQ.5) THEN
              W=0.0d0 !no flux
            ELSE
              ERROR='>>Invalid IVAL number'
              GOTO 9999
            ENDIF
          ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !Mulitple circles
            circle=nr !the circle number/region the node is in.
            IF(circle.EQ.1) THEN !in the first circle
              IF(IVAL.EQ.1) THEN !value
                V=-(ANAL_A*DCOS(THETA)+ANAL_B*DSIN(THETA))/(2.0d0*
     '            PI*SIGMA(1)*R)
                U=DIPOLE_COEFF(1,1)*R*(ANAL_A*DCOS(THETA)+
     '            ANAL_B*DSIN(THETA))
              ELSE IF(IVAL.EQ.2) THEN !s derivative
                V=(-ANAL_A*DSIN(THETA)+ANAL_B*DCOS(THETA))/(2.0d0*
     '            PI*SIGMA(1)*R**2)
                U=DIPOLE_COEFF(1,1)*(ANAL_A*DSIN(THETA)-
     '            ANAL_B*DCOS(THETA))
              ELSE IF(IVAL.EQ.4) THEN !normal derivative
                V=(ANAL_A*DCOS(THETA)+ANAL_B*DSIN(THETA))/(2.0d0*
     '            PI*SIGMA(1)*R**2)
                U=DIPOLE_COEFF(1,1)*(ANAL_A*DCOS(THETA)+
     '            ANAL_B*DSIN(THETA))
              ELSE IF(IVAL.EQ.5) THEN !s deriv of normal derivative
                V=(ANAL_A*DSIN(THETA)-ANAL_B*DCOS(THETA))/(2.0d0*
     '            PI*SIGMA(1)*R**3)
                U=DIPOLE_COEFF(1,1)/R*(ANAL_A*DSIN(THETA)-
     '            ANAL_B*DCOS(THETA))
              ELSE
                ERROR='>>Invalid IVAL number'
                GOTO 9999
              ENDIF
            ELSE !not in the first circle
              b_coeff=2*(circle-1)
              c_coeff=2*(circle-1)+1
              IF(IVAL.EQ.1) THEN !value
                V=-(ANAL_A*DCOS(THETA)+ANAL_B*DSIN(THETA))/(2.0d0*
     '            PI*SIGMA(1)*R)
                U=(DIPOLE_COEFF(b_coeff,1)*R+DIPOLE_COEFF(c_coeff,1)/
     '            R)*(ANAL_A*DCOS(THETA)+ANAL_B*DSIN(THETA))
              ELSE IF(IVAL.EQ.2) THEN !s derivative
                V=(-ANAL_A*DSIN(THETA)+ANAL_B*DCOS(THETA))/(2.0d0*
     '            PI*SIGMA(1)*R**2)
                U=(DIPOLE_COEFF(b_coeff,1)+DIPOLE_COEFF(c_coeff,1)/
     '            (R**2))*(ANAL_A*DSIN(THETA)-ANAL_B*DCOS(THETA))
              ELSE IF(IVAL.EQ.4) THEN !normal derivative
                V=(ANAL_A*DCOS(THETA)+ANAL_B*DSIN(THETA))/(2.0d0*
     '            PI*SIGMA(1)*R**2)
                U=(DIPOLE_COEFF(b_coeff,1)-DIPOLE_COEFF(c_coeff,1)/
     '            (R**2))*(ANAL_A*DCOS(THETA)+ANAL_B*DSIN(THETA))
              ELSE IF(IVAL.EQ.5) THEN !s deriv of normal derivative
                V=(ANAL_A*DSIN(THETA)-ANAL_B*DCOS(THETA))/(2.0d0*
     '            PI*SIGMA(1)*R**3)
                U=(DIPOLE_COEFF(b_coeff,1)/R-DIPOLE_COEFF(c_coeff,1)/
     '            (R**3))*(ANAL_A*DSIN(THETA)-ANAL_B*DCOS(THETA))
              ELSE
                ERROR='>>Invalid IVAL number'
                GOTO 9999
              ENDIF
            ENDIF
            W=V-U
            IF(IVAL.GE.4) W=SIGMA(circle)*W !adjust normal derivatives
C                                            for conductivity
          ENDIF
        ELSE IF(NJT.EQ.3) THEN
          IF(ANAL_CHOICE(nr).EQ.12) THEN !Single sphere
            ZZ(1)=X
            ZZ(2)=Y
            ZZ(3)=Z
            CALL ZX(3,ZZ,XX)
            R=XX(1)
            THETA=XX(2)
            PHI=XX(3)
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(DIPOLE_EVALUATE_2)
              WRITE(OP_STRING,'('' Nodal coordinates:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' X='',D12.4,'', Y='','
     '          //'D12.4,'', Z='',D12.4)') X,Y,Z
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' R='',D12.4,'', THETA='',D12.4,'
     '          //''', PHI='',D12.4)') R,THETA,PHI
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(DIPOLE_EVALUATE_2)
            ENDIF
            IF(IVAL.EQ.1) THEN !value
              W=3.0d0/(4.0d0*PI*SIGMA(1)*R**2)*(ANAL_A*
     '          DCOS(THETA)*DCOS(PHI)+ANAL_B*DSIN(THETA)*
     '          DCOS(PHI)+ANAL_C*DSIN(PHI))
            ELSE IF(IVAL.EQ.2) THEN !s1 derivative
              W=3.0d0/(4.0d0*PI*SIGMA(1)*R**3)*(-ANAL_A*
     '          DSIN(THETA)+ANAL_B*DCOS(THETA))
            ELSE IF(IVAL.EQ.3) THEN !s2 derivative
              W=3.0d0/(4.0d0*PI*SIGMA(1)*R**3)*(-ANAL_A*
     '          DCOS(THETA)*DSIN(PHI)-ANAL_B*DSIN(THETA)*
     '          DSIN(PHI)+ANAL_C*DCOS(PHI))
            ELSE IF(IVAL.EQ.4) THEN !normal derivative
              W=0.0d0 !no flux
            ELSE IF(IVAL.EQ.5) THEN !s1 deriv of normal derivative
              W=0.0d0 !no flux
            ELSE IF(IVAL.EQ.6) THEN !s2 deriv of normal derivative
              W=0.0d0 !no flux
            ELSE IF(IVAL.EQ.7) THEN !cross deriv of value
              W=0.0d0
            ELSE IF(IVAL.EQ.8) THEN !cross deriv of normal deriv
              W=0.0d0 !no flux
            ELSE
              ERROR='>>Invalid IVAL number'
              GOTO 9999
            ENDIF
          ELSE IF(ANAL_CHOICE(nr).EQ.13) THEN !Multiple spheres
            ZZ(1)=X
            ZZ(2)=Y
            ZZ(3)=Z
            CALL ZX(3,ZZ,XX)
            R=XX(1)
            THETA=XX(2)
            PHI=XX(3)
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(DIPOLE_EVALUATE_3)
              WRITE(OP_STRING,'('' Nodal coordinates:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' X='',D12.4,'', Y='','
     '          //'D12.4,'', Z='',D12.4)') X,Y,Z
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' R='',D12.4,'', THETA='',D12.4,'
     '          //''', PHI='',D12.4)') R,THETA,PHI
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(DIPOLE_EVALUATE_3)
            ENDIF
            sphere=nr !the sphere number/region the node is in.
            IF(sphere.EQ.1) THEN !in the first sphere
              IF(IVAL.EQ.1) THEN !value
                V=1.0d0/(4.0d0*PI*SIGMA(1)*R**2)*(ANAL_A*
     '            DCOS(THETA)*DCOS(PHI)+ANAL_B*DSIN(THETA)*
     '            DCOS(PHI)+ANAL_C*DSIN(PHI))
                U=DIPOLE_COEFF(1,1)*R*(ANAL_A*DCOS(THETA)*DCOS(PHI)+
     '            ANAL_B*DSIN(THETA)*DCOS(PHI)+ANAL_C*DSIN(PHI))
              ELSE IF(IVAL.EQ.2) THEN !s1 derivative
                V=1.0d0/(4.0d0*PI*SIGMA(1)*R**3)*(-ANAL_A*
     '            DSIN(THETA)+ANAL_B*DCOS(THETA))
                U=DIPOLE_COEFF(1,1)*(-ANAL_A*DSIN(THETA)+
     '            ANAL_B*DCOS(THETA))
              ELSE IF(IVAL.EQ.3) THEN !s2 derivative
                V=1.0d0/(4.0d0*PI*SIGMA(1)*R**3)*(-ANAL_A*
     '            DCOS(THETA)*DSIN(PHI)-ANAL_B*DSIN(THETA)*
     '            DSIN(PHI)+ANAL_C*DCOS(PHI))
                U=DIPOLE_COEFF(1,1)*(-ANAL_A*DCOS(THETA)*
     '            DSIN(PHI)-ANAL_B*DSIN(THETA)*DSIN(PHI)+ANAL_C*
     '            DCOS(PHI))
              ELSE IF(IVAL.EQ.4) THEN !normal derivative
                V=-1.0d0/(2.0d0*PI*SIGMA(1)*R**3)*(ANAL_A*
     '            DCOS(THETA)*DCOS(PHI)+ANAL_B*DSIN(THETA)*
     '            DCOS(PHI)+ANAL_C*DSIN(PHI))
                U=DIPOLE_COEFF(1,1)*(ANAL_A*DCOS(THETA)*
     '            DCOS(PHI)+ANAL_B*DSIN(THETA)*DCOS(PHI)+
     '            ANAL_C*DSIN(PHI))
              ELSE IF(IVAL.EQ.5) THEN !s1 deriv of normal deriv
                V=1.0d0/(2.0d0*PI*SIGMA(1)*R**4)*(ANAL_A*
     '            DSIN(THETA)-ANAL_B*DCOS(THETA))
                U=DIPOLE_COEFF(1,1)/R*(-ANAL_A*DSIN(THETA)
     '            +ANAL_B*DCOS(THETA))
              ELSE IF(IVAL.EQ.6) THEN !s2 deriv of normal deriv
                V=1.0d0/(2.0d0*PI*SIGMA(1)*R**4)*(ANAL_A*
     '            DCOS(THETA)*DSIN(PHI)+ANAL_B*DSIN(THETA)*
     '            DSIN(PHI)-ANAL_C*DCOS(PHI))
                U=DIPOLE_COEFF(1,1)/R*
     '            (-ANAL_A*DCOS(THETA)*DSIN(PHI)-
     '            ANAL_B*DSIN(THETA)*DSIN(PHI)+ANAL_C*DCOS(PHI))
              ELSE IF(IVAL.EQ.7) THEN ! s1&s2 derivative
                V=0.0d0
                U=0.0d0
              ELSE IF(IVAL.EQ.8) THEN ! s1&s2 deriv of normal deriv
                V=0.0d0
                U=0.0d0
              ELSE
                ERROR='>>Invalid IVAL number'
                GOTO 9999
              ENDIF
            ELSE !not in the first sphere
              b_coeff=2*(sphere-1)
              c_coeff=2*(sphere-1)+1
              IF(IVAL.EQ.1) THEN !value
                V=1.0d0/(4.0d0*PI*SIGMA(1)*R**2)*(ANAL_A*DCOS(THETA)*
     '            DCOS(PHI)+ANAL_B*DSIN(THETA)*DCOS(PHI)+ANAL_C*
     '            DSIN(PHI))
                U=(DIPOLE_COEFF(b_coeff,1)*R+DIPOLE_COEFF(c_coeff,1)/
     '            R**2)*(ANAL_A*DCOS(THETA)*DCOS(PHI)+
     '            ANAL_B*DSIN(THETA)*DCOS(PHI)+ANAL_C*DSIN(PHI))
              ELSE IF(IVAL.EQ.2) THEN !s1 derivative
                V=1.0d0/(4.0d0*PI*SIGMA(1)*R**3)*(-ANAL_A*
     '            DSIN(THETA)+ANAL_B*DCOS(THETA))
                U=(DIPOLE_COEFF(b_coeff,1)+DIPOLE_COEFF(c_coeff,1)/
     '            R**3)*(-ANAL_A*DSIN(THETA)
     '            +ANAL_B*DCOS(THETA))
              ELSE IF(IVAL.EQ.3) THEN !s2 derivative
                V=1.0d0/(4.0d0*PI*SIGMA(1)*R**3)*(-ANAL_A*
     '            DCOS(THETA)*DSIN(PHI)-ANAL_B*DSIN(THETA)*
     '            DSIN(PHI)+ANAL_C*DCOS(PHI))
                U=(DIPOLE_COEFF(b_coeff,1)+DIPOLE_COEFF(c_coeff,1)/
     '            R**3)*(-ANAL_A*DCOS(THETA)*
     '            DSIN(PHI)-ANAL_B*DSIN(THETA)*DSIN(PHI)+ANAL_C*
     '            DCOS(PHI))
              ELSE IF(IVAL.EQ.4) THEN !normal derivative
                V=-1.0d0/(2.0d0*PI*SIGMA(1)*R**3)*(ANAL_A*
     '            DCOS(THETA)*DCOS(PHI)+ANAL_B*DSIN(THETA)*
     '            DCOS(PHI)+ANAL_C*DSIN(PHI))
                U=(DIPOLE_COEFF(b_coeff,1)-
     '            2.0d0*DIPOLE_COEFF(c_coeff,1)/R**3)*
     '            (ANAL_A*DCOS(THETA)*DCOS(PHI)+ANAL_B*DSIN(THETA)*
     '            DCOS(PHI)+ANAL_C*DSIN(PHI))
              ELSE IF(IVAL.EQ.5) THEN !s1 deriv of normal deriv
                V=1.0d0/(2.0d0*PI*SIGMA(1)*R**4)*(ANAL_A*
     '            DSIN(THETA)-ANAL_B*DCOS(THETA))
                U=(DIPOLE_COEFF(b_coeff,1)/R-
     '            2.0d0*DIPOLE_COEFF(c_coeff,1)/
     '            R**4)*(-ANAL_A*DSIN(THETA)
     '            +ANAL_B*DCOS(THETA))
              ELSE IF(IVAL.EQ.6) THEN !s2 deriv of normal deriv
                V=1.0d0/(2.0d0*PI*SIGMA(1)*R**4)*(ANAL_A*
     '            DCOS(THETA)*DSIN(PHI)+ANAL_B*DSIN(THETA)*
     '            DSIN(PHI)-ANAL_C*DCOS(PHI))
                U=(DIPOLE_COEFF(b_coeff,1)/R-
     '            2.0d0*DIPOLE_COEFF(c_coeff,1)/
     '            R**4)*(-ANAL_A*DCOS(THETA)*DSIN(PHI)-
     '            ANAL_B*DSIN(THETA)*DSIN(PHI)+ANAL_C*DCOS(PHI))
              ELSE IF(IVAL.EQ.7) THEN ! s1&s2 derivative
                V=0.0d0
                U=0.0d0
              ELSE IF(IVAL.EQ.8) THEN ! s1&s2 deriv of normal deriv
                V=0.0d0
                U=0.0d0
              ELSE
                ERROR='>>Invalid IVAL number'
                GOTO 9999
              ENDIF
            ENDIF
            W=V-U
            IF(nr.GT.1) THEN
              IF((IVAL.GE.4).AND.(DABS(R-MESH3_RAD(sphere-1))/
     '          MESH3_RAD(sphere-1).LE.SPHERE_RAD_TOL)) THEN
                W=-W !node is on the interface where the normal
C                     is reversed so change the sign of W.
              ENDIF
            ENDIF
            IF(IVAL.GE.4) W=SIGMA(sphere)*W !adjust normal derivs
C                                            for conductivity
          ELSE IF(ANAL_CHOICE(nr).EQ.15) THEN !Multiple spheres,ecc
            IF(DIPOLE_AXIS.EQ.1) THEN
              ZZ(1)=-Z
              ZZ(2)=Y
              ZZ(3)=X
            ELSE IF(DIPOLE_AXIS.EQ.2) THEN
              ZZ(1)=X
              ZZ(2)=-Z
              ZZ(3)=Y
            ELSE IF(DIPOLE_AXIS.EQ.3) THEN
              ZZ(1)=X
              ZZ(2)=Y
              ZZ(3)=Z
            ELSE
              ERROR='>>Invalid dipole axis'
              GOTO 9999
            ENDIF
            CALL ZX(3,ZZ,XX)
            R=XX(1)
            THETA=XX(2)
C           Measure PHI from the z axis
            PHI=PI/2.0d0-XX(3)
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(DIPOLE_EVALUATE_4)
              WRITE(OP_STRING,'('' Nodal coordinates:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' X='',D12.4,'', Y='','
     '          //'D12.4,'', Z='',D12.4)') X,Y,Z
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' R='',D12.4,'', THETA='',D12.4,'
     '          //''', PHI='',D12.4)') R,THETA,PHI
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(DIPOLE_EVALUATE_4)
            ENDIF
            sphere=nr !the sphere number/region the node is in.
            a_coeff=2*sphere-2
            b_coeff=2*sphere-1
            CALL ASSERT(IVAL.GE.1.AND.IVAL.LE.7,
     '        '>>Invalid IVAL',ERROR,*9999)

            IF((IVAL.EQ.1).OR.(IVAL.EQ.4)) THEN
C MLB what on earth was this trying to do - it didn't work anyway
C            IF(IVAL.NE.1.OR.IVAL.NE.4) THEN

            ELSE
C             Calculate the s directions according to coord
C             system with dipole on z axis
              IF(DIPOLE_AXIS.EQ.1) THEN
                s1_dir(1)=0.0d0
                s1_dir(2)=DCOS(THETA)
                s1_dir(3)=DSIN(THETA)
                s2_dir(1)=DSIN(PHI)
                s2_dir(2)=-DCOS(PHI)*DSIN(THETA)
                s2_dir(3)=DCOS(PHI)*DCOS(THETA)
              ELSE IF(DIPOLE_AXIS.EQ.3) THEN
                s1_dir(1)=-DSIN(THETA)
                s1_dir(2)=DCOS(THETA)
                s1_dir(3)=0.0d0
                s2_dir(1)=-DCOS(PHI)*DCOS(THETA)
                s2_dir(2)=-DCOS(PHI)*DSIN(THETA)
                s2_dir(3)=DSIN(PHI)
              ELSE
                ERROR='>>Invalid dipole axis'
                GOTO 9999
              ENDIF
C             Deriv coefficients are dot products
              s1_coeff=0.0d0
              s2_coeff=0.0d0
              IF(IVAL.LE.4) THEN
                DO nj=1,3 !not general
                  s1_coeff=s1_coeff+XP(IVAL,nv,nj)*s1_dir(nj)
                  s2_coeff=s2_coeff+XP(IVAL,nv,nj)*s2_dir(nj)
                ENDDO !nj
              ELSE
                DO nj=1,3 !not general
                  s1_coeff=s1_coeff+XP(IVAL-3,nv,nj)*s1_dir(nj)
                  s2_coeff=s2_coeff+XP(IVAL-3,nv,nj)*s2_dir(nj)
                ENDDO !nj
              ENDIF
            ENDIF
C$OMP PARALLEL DO
C$OMP&  PRIVATE(NSERIES,MSERIES,W_TEMP,index,first_compo,
C$OMP&          second_compo,nj,first_compo_s1,first_compo_s2)
C$OMP&  SHARED(R,THETA,PHI,sphere,a_coeff,b_coeff,s1_dir,s2_dir,
C$OMP&         s1_coeff,s2_coeff,IVAL)
C$OMP&  REDUCTION(+:W)
            DO NSERIES=1,ANAL_N
              DO MSERIES=0,1
                W_TEMP = 0.0d0 !temp
                index=NSERIES*2+MSERIES
C               Potential or flux are simple, gradients require
C               3 components
                IF(IVAL.EQ.1.OR.IVAL.EQ.4) THEN !value
                  IF(MSERIES.EQ.0) THEN
                    first_compo=DBLE(NSERIES)*DIPOLE_STRENGTH(3)
                  ELSE
                    first_compo=-DIPOLE_STRENGTH(1)*DCOS(THETA)-
     '                DIPOLE_STRENGTH(2)*DSIN(THETA)
                  ENDIF
                  IF(IVAL.EQ.1) THEN !value
                    IF(sphere.EQ.1) THEN
                      second_compo=
     '                  DIPOLE_F**(NSERIES-1)/(R**(NSERIES+1))+
     '                  DIPOLE_COEFF(b_coeff,index)*R**NSERIES
                    ELSE
                      second_compo=
     '                  DIPOLE_COEFF(a_coeff,index)/R**(NSERIES+1)+
     '                  DIPOLE_COEFF(b_coeff,index)*R**NSERIES
                    ENDIF
                  ELSE !normal derivative
                    IF(sphere.EQ.1) THEN
                      second_compo=-DBLE(NSERIES+1)*
     '                  DIPOLE_F**(NSERIES-1)/R**(NSERIES+2)+
     '                  DBLE(NSERIES)*DIPOLE_COEFF(b_coeff,index)*
     '                  R**(NSERIES-1)
                    ELSE
                      second_compo=-DBLE(NSERIES+1)*
     '                  DIPOLE_COEFF(a_coeff,index)/R**(NSERIES+2)+
     '                  DBLE(NSERIES)*DIPOLE_COEFF(b_coeff,index)*
     '                  R**(NSERIES-1)
                    ENDIF
                  ENDIF
                  W_TEMP=first_compo*second_compo*
     '              PLGNDR(NSERIES,MSERIES,DCOS(PHI))/
     '              (4.0d0*PI*SIGMA(1))
                ELSE IF(IVAL.LE.6) THEN !derivative (d by dphi)
                  IF(MSERIES.EQ.0) THEN
                    first_compo_s1=0.0d0
                    first_compo_s2=DBLE(NSERIES)*DIPOLE_STRENGTH(3)
                  ELSE
                    first_compo_s1=DIPOLE_STRENGTH(1)*DSIN(THETA)-
     '                DIPOLE_STRENGTH(2)*DCOS(THETA)
                    first_compo_s2=-DIPOLE_STRENGTH(1)*DCOS(THETA)-
     '                DIPOLE_STRENGTH(2)*DSIN(THETA)
                  ENDIF
                  IF(IVAL.LT.4) THEN
                    IF(sphere.EQ.1) THEN
                      second_compo=
     '                  DIPOLE_F**(NSERIES-1)/(R**(NSERIES+2))+
     '                  DIPOLE_COEFF(b_coeff,index)*R**(NSERIES-1)
                    ELSE
                      second_compo=
     '                  DIPOLE_COEFF(a_coeff,index)/R**(NSERIES+2)+
     '                  DIPOLE_COEFF(b_coeff,index)*R**(NSERIES-1)
                    ENDIF
                  ELSE
                    IF(sphere.EQ.1) THEN
                      second_compo=-DBLE(NSERIES+1)*
     '                  DIPOLE_F**(NSERIES-1)/(R**(NSERIES+3))+
     '                  DBLE(NSERIES)*DIPOLE_COEFF(b_coeff,index)*
     '                  R**(NSERIES-2)
                    ELSE
                      second_compo=-DBLE(NSERIES+1)*
     '                  DIPOLE_COEFF(a_coeff,index)/R**(NSERIES+3)+
     '                  DBLE(NSERIES)*DIPOLE_COEFF(b_coeff,index)*
     '                  R**(NSERIES-2)
                    ENDIF
                  ENDIF
                  W_TEMP=(second_compo/(4.0d0*PI*SIGMA(1)))*
     '              (s1_coeff*first_compo_s1*
     '              M_PLGNDR_SIN(NSERIES,MSERIES,DCOS(PHI))-
     '              s2_coeff*first_compo_s2*
     '              (PLGNDR(NSERIES,MSERIES+1,DCOS(PHI))+
     '              DCOS(PHI)*
     '              M_PLGNDR_SIN(NSERIES,MSERIES,DCOS(PHI))))
                ELSE IF(IVAL.EQ.7) THEN !CPB adding cross-deriv
                  IF(MSERIES.EQ.0) THEN
                    W_TEMP=0.0d0
                  ELSE
                    first_compo=DIPOLE_STRENGTH(1)*DSIN(THETA)-
     '                DIPOLE_STRENGTH(2)*DCOS(THETA)
                    IF(sphere.EQ.1) THEN
                      second_compo=
     '                  DIPOLE_F**(NSERIES-1)/(R**(NSERIES+3))+
     '                  DIPOLE_COEFF(b_coeff,index)*R**(NSERIES-2)
                    ELSE
                      second_compo=
     '                  DIPOLE_COEFF(a_coeff,index)/R**(NSERIES+3)+
     '                  DIPOLE_COEFF(b_coeff,index)*R**(NSERIES-2)
                    ENDIF
                    IF(DABS(DSIN(PHI)).GT.ZERO_TOL) THEN
                      W_TEMP=(second_compo/(4.0d0*PI*SIGMA(1)))*(
     '                  first_compo*(M_PLGNDR_SIN(NSERIES,MSERIES+1,
     '                  DCOS(PHI))/DBLE(MSERIES+1)+
     '                  (DCOS(PHI)/DSIN(PHI))*
     '                  M_PLGNDR_SIN(NSERIES,MSERIES,DCOS(PHI))))
                    ELSE
                      W_TEMP=0.0d0
                    ENDIF
                  ENDIF
                ENDIF
                IF(NP_INTERFACE(np,1).NE.nr) THEN
                  IF((IVAL.GE.4.AND.IVAL.LE.6).OR.IVAL.EQ.8) THEN
                    W_TEMP = -W_TEMP !Reverse the normal
                  ENDIF
                ENDIF
                IF((IVAL.GE.4.AND.IVAL.LE.6).OR.IVAL.EQ.8)
     '            W_TEMP=SIGMA(sphere)*W_TEMP !adjust norm derivatives
C                                                for conductivity
                W=W+W_TEMP
              ENDDO !MSERIES
            ENDDO !NSERIES
          ENDIF !anal_choice
        ENDIF !NJT=2/3
      ELSE !ANAL_CHOICE=14 (eccentric dipole)
        IF(IVAL.LE.3) THEN
          CALL ASSERT(NJT.EQ.3,'>>Must be 3D for this solution',
     '      ERROR,*9999)
          IF(ITYP3(nr,nx).EQ.2) THEN !generalised Laplace
            COND=CE(1)
          ELSE
            COND=1.0d0
          ENDIF
C           find the dipole
          CALL ASSERT(NDIPOLES(nr).EQ.1,
     '      '>>Must have only one dipole',ERROR,*9999)
          ndip=1
          CALL ASSERT(DIPOLE_CEN_NTIME(ndip,nr).EQ.0,
     '      '>>No time dependence for position allowed',ERROR,*9999)
          CALL ASSERT(DIPOLE_DIR_NTIME(ndip,nr).EQ.0,
     '      '>>No time dependence for direction allowed',ERROR,*9999)
          POSITION(1)=X
          POSITION(2)=Y
          POSITION(3)=Z
          CALL ECCENTRIC_DIPOLE(COND,DIPOLE_DIR(1,0,ndip,nr),
     '      DIPOLE_CEN(1,0,ndip,nr),POSITION,POTENTIAL,
     '      GRAD_POTENTIAL,ERROR,*9999)
          IF(IVAL.EQ.1) THEN
            W=POTENTIAL
          ELSEIF((IVAL.EQ.2).OR.(IVAL.EQ.3)) THEN
            W=0.0d0
            DO nj=1,3 !Not general
              W=W+XP(IVAL,nv,nj)*GRAD_POTENTIAL(nj) !s1, s2 deriv.
            ENDDO
          ENDIF !IVAL
        ELSE
          W=0.0d0
        ENDIF !IVAL
      ENDIF !not single sphere

C cpb 30/10/96 Adding normal correction

      IF((NP_INTERFACE(np,0).GT.1).AND.(nr.NE.NP_INTERFACE(np,1))) THEN
C       Interface node in slave region
        IF(NJT.EQ.2.AND.IVAL.EQ.4) THEN
C         Reverse sign of normal at the interface
          W=-W
        ELSE IF(NJT.EQ.3.AND.(IVAL.EQ.4.OR.IVAL.EQ.6)) THEN
C         Reverse sign of normal and s2 tangent at the interface
          W=-W
        ENDIF ! njt
      ENDIF ! np_interface


      CALL EXITS('DIPOLE_EVALUATE')
      RETURN
 9999 CALL ERRORS('DIPOLE_EVALUATE',ERROR)
      CALL EXITS('DIPOLE_EVALUATE')
      RETURN 1
      END


