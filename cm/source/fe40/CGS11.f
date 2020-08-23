      SUBROUTINE CGS11(nr,nw,nx,CG,EG,RT,TG,ZG,ERROR,*)

C#### Subroutine: CGS11
C###  Description:
C###    CGS11 evaluates stress (TG) and strain (EG) for plane stress
C###    and plane strain.

      IMPLICIT NONE
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nr,nw,nx
      REAL*8 CG(NMM),EG(3,3),RT(3,3),TG(3,3),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,k,NU1(0:3)
      REAL*8 CC(3,3),DENS,DIL,FACTOR,POIS,POIS1,POIS2,RLAMDA,RMU,
     '  YMOD,YMOD1,YMOD2

      DATA NU1/1,2,4,7/

      CALL ENTERS('CGS11',*9999)

C cpb 13/5/96 Old way
C      EG(1,1)=dX_dNu(1,1)*ZG(1,NU1(1))+dX_dNu(2,1)*ZG(2,NU1(1))
C      EG(2,2)=dX_dNu(1,2)*ZG(1,NU1(2))+dX_dNu(2,2)*ZG(2,NU1(2))
C      EG(1,2)=0.5d0*(dX_dNu(1,1)*ZG(1,NU1(2))+dX_dNu(2,1)*ZG(2,NU1(2))
C     '             + dX_dNu(1,2)*ZG(1,NU1(1))+dX_dNu(2,2)*ZG(2,NU1(1)))
C      EG(2,1)=EG(1,2)
      DO i=1,2
        DO j=1,2
          EG(i,j)=0.0d0
          DO k=1,2
            EG(i,j)=EG(i,j)+0.5d0*
     '        (RT(k,i)*ZG(k,NU1(j))+RT(k,j)*ZG(k,NU1(i)))
          ENDDO !k
        ENDDO !j
      ENDDO !i

      IF(nw.EQ.11) THEN !plane stress
        IF(IMT(nw).EQ.1) THEN !isotropic
          YMOD=CG(1)  !is Young's modulus
          POIS=CG(2)  !is Poisson's ratio nu
          RLAMDA=YMOD*POIS/((1.d0+POIS)*(1.d0-2.d0*POIS))
          RMU   =YMOD/(2.d0*(1.d0+POIS))
C CPB 26/9/93 moved EG(3,3) calc before DIL calc. ???
          EG(3,3)=-POIS/(1.d0-POIS)*(EG(1,1)+EG(2,2))
          DIL=EG(1,1)+EG(2,2)+EG(3,3)
          TG(1,1)=RLAMDA*DIL+2.d0*RMU*EG(1,1)
          TG(2,2)=RLAMDA*DIL+2.d0*RMU*EG(2,2)
          TG(1,2)=2.d0*RMU*EG(1,2)

        ELSE IF(IMT(nw).EQ.2) THEN !transversely isotropic 1
          CALL CGC11(nr,nw,nx,CC,CG,DENS,ERROR,*9999)
          TG(1,1)=CC(1,1)*EG(1,1)+CC(1,2)*EG(2,2)+CC(1,3)*EG(1,2)
          TG(2,2)=CC(2,1)*EG(1,1)+CC(2,2)*EG(2,2)+CC(2,3)*EG(1,2)
          TG(1,2)=CC(3,1)*EG(1,1)+CC(3,2)*EG(2,2)+CC(3,3)*EG(1,2)
          YMOD1=CG(1) !is Young's modulus E1
          YMOD2=CG(2) !is Young's modulus E2
C          SMOD1=CG(3) !is shear modulus G1
          POIS1=CG(4) !is Poisson's ratio nu1
          POIS2=CG(5) !is Poisson's ratio nu2
          EG(3,3)=-POIS1/YMOD1*TG(1,1)-POIS2/YMOD2*TG(2,2)

        ELSE IF(IMT(nw).EQ.3) THEN !transversely isotropic 2
          CALL CGC11(nr,nw,nx,CC,CG,DENS,ERROR,*9999)
          TG(1,1)=CC(1,1)*EG(1,1)+CC(1,2)*EG(2,2)+CC(1,3)*EG(1,2)
          TG(2,2)=CC(2,1)*EG(1,1)+CC(2,2)*EG(2,2)+CC(2,3)*EG(1,2)
          TG(1,2)=CC(3,1)*EG(1,1)+CC(3,2)*EG(2,2)+CC(3,3)*EG(1,2)
          YMOD1=CG(1) !is Young's modulus E1
          YMOD2=CG(2) !is Young's modulus E2
C          SMOD1=CG(3) !is shear modulus G1
          POIS1=CG(4) !is Poisson's ratio nu1
          POIS2=CG(5) !is Poisson's ratio nu2
          EG(3,3)=-POIS1/YMOD1*TG(1,1)-POIS2/YMOD2*TG(2,2)

        ELSE IF(IMT(nw).EQ.4) THEN !orthotropic
        ENDIF
        TG(3,3)=0.d0

      ELSE IF(nw.EQ.12) THEN  !!!NEEDS FIXING FOR TRANSVERSE ISO ETC
        YMOD=CG(1)
        POIS=CG(2)
C cpb 26/3/96 Adding in plain strain
C        RLAMDA=YMOD*POIS/((1.d0+POIS)*(1.d0-2.d0*POIS))
C        RMU   =YMOD/(2.d0*(1.d0+POIS))
C        EG(3,3)=0.d0
C        DIL=EG(1,1)+EG(2,2)+EG(3,3)
C        TG(1,1)=RLAMDA*DIL+2.d0*RMU*EG(1,1)
C        TG(2,2)=RLAMDA*DIL+2.d0*RMU*EG(2,2)
C        TG(3,3)=POIS*(TG(1,1)+TG(2,2))
        IF(IMT(nw).EQ.1) THEN !isotropic
          FACTOR=YMOD/((1.0d0+POIS)*(1.0d0-2.0d0*POIS))
          TG(1,1)=FACTOR*(1.0d0-POIS)*EG(1,1)+FACTOR*POIS*EG(2,2)
          TG(2,2)=FACTOR*POIS*EG(1,1)+FACTOR*(1.0d0-POIS)*EG(2,2)
          TG(1,2)=FACTOR*(1.0d0-2.0d0*POIS)/2.0d0*EG(1,2)
          TG(3,3)=POIS*(TG(1,1)+TG(2,2))
        ELSE
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF
      ENDIF
      TG(2,1)=TG(1,2)
      TG(3,1)=0.d0
      TG(1,3)=0.d0
      TG(3,2)=0.d0
      TG(2,3)=0.d0

      CALL EXITS('CGS11')
      RETURN
 9999 CALL ERRORS('CGS11',ERROR)
      CALL EXITS('CGS11')
      RETURN 1
      END


