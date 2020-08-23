      SUBROUTINE D_ZGTG53(PARAMTYPE,nb,NMNO,nr,nx,
     '  AXU,AZ,AZL,AZU,CG,D_AZ,D_AZL,D_AZU,D_EG,D_RI3,D_TG,D_ZG,
     '  EG,ZG,ERROR,*)

C#### Subroutine: D_ZGTG53
C###  Description:
C###    D_ZGTG53 evaluates derivatives of components of 2nd
C###    Piola-Kirchhoff stress tensor TG wrt material parameters or
C###    geometric variables at current Gauss point for 3-dimensional
C###    problems (ktyp51(nr)=3).

C**** The following code was copied from ZGTG53 and altered and should
C**** be kept in sych with ZGTG53.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER nb,NMNO,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),D_AZ,D_AZL(3,3),
     '  D_AZU(3,3),D_EG(3,3),D_RI3,D_TG(3,3),D_ZG(NHM,NUM),EG(3,3),
     '  ZG(NHM,NUM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER i,j,m,n
      REAL*8 AXL(3,3),AZL_tmp(3,3),D_DW(6),D_W3,Fgrowth(3,3),W3

      CALL ENTERS('D_ZGTG53',*9999)

      CALL ASSERT(KTYP55(nr).EQ.3,'Stresses must be a function of'
     '  //' fibre and transverse strains in the constitutive law',
     '  ERROR,*9999)
      CALL ASSERT(KTYP53(nr).GT.1,'>>Stresses must be wrt body/fibre '
     '  //'coords',ERROR,*9999)

C     Compute undeformed metric tensors
      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
          AXL(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
        AXL(i,i)=1.0d0
      ENDDO

C!!! NOTE can only use resid strains for pole zero law until init extn
C!!!      mat params are set up for other problem types
      IF(KTYP56(nr).EQ.3) THEN !pole zero law
C       Calc growth defm tens for resid strain and copy AZL to AZL_tmp
        DO i=1,3
          DO j=1,3
            Fgrowth(i,j)=0.0d0 !growth defm tensor for resid strains
            AZL_tmp(i,j)=AZL(i,j)
          ENDDO !j
          Fgrowth(i,i)=CG(27+i) !NOTE init extns CG(28),CG(29),CG(30)
        ENDDO !i
C       Apply growth defm to deformed covariant metric tensor AZL
        DO i=1,3
          DO j=1,3
            AZL(i,j)=0.0d0
            DO m=1,3
              DO n=1,3
                AZL(i,j)=AZL(i,j)+
     '            Fgrowth(i,m)*AZL_tmp(m,n)*Fgrowth(n,j)
              ENDDO !n
            ENDDO !m
          ENDDO !j
        ENDDO !i
C       Recompute AZU,AZ from transformed AZL
        CALL INVERT(NIT(nb),AZL,AZU,AZ)
      ENDIF !pole zero law
C     Compute Green strain tensor wrt material fibre coords
      DO i=1,3
        DO j=1,3
          EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
        ENDDO !j
      ENDDO !i
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' EG: '',3D12.4,/(5X,3D12.4))')
     '    ((EG(i,j),j=1,3),i=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
        D_RI3=0.0d0
C       Calc derivs of strain energy fn wrt material parameters
        CALL D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,EG(1,1),EG(2,2),EG(3,3),
     '    EG(1,2),EG(1,3),EG(2,3),ERROR,*9999)
C       Calculate derivs of TG wrt material parameters
        DO i=1,3
          D_TG(i,i)=D_DW(i)/AXL(i,i)
        ENDDO
        D_TG(1,2)=D_DW(4)/DSQRT(AXL(1,1)*AXL(2,2))
        D_TG(1,3)=D_DW(5)/DSQRT(AXL(1,1)*AXL(3,3))
        D_TG(2,3)=D_DW(6)/DSQRT(AXL(2,2)*AXL(3,3))
        D_TG(2,1)=D_TG(1,2)
        D_TG(3,1)=D_TG(1,3)
        D_TG(3,2)=D_TG(2,3)

      ELSE IF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN

        IF(KTYP56(nr).EQ.3) THEN !pole zero law
          CALL ASSERT(CG(28).EQ.1.0d0.AND.CG(29).EQ.1.0d0.AND.
     '      CG(30).EQ.1.0d0,'>>Needs updating for growth tensor',
     '      ERROR,*9999)
        ENDIF

C       Set the derivative of the third strain invariant wrt the
C       current geometric variable
        D_RI3=D_AZ
C       Calculate the derivative of the strain matrix, D_EG, wrt the
C       current geometric variable at the current Gauss point
C       Here D_AZL is the derivative of deformed metric tensor
C       (Nu coords) wrt the current geometric variable
        DO i=1,3
          DO j=1,3
            D_EG(i,j)=0.50d0*D_AZL(i,j)/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' D_EG: '',9D12.4)')
     '      ((D_EG(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C       Get deriv's of strain energy fn wrt fibre strains
        CALL D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,EG(1,1),EG(2,2),EG(3,3),
     '    EG(1,2),EG(1,3),EG(2,3),ERROR,*9999)

        IF(KTYP52(nr).EQ.1.OR.KTYP52(nr).EQ.4.OR.KTYP52(nr).EQ.6) THEN !compressible
          W3=0.0d0   !Check this
          D_W3=0.0d0 !Check this
        ELSE                                        !incompressible
          W3=ZG(NH_LOC(0,nx),1)
          D_W3=D_ZG(NH_LOC(0,nx),1)
        ENDIF

C       Get deriv's of stress tensor wrt geometric variables
        DO i=1,3
          D_TG(i,i)=D_DW(i)*D_EG(i,i)/AXL(i,i)
     '            + 2.0d0*D_W3*AZU(i,i) + 2.0d0*W3*D_AZU(i,i)
        ENDDO
        D_TG(1,2)=D_DW(4)*D_EG(1,2)/DSQRT(AXL(1,1)*AXL(2,2))
     '            + 2.0d0*D_W3*AZU(1,2) + 2.0d0*W3*D_AZU(1,2)
        D_TG(1,3)=D_DW(5)*D_EG(1,3)/DSQRT(AXL(1,1)*AXL(3,3))
     '            + 2.0d0*D_W3*AZU(1,3) + 2.0d0*W3*D_AZU(1,3)
        D_TG(2,3)=D_DW(6)*D_EG(2,3)/DSQRT(AXL(2,2)*AXL(3,3))
     '            + 2.0d0*D_W3*AZU(2,3) + 2.0d0*W3*D_AZU(2,3)
        D_TG(2,1)=D_TG(1,2)
        D_TG(3,1)=D_TG(1,3)
        D_TG(3,2)=D_TG(2,3)
      ENDIF

      CALL EXITS('D_ZGTG53')
      RETURN
 9999 CALL ERRORS('D_ZGTG53',ERROR)
      CALL EXITS('D_ZGTG53')
      RETURN 1
      END


