      SUBROUTINE D_ZGTG54(PARAMTYPE,NMNO,nr,
     '  AXU,AZ,AZL,AZU,CG,D_AZ,D_AZL,D_AZU,D_EG,D_RI3,D_TG,
     '  EG,YG,ERROR,*)

C#### Subroutine: D_ZGTG54
C###  Description:
C###    D_ZGTG54 evaluates derivatives of components of 2nd
C###    Piola-Kirchhoff stress tensor TG wrt material parameters or
C###    geometric variables at current Gauss point for membrane
C###    problems (ktyp51(nr)=4).

C**** The following code was copied from ZGTG54 and altered and should
C**** be kept in sych with ZGTG54.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp50.cmn'
!     Parameter List
      INTEGER NMNO,nr
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),D_AZ,D_AZL(3,3),
     '  D_AZU(3,3),D_EG(3,3),D_RI3,D_TG(3,3),EG(3,3),YG(NIYGM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER i,ie,j,je
      REAL*8 AXL(3,3),D_DW(6),D_H,DW(6),H

      CALL ENTERS('D_ZGTG54',*9999)

      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
          AXL(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
        AXL(i,i)=1.0d0
      ENDDO

      CALL ASSERT(KTYP55(nr).EQ.3,'Stresses must be a function of'
     '  //' fibre  and transverse strains in the constitutive law',
     '  ERROR,*9999)

      IF(KTYP53(nr).EQ.1) THEN      !stresses in theta coords
        ERROR='>>Stresses must be wrt body/fibre coords'
        GO TO 9999
      ELSE IF(KTYP53(nr).GT.1) THEN !stresses in body/fibre Nu coords
        AZL(3,3)=1.0d0/AZ !AZ=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        AZU(3,3)=AZ
        AZL(1,3)=0.0d0
        AZL(2,3)=0.0d0
        AZL(3,1)=0.0d0
        AZL(3,2)=0.0d0
        IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
          D_AZL(3,3)=0.0d0
          D_AZU(3,3)=0.0d0
        ELSEIF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
          D_AZL(3,3)=-D_AZ/(AZ*AZ)
          D_AZU(3,3)=D_AZ
        ENDIF
        D_AZL(1,3)=0.0d0
        D_AZL(2,3)=0.0d0
        D_AZL(3,1)=0.0d0
        D_AZL(3,2)=0.0d0
        D_RI3=0.0d0  !wrt both material and geometric parameters
      ENDIF

      !Calculate strain matrix, EG, at current Gauss point
      !Here AZL is the deformed metric tensor wrt Nu coords
      DO ie=1,3
        DO je=1,3
          EG(ie,je)=0.50d0*(AZL(ie,je)-AXL(ie,je))/
     '      DSQRT(AXL(ie,ie)*AXL(je,je))
        ENDDO
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' EG: '',9D12.4)')
     '    ((EG(ie,je),je=1,3),ie=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(PARAMTYPE.EQ.'MATERIAL_PARAMETERS') THEN
        !Calc derivs of strain energy fn wrt material parameters
        CALL D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,EG(1,1),EG(2,2),EG(3,3),
     '    EG(1,2),EG(1,3),EG(2,3),ERROR,*9999)
        !Calculate derivs of TG wrt material parameters
        IF(KTYP52(nr).EQ.1.OR.KTYP52(nr).EQ.4.OR.KTYP52(nr).EQ.6) THEN !compressible
          !22Jan89: Compressible case should also be included here
        ELSE                                        !incompressible
          D_H=-D_DW(3)/AZU(3,3)
        ENDIF
        D_TG(1,1)=D_DW(1)/AXL(1,1)                +D_H*AZU(1,1)
        D_TG(2,2)=D_DW(2)/AXL(2,2)                +D_H*AZU(2,2)
        D_TG(1,2)=D_DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+D_H*AZU(1,2)
        D_TG(2,1)=D_TG(1,2)
        D_TG(1,3)=0.0d0
        D_TG(2,3)=0.0d0
        D_TG(3,1)=0.0d0
        D_TG(3,2)=0.0d0

      ELSEIF(PARAMTYPE.EQ.'GEOMETRIC_PARAMETERS') THEN
        !Calculate the derivative of the strain matrix, D_EG, wrt the
        !current geometric variable at the current Gauss point
        !Here D_AZL is the derivative of deformed metric tensor
        !(Nu coords) wrt the current geometric variable
        DO ie=1,3
          DO je=1,3
            D_EG(ie,je)=0.50d0*D_AZL(ie,je)/DSQRT(AXL(ie,ie)*AXL(je,je))
          ENDDO
        ENDDO
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' D_EG: '',9D12.4)')
     '      ((D_EG(ie,je),je=1,3),ie=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
C news VJ 10Dec2003 Added YG to ENERGY param list
        CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     '    EG(2,3),YG,ERROR,*9999)

        !Get deriv's of strain energy fn wrt fibre strains
        CALL D_ENERGY(PARAMTYPE,NMNO,nr,CG,D_DW,EG(1,1),EG(2,2),EG(3,3),
     '    EG(1,2),EG(1,3),EG(2,3),ERROR,*9999)

        IF(KTYP52(nr).EQ.1) THEN !compressible
          !22Jan89: Compressible case should also be included here
        ELSE                 !incompressible
          H= -DW(3)/AZU(3,3)
          D_H= -D_DW(3)*D_EG(3,3)/AZU(3,3)
     '         +DW(3)*D_AZU(3,3)/(AZU(3,3)*AZU(3,3))
        ENDIF

        !Get deriv's of stress tensor wrt geometric variables
        D_TG(1,1)=D_DW(1)*D_EG(1,1)/AXL(1,1)
     '    + D_H*AZU(1,1) + H*D_AZU(1,1)
        D_TG(2,2)=D_DW(2)*D_EG(2,2)/AXL(2,2)
     '    + D_H*AZU(2,2) + H*D_AZU(2,2)
        D_TG(1,2)=D_DW(4)*D_EG(1,2)/DSQRT(AXL(1,1)*AXL(2,2))
     '    + D_H*AZU(1,2) + H*D_AZU(1,2)
        D_TG(2,1)=D_TG(1,2)
        D_TG(1,3)=0.0d0
        D_TG(2,3)=0.0d0
        D_TG(3,1)=0.0d0
        D_TG(3,2)=0.0d0
      ENDIF

      CALL EXITS('D_ZGTG54')
      RETURN
 9999 CALL ERRORS('D_ZGTG54',ERROR)
      CALL EXITS('D_ZGTG54')
      RETURN 1
      END


