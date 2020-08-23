      SUBROUTINE ZGTG54(nb,nr,AXU,AZ,AZL,AZU,CG,EG,RI1,RI2,RI3,
     '  TG,YG,ERROR,*)

C#### Subroutine: ZGTG54
C###  Description:
C###    ZGTG54 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for membrane problems
C###    (KTYP51(nr)=4).
C**** Stress components must always be with respect to
C**** (undeformed) Nu (body/fibre) coords (KTYP53(nr)>1).
C**** Note: For transversely isotropic strain energy functions of K1,K2
C****       the isotropic plane is transverse to Nu(1) (KTYP53(nr)>1)
C****       For transv isotropic strain energy functions of strains EG,
C****       the isotropic plane is transverse to Nu(1) (KTYP53(nr)>1)
C**** TGP are principal second Piola Kirchhoff stresses

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER nb,nr
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),EG(3,3),RI1,RI2,RI3,
     '  TG(3,3),YG(NIYGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,m,mi,mj,n,ni,NITB,nj
      REAL*8 AXL(3,3),AZL_tmp(3,3),BG(3,3),DW(6),EG12,EG13,EVAL(3),
     '  Fgrowth(3,3),H,RK1,RK2,RL1,RL2,RL3,RM(3,3),
     '  TGP(3,3),VALTMP,W3

      CALL ENTERS('ZGTG54',*9999)
      NITB=NIT(nb)

      DO i=1,3
        DO j=1,3
          AXL(i,j)=0.0d0
          AXU(i,j)=0.0d0
          TG(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
        AXL(i,i)=1.0d0
      ENDDO

      IF(KTYP53(nr).EQ.1) THEN      !stresses in theta coords
        ERROR='>>Stresses must be wrt body/fibre coords'
        GO TO 9999
      ELSE IF(KTYP53(nr).GT.1) THEN !stresses in body/fibre Nu coords
        AZL(3,3)=1.0d0/AZ !AZ=azl(1,1)*azl(2,2)-azl(1,2)*azl(2,1)
        AZU(3,3)=AZ
        AZL(1,3)=0.0d0
        AZL(2,3)=0.0d0
        AZL(3,1)=0.0d0
        AZL(3,2)=0.0d0
C       RI3= AZ    !????? I3=det(AZL)=AZ/AZ=1.0d0
        RI3=1.0d0    !13Aug90
        RI1= AZL(1,1)+AZL(2,2)+AZL(3,3)
        RI2=(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3
      ENDIF

      IF(KTYP55(nr).EQ.1) THEN !TG is fn of the principal invariants
        IF(KTYP53(nr).GT.1) THEN
C         Invariants K1,K2 are fns of Physical strain that are
C         isotropic in the plane transverse to the Nu(1) axis
          EG13=AZL(1,3)/2.0d0
          EG12=AZL(1,2)/2.0d0
          RK1=(AZL(1,1)-1.0d0)/2.0d0
          RK2=EG13*EG13+EG12*EG12
          BG(1,1)=RI1-AZL(1,1)
          BG(2,2)=RI1-AZL(2,2)
          BG(3,3)=RI1-AZL(3,3)
          BG(2,1)=   -AZL(2,1)
          BG(3,1)=   -AZL(3,1)
          BG(3,2)=   -AZL(3,2)
          CALL ENERGY(nr,CG,DW,RI1,RI2,RI3,RK1,RK2,0.0d0,YG,ERROR,*9999)
          IF(KTYP52(nr).EQ.1) THEN !compressible
            W3=RI3*DW(3)
          ELSE                 !incompressible
            W3=-(DW(1)+DW(2)*BG(3,3))*AZL(3,3)  !hydro press=2*W3
          ENDIF
          DO mj=1,NITB
            DO nj=1,mj
              TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)+W3
     '          *AZU(mj,nj))
            ENDDO
          ENDDO
          TG(3,3)=0.0d0
          TG(3,1)=TG(3,1)+EG13*DW(5)
          TG(2,1)=TG(2,1)+EG12*DW(5)
          TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
          TG(1,2)=TG(2,1)
          TG(1,3)=TG(3,1)
          TG(2,3)=TG(3,2)
        ENDIF

      ELSE IF(KTYP55(nr).EQ.2) THEN !TG is fn of the principal stretches
        DO i=1,3
          DO j=1,3
            EG(i,j)=0.5d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO !j
        ENDDO !i
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' EG:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((EG(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C        CALL F02ABF(EG,3,NITB,EVAL,RM,3,WK1_LOCAL,IFAIL)
C          DO ni=1,3
C            DO mi=1,3
C              RM(ni,mi)=EG(ni,mi)
C            ENDDO
C          ENDDO
C MLB 19/3/97
C This may not give evectors as accurately as NAG
C          CALL DSYEV('V','L',NITB,RM,3,EVAL,WK1_LOCAL,3,IFAIL)

        CALL ESOLVE(.TRUE.,NITB,EG,3,EVAL,RM,ERROR,*9999)

C news MPN 30-Jun-94
C ***   Order EVAL from max->min and change cols of RM accordingly
        IF(EVAL(2).GT.EVAL(1)) THEN
          VALTMP=EVAL(1)
          EVAL(1)=EVAL(2)
          EVAL(2)=VALTMP
          DO ni=1,NITB
            VALTMP=RM(ni,1)
            RM(ni,1)=RM(ni,2)
            RM(ni,2)=VALTMP
          ENDDO
        ENDIF
        IF(NITB.EQ.3) THEN
          DO mi=1,2
            IF(EVAL(3).GT.EVAL(MI)) THEN
              VALTMP=EVAL(MI)
              EVAL(MI)=EVAL(3)
              EVAL(3)=VALTMP
              DO ni=1,NITB
                VALTMP=RM(ni,mi)
                RM(ni,mi)=RM(ni,3)
                RM(ni,3)=VALTMP
              ENDDO
            ENDIF
          ENDDO
        ENDIF
C newe
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' EVAL:'',7X,(3X,3D11.3))')
     '      (EVAL(i),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RM:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((RM(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        RL1=DSQRT(2.0d0*EVAL(1)+1.0d0)
        RL2=DSQRT(2.0d0*EVAL(2)+1.0d0)
        IF(KTYP52(nr).EQ.1) THEN      !compressible
          RL3=1.0d0
        ELSE IF(KTYP52(nr).GT.1) THEN !incompressible
          RL3=1.0d0/(RL1*RL2)
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' RL1,RL2,RL3:'',(3X,3D11.3))')
     '      RL1,RL2,RL3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        DO i=1,3
          DO j=1,3
            TGP(i,j)=0.0d0
          ENDDO
        ENDDO
        CALL ENERGY(nr,CG,DW,RL1,RL2,RL3,0.0d0,0.0d0,0.0d0,YG,
     '    ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN
          TGP(1,1)=DW(1)/RL1
          TGP(2,2)=DW(2)/RL2
          TGP(3,3)=DW(3)/RL3
        ELSE IF(KTYP52(nr).GT.1)  THEN
          RL3=1.0d0/(RL1*RL2)
          TGP(1,1)=DW(1)/RL1
          TGP(2,2)=DW(2)/RL2
          TGP(3,3)=DW(3)/RL3
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' TGP:'',12X,3D11.3,/(16X,3D11.3))')
     '       ((TGP(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        CALL TRAN(3,TGP,TG,RM)

      ELSE IF(KTYP55(nr).EQ.3) THEN !TG is fn of the fibre trans strains
C!!! NOTE can only use resid strains for pole zero law until init extn
C!!!      mat params are set up for other problem types
        IF(KTYP56(nr).EQ.3) THEN !pole zero law
C         Calc growth defm tens for resid strain and copy AZL to AZL_tmp
          DO i=1,3
            DO j=1,3
              Fgrowth(i,j)=0.0d0 !growth defm tensor for resid strains
              AZL_tmp(i,j)=AZL(i,j)
            ENDDO !j
            Fgrowth(i,i)=CG(27+i) !NOTE init extns CG(28),CG(29),CG(30)
          ENDDO !i
C         Apply growth defm to deformed covariant metric tensor AZL
          DO i=1,3
            DO j=1,3
              AZL(i,j)=0.0d0
              DO m=1,3
                DO n=1,3
                  AZL(i,j)=AZL(i,j)+
     '              Fgrowth(i,m)*AZL_tmp(m,n)*Fgrowth(n,j)
                ENDDO !n
              ENDDO !m
            ENDDO !j
          ENDDO !i
C         Recompute AZU,AZ from transformed AZL
          CALL INVERT(NIT(nb),AZL,AZU,AZ)
        ENDIF !KTYP56=3 (pole-zero law)
        DO i=1,3
          DO j=1,3
            EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO !j
        ENDDO !i
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' EG:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((EG(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(KTYP56(nr).EQ.6) THEN !linear viscous relation
          write(*,'('' linear viscous relation'')')
          TG(1,1)=CG(1)*EG(1,1)
          TG(2,2)=CG(2)*EG(2,2)
          TG(1,2)=CG(3)*EG(1,2)
          TG(2,1)=TG(1,2)
          TG(1,3)=0.0d0
          TG(2,3)=0.0d0
          TG(3,1)=0.0d0
          TG(3,2)=0.0d0
          TG(3,3)=0.0d0

        ELSE IF(KTYP56(nr).NE.6) THEN !all others
          CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     '      EG(2,3),YG,ERROR,*9999)
          IF(KTYP52(nr).EQ.1) THEN !compressible
C ***     22Jan89: Compressible case should also be included here
          ELSE !incompressible
            H=-DW(3)/AZU(3,3)
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
            WRITE(OP_STRING,'('' DW(1..4):'',4D11.3,'' H='',D11.3)')
     '        (DW(i),i=1,4),H
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
          ENDIF
          TG(1,1)=DW(1)/AXL(1,1)                +H*AZU(1,1)
          TG(2,2)=DW(2)/AXL(2,2)                +H*AZU(2,2)
          TG(1,2)=DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+H*AZU(1,2)
          TG(2,1)=TG(1,2)
          TG(1,3)=0.0d0
          TG(2,3)=0.0d0
          TG(3,1)=0.0d0
          TG(3,2)=0.0d0
        ENDIF !ktyp56
      ENDIF !ktyp55

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' TG:'',12X,3D12.4,/(16X,3D12.4))')
     '    ((TG(i,j),j=1,3),i=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZGTG54')
      RETURN
 9999 CALL ERRORS('ZGTG54',ERROR)
      CALL EXITS('ZGTG54')
      RETURN 1
      END


