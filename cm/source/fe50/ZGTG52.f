      SUBROUTINE ZGTG52(nb,nr,nx,AXU,AZ,AZL,AZU,CG,RI1,RI2,RI3,
     '  TG,YG,ZG,ERROR,*)

C#### Subroutine: ZGTG52
C###  Description:
C###    ZGTG52 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for plane stress problems
C###    (KTYP51(nr)=1).
C**** Note: Stress components should be referred to (undeformed)
C****       Nu (body/fibre) coords (KTYP53(nr)>1).
C**** TGP are principal second Piola Kirchhoff stresses

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER nb,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),RI1,RI2,RI3,
     '  TG(3,3),YG(NIYGM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,mi,mj,ni,NITB,nj
      REAL*8 AXL(3,3),BG(3,3),DW(6),EDG(3,3),EG(3,3),EG12,EG13,EVAL(3),
     '  RK1,RK2,RL1,RL2,RL3,RM(3,3),TGP(3,3),VALTMP,W3

      CALL ENTERS('ZGTG52',*9999)

      IF(KTYP53(nr).EQ.1) THEN
        ERROR='Stresses must be referred to body/fibre coords'
        GOTO 9999
      ENDIF

      NITB=NIT(nb)
      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
      ENDDO

      RI1=AZL(1,1)+AZL(2,2)
      RI2=AZ
C     RI3=(AZU(1,1)+AZU(2,2))*RI3

      IF(KTYP55(nr).EQ.1) THEN
C ***   Invariants K1,K2 are fns of Physical strain that are
C ***   isotropic in the plane transverse to the Nu(1) axis
        EG13=0.0d0
        EG12=AZL(1,2)/2.0d0
        RK1=(AZL(1,1)-1.0d0)/2.0d0
        RK2=EG13*EG13+EG12*EG12
        BG(1,1)=RI1-AZL(1,1)
        BG(2,2)=RI1-AZL(2,2)
        BG(3,3)=RI1-1.0d0
        BG(2,1)=-AZL(2,1)
        BG(3,1)=0.0d0
        BG(3,2)=0.0d0
        CALL ENERGY(nr,CG,DW,RI1,RI2,RI3,RK1,RK2,0.0d0,YG,ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN
          W3=RI3*DW(3)
        ELSE IF(KTYP52(nr).GT.1) THEN
          W3=ZG(NH_LOC(0,nx),1)
        ENDIF
        DO mj=1,3
          DO nj=1,mj
            TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj)
     '        +W3*AZU(mj,nj))
          ENDDO
        ENDDO
        TG(3,1)=TG(3,1)+EG13*DW(5)
        TG(2,1)=TG(2,1)+EG12*DW(5)
        TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
        TG(1,2)=TG(2,1)
        TG(1,3)=TG(3,1)
        TG(2,3)=TG(3,2)

      ELSE IF(KTYP55(nr).EQ.2) THEN
C ***   TG is a function of the principal stretches.
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
          ENDDO
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO
        DO i=1,2
          DO j=1,2
            EDG(i,j)=0.5d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '('' EDG: '',9D12.4)') ((EDG(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C        CALL F02ABF(EDG,3,NITB,EVAL,RM,3,WK1_LOCAL,IFAIL)
C          DO ni=1,3
C            DO mi=1,3
C              RM(ni,mi)=EDG(ni,mi)
C            ENDDO
C          ENDDO
C MLB 19/3/97
C This may not give evectors as accurately as NAG
C          CALL DSYEV('V','L',NITB,RM,3,EVAL,WK1_LOCAL,3,IFAIL)

        CALL ESOLVE(.TRUE.,NITB,EDG,3,EVAL,RM,ERROR,*9999)

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
        RL3=1.0d0
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
        ELSE IF(KTYP52(nr).GT.1) THEN
          TGP(1,1)=DW(1)/RL1+2.0d0*ZG(NH_LOC(0,nx),1)/RL1**2
          TGP(2,2)=DW(2)/RL2+2.0d0*ZG(NH_LOC(0,nx),1)/RL2**2
          TGP(3,3)=DW(3)/RL3+2.0d0*ZG(NH_LOC(0,nx),1)/RL3**2
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' TGP:'',12X,3D11.3,/(16X,3D11.3))')
     '      ((TGP(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        CALL TRAN(3,TGP,TG,RM)

      ELSE IF(KTYP55(nr).EQ.3) THEN
C ***   TG is a function of the fibre and transverse strains
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
            EG(i,j)=0.0d0
          ENDDO
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO
        DO i=1,2
          DO j=1,2
            EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' EG: '',9D12.4)')
     '      ((EG(i,j),j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     '    EG(2,3),YG,ERROR,*9999)
        IF(KTYP52(nr).EQ.1) THEN
          W3=0.0d0 !Check this
        ELSE
          W3=ZG(NH_LOC(0,nx),1)
        ENDIF
        DO i=1,3
          TG(i,i)=DW(i)/AXL(i,i)+2.0d0*W3*AZU(i,i)
        ENDDO
        TG(2,1)=DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+2.0d0*W3*AZU(2,1)
        TG(3,1)=DW(5)/DSQRT(AXL(1,1)*AXL(3,3))+2.0d0*W3*AZU(3,1)
        TG(3,2)=DW(6)/DSQRT(AXL(2,2)*AXL(3,3))+2.0d0*W3*AZU(3,2)
        TG(1,2)=TG(2,1)
        TG(1,3)=TG(3,1)
        TG(2,3)=TG(3,2)
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' TG:'',12X,3D12.4,/(16X,3D12.4))')
     '    ((TG(i,j),J=1,3),I=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(4,'('' TG:'',12X,3D12.4,/(16X,3D12.4))')
     '    ((TG(i,j),J=1,3),I=1,3)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZGTG52')
      RETURN
 9999 CALL ERRORS('ZGTG52',ERROR)
      CALL EXITS('ZGTG52')
      RETURN 1
      END


