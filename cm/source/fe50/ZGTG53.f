      SUBROUTINE ZGTG53(STRESSTYPE,nb,nr,nx,AXU,
     '  AZ,AZL,AZU,CG,EG,RI1,RI2,RI3,TG,XG,YG,ZG,ERROR,*)

C#### Subroutine: ZGTG53
C###  Description:
C###    ZGTG53 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for 3-dimensional problems
C###    (KTYP51(nr)=3). Stress components are wrt (undeformed) theta
C###    coordinates (KTYP53(nr)=1) or (undeformed) Nu (body/fibre)
C###    coordinates (KTYP53(nr)>1).

C**** Note: For transversely isotropic strain energy functions of K1,K2
C****       the isotropic plane is transverse to theta(3) (KTYP53(nr)=1)
C****       or to Nu(1) (KTYP53(nr)>1)
C****       For transv. isotropic strain energy functions of strains
C****       EG, the isotropic plane is transverse to theta(1)
C****       (KTYP53(nr)=1) or to Nu(1) (KTYP53(nr)>1)
C**** TGP are principal second Piola Kirchhoff stresses
      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER nb,nr,nx
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CG(NMM),EG(3,3),RI1,RI2,RI3,
     '  TG(3,3),XG(NJM,NUM),YG(NIYGM),ZG(NHM,NUM),DW(6)
      CHARACTER ERROR*(*),STRESSTYPE*(*)
!     Local Variables
      INTEGER i,j,m,mi,mj,n,nh_pressure,ni,NITB,nj
      REAL*8 AA,AXL(3,3),AZL_tmp(3,3),BG(3,3),D(3,3),DD1(3,3),DD2(3,3),
     $     dNudXi(3,3),EDG(3,3),EG12,EG13,EG31,EG32,EVAL(3),
     $     Fgrowth(3,3),fibre_angle,G1,G3,hyd_press,RC,RK1,RK2,RL1,RL2,
     $     RL3,RL33,RM(3,3),RR,SLX,SMX,TGP(3,3),VALTMP,W3,W4

      CALL ENTERS('ZGTG53',*9999)
C VJ 10Dec2003 added YG in parameter list of ZGTG53
      NITB=NIT(nb)

      DO i=1,3
        DO j=1,3
          AXU(i,j)=0.0d0
        ENDDO
        AXU(i,i)=1.0d0
      ENDDO

      IF(KTYP53(nr).EQ.1) THEN !stresses referred to theta coords
        IF(ITYP10(nr).EQ.1) THEN
          RI3= AZ
          RI1= AZL(1,1)+AZL(2,2)+AZL(3,3)
          RI2=(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3
        ELSE IF(ITYP10(nr).EQ.2) THEN
          RR=XG(1,1)**2
          RI3=AZ/RR
          RI1= AZL(1,1)+AZL(2,2)/RR+AZL(3,3)
          RI2=(AZU(1,1)+AZU(2,2)*RR+AZU(3,3))*RI3
          AXU(2,2)=1.0d0/RR
        ELSE IF(ITYP10(nr).EQ.3) THEN
          RR=XG(1,1)**2
          RC=RR*DCOS(XG(3,1))**2
          RI3=AZ/(RR*RC)
          RI1= AZL(1,1)+AZL(2,2)/RC+AZL(3,3)/RR
          RI2=(AZU(1,1)+AZU(2,2)*RC+AZU(3,3)*RR)*RI3
          AXU(2,2)=1.0d0/RC
          AXU(3,3)=1.0d0/RR
        ELSE IF(ITYP10(nr).EQ.4) THEN
          AA=FOCUS*FOCUS
          SLX=DSINH(XG(1,1))
          SMX=DSIN (XG(2,1))
          G1=AA*(SLX*SLX+SMX*SMX)
          G3=AA* SLX*SLX*SMX*SMX
          RI3=AZ/(G1*G1*G3)
          RI1= (AZL(1,1)+AZL(2,2))/G1+AZL(3,3)/G3
          RI2=((AZU(1,1)+AZU(2,2))*G1+AZU(3,3)*G3)*RI3
          AXU(1,1)=1.0d0/G1
          AXU(2,2)=1.0d0/G1
          AXU(3,3)=1.0d0/G3
        ELSE IF(ITYP10(nr).EQ.5) THEN
        ENDIF

      ELSE IF(KTYP53(nr).GT.1) THEN
        !stress referred to body/fibre Nu coords
        RI3= AZ
        RI1= AZL(1,1)+AZL(2,2)+AZL(3,3)
        RI2=(AZU(1,1)+AZU(2,2)+AZU(3,3))*RI3
      ENDIF


C-TG=fn of principal invariants--------------------------------------

      IF(KTYP55(nr).EQ.1) THEN !TG is a fn of the principal invariants
        IF(KTYP53(nr).EQ.1) THEN !theta coords (isotropic)
          !Invariants K1,K2 are fns of Physical strain that are
          !isotropic in the plane transverse to the theta(3) axis
          IF(ITYP10(nr).EQ.1) THEN
            EG31=AZL(3,1)/2.0d0
            EG32=AZL(3,2)/2.0d0
            RK1=(AZL(3,3)-1.0d0)/2.0d0
            RK2=EG31*EG31+EG32*EG32
            BG(1,1)=RI1-AZL(1,1)
            BG(2,2)=RI1-AZL(2,2)
            BG(3,3)=RI1-AZL(3,3)
            BG(2,1)=   -AZL(2,1)
            BG(3,1)=   -AZL(3,1)
            BG(3,2)=   -AZL(3,2)
          ELSE IF(ITYP10(nr).EQ.2) THEN
            EG31=AZL(3,1)/2.0d0
            EG32=AZL(3,2)/(2.0d0*RR)
            RK1=(AZL(3,3)-1.0d0)/2.0d0
            RK2=(AZL(1,3)*EG31+AZL(2,3)*EG32)/2.0d0
            BG(1,1)= RI1-AZL(1,1)
            BG(2,2)=(RI1-AZL(2,2)/RR)/RR
            BG(3,3)= RI1-AZL(3,3)
            BG(2,1)=    -AZL(2,1)/RR
            BG(3,1)=    -AZL(3,1)
            BG(3,2)=    -AZL(3,2)/RR
          ELSE IF(ITYP10(nr).EQ.3) THEN
            EG31=AZL(3,1)/(2.0d0*RR)
            EG32=AZL(3,2)/(2.0d0*RR*RC)
            RK1=(AZL(3,3)-RR)/(2.0d0*RR)
            RK2=(AZL(1,3)*EG31+AZL(2,3)*EG32)/2.0d0
            BG(1,1)= RI1-AZL(1,1)
            BG(2,2)=(RI1-AZL(2,2)/RC)/RC
            BG(3,3)=(RI1-AZL(3,3)/RR)/RR
            BG(2,1)=    -AZL(2,1)/RC
            BG(3,1)=    -AZL(3,1)/RR
            BG(3,2)=    -AZL(3,2)/(RR*RC)
          ELSE IF(ITYP10(nr).EQ.4) THEN
            EG31=AZL(3,1)/(2.0d0*G1*G3)
            EG32=AZL(3,2)/(2.0d0*G1*G3)
            RK1=(AZL(3,3)-G3)/(2.0d0*G3)
            RK2=(AZL(1,3)*EG31+AZL(2,3)*EG32)/2.0d0
            BG(1,1)=(RI1-AZL(1,1)/G1)/G1
            BG(2,2)=(RI1-AZL(2,2)/G1)/G1
            BG(3,3)=(RI1-AZL(3,3)/G3)/G3
            BG(2,1)=    -AZL(2,1)/(G1*G1)
            BG(3,1)=    -AZL(3,1)/(G1*G3)
            BG(3,2)=    -AZL(3,2)/(G1*G3)
          ELSE IF(ITYP10(nr).EQ.5) THEN
          ENDIF
          CALL ENERGY(nr,CG,DW,RI1,RI2,RI3,RK1,RK2,0.0d0,YG,ERROR,*9999
     $         )
          IF(KTYP52(nr).EQ.1) THEN
            W3=RI3*DW(3)
          ELSE IF(KTYP52(nr).GT.1) THEN
            IF(KTYP51(nr).EQ.3) THEN
              W3=ZG(NH_LOC(0,nx),1)
            ELSE IF(KTYP51(nr).EQ.4) THEN
              RL33=1.0d0 !!!!!!!!!!!! Check this PJH 15Jan90
              W3=-(DW(1)+DW(2)*BG(3,3))*RL33
            ENDIF
          ENDIF
          DO mj=1,NITB
            DO nj=1,mj
               TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj) +W3
     $              *AZU(mj,nj))
             ENDDO
          ENDDO
          TG(3,1)=TG(3,1)+EG31*DW(5)
          TG(3,2)=TG(3,2)+EG32*DW(5)
          TG(3,3)=TG(3,3)+AXU(3,3)*DW(4)
          TG(1,2)=TG(2,1)
          TG(1,3)=TG(3,1)
          TG(2,3)=TG(3,2)
          IF(KTYP52(nr).EQ.4) THEN
            nh_pressure=NH_LOC(NH_LOC(0,nx),nx)
            TG(1,1)=TG(1,1)-ZG(nh_pressure,1)*AXU(1,1)
            TG(1,2)=TG(1,2)-ZG(nh_pressure,1)*AXU(1,2)
            TG(1,3)=TG(1,3)-ZG(nh_pressure,1)*AXU(1,3)
            TG(2,1)=TG(2,1)-ZG(nh_pressure,1)*AXU(2,1)
            TG(2,2)=TG(2,2)-ZG(nh_pressure,1)*AXU(2,2)
            TG(2,3)=TG(2,3)-ZG(nh_pressure,1)*AXU(2,3)
            TG(3,1)=TG(3,1)-ZG(nh_pressure,1)*AXU(3,1)
            TG(3,2)=TG(3,2)-ZG(nh_pressure,1)*AXU(3,2)
            TG(3,3)=TG(3,3)-ZG(nh_pressure,1)*AXU(3,3)
          ENDIF

        ELSE IF(KTYP53(nr).GT.1) THEN !fibre coords (aeolotropic)
          !Invariants K1,K2 are fns of Physical strain that are
          !isotropic in the plane transverse to the Nu(1) axis
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
          IF(KTYP52(nr).EQ.1) THEN      !compressible
            W3=RI3*DW(3)
          ELSE IF(KTYP52(nr).EQ.6) THEN !comp + fluid for lung
            W3=RI3*DW(3)
C TVK 05/06/2000 Check next line
          ELSE IF(KTYP52(nr).EQ.4) THEN !compressible + fluid
            W3=RI3*DW(3)*ZG(NH_LOC(0,nx),1)
          ELSE IF(KTYP52(nr).GT.1) THEN !incompressible
            IF(KTYP52(nr).NE.5) THEN !not inextensible
              W3=ZG(NH_LOC(0,nx),1) !is pressure term
            ELSE IF(KTYP52(nr).EQ.5) THEN !incomp+inext
              W3=ZG(NH_LOC(0,nx)-1,1) !is hyd pressure term
              W4=ZG(NH_LOC(0,nx),1)   !is fibre stress term
              IF(DOP) THEN
                 WRITE(OP_STRING,'('' W3='',D11.3,'' W4='',D11.3)') W3
     $                ,W4
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF !KTYP52
          ENDIF !KTYP52
          DO mj=1,NITB
            DO nj=1,mj
               TG(mj,nj)=2.0d0*(DW(1)*AXU(mj,nj)+DW(2)*BG(mj,nj) + W3
     $              *AZU(mj,nj))
            ENDDO
          ENDDO
          IF(KTYP52(nr).EQ.5) THEN !incomp+inext
            TG(1,1)=TG(1,1)+W4*AZU(1,1) !is inext fibre stress term
            IF(DOP) THEN
               WRITE(OP_STRING,'('' AZU(1,1)='',D11.3
     $              ,' //''' TG(1,1)='',D11.3)') AZU(1,1),TG(1,1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
C         Note: DW(4) & DW(5) are zero for isotropic case
          TG(3,1)=TG(3,1)+EG13*DW(5)
          TG(2,1)=TG(2,1)+EG12*DW(5)
          TG(1,1)=TG(1,1)+AXU(1,1)*DW(4)
          TG(1,2)=TG(2,1)
          TG(1,3)=TG(3,1)
          TG(2,3)=TG(3,2)
          IF(KTYP52(nr).EQ.4) THEN
            nh_pressure=NH_LOC(NH_LOC(0,nx),nx)
            TG(1,1)=TG(1,1)-ZG(nh_pressure,1)*AXU(1,1)
            TG(1,2)=TG(1,2)-ZG(nh_pressure,1)*AXU(1,2)
            TG(1,3)=TG(1,3)-ZG(nh_pressure,1)*AXU(1,3)
            TG(2,1)=TG(2,1)-ZG(nh_pressure,1)*AXU(2,1)
            TG(2,2)=TG(2,2)-ZG(nh_pressure,1)*AXU(2,2)
            TG(2,3)=TG(2,3)-ZG(nh_pressure,1)*AXU(2,3)
            TG(3,1)=TG(3,1)-ZG(nh_pressure,1)*AXU(3,1)
            TG(3,2)=TG(3,2)-ZG(nh_pressure,1)*AXU(3,2)
            TG(3,3)=TG(3,3)-ZG(nh_pressure,1)*AXU(3,3)
C TVK 31/05/2000 Changes for ktyp52=6
          ELSEIF(KTYP52(nr).EQ.6) THEN
            nj_pressure=NJ_LOC(NJL_FIEL,1,nr)
            TG(1,1)=TG(1,1)-XG(nj_pressure,1)*AXU(1,1)
            TG(1,2)=TG(1,2)-XG(nj_pressure,1)*AXU(1,2)
            TG(1,3)=TG(1,3)-XG(nj_pressure,1)*AXU(1,3)
            TG(2,1)=TG(2,1)-XG(nj_pressure,1)*AXU(2,1)
            TG(2,2)=TG(2,2)-XG(nj_pressure,1)*AXU(2,2)
            TG(2,3)=TG(2,3)-XG(nj_pressure,1)*AXU(2,3)
            TG(3,1)=TG(3,1)-XG(nj_pressure,1)*AXU(3,1)
            TG(3,2)=TG(3,2)-XG(nj_pressure,1)*AXU(3,2)
            TG(3,3)=TG(3,3)-XG(nj_pressure,1)*AXU(3,3)
          ENDIF
        ENDIF
!news  21-MAY-1991 Printing strain invariants for transv isotropy. JSW
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
           WRITE(OP_STRING,'('' RK1= '',D12.4,'' RK2= '',D12.4)') RK1
     $          ,RK2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(KTYP53(nr).EQ.1) THEN
            WRITE(OP_STRING,'('' G3= '',D12.4)')G3
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
CC$        call mp_unsetlock()
        ENDIF !DOP
!newe

C-TG=fn of principal stretches---------------------------------------

      ELSE IF(KTYP55(nr).EQ.2) THEN
C       TG is a function of the princ stretches
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
          ENDDO
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO
        DO i=1,3
          DO j=1,3
            EDG(i,j)=0.5d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO
        ENDDO
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
           WRITE(OP_STRING, '('' EDG: '',9D12.4)') ((EDG(i,j),J=1,3),I=1
     $          ,3)
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
           WRITE(OP_STRING,'('' EVAL:'',7X,(3X,3D11.3))') (EVAL(i),I=1,3
     $          )
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' RM:'',12X,3D11.3,/(16X,3D11.3))') ((RM(i
     $         ,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        RL1=DSQRT(2.0d0*EVAL(1)+1.0d0)
        RL2=DSQRT(2.0d0*EVAL(2)+1.0d0)
        RL3=DSQRT(2.0d0*EVAL(3)+1.0d0)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
           WRITE(OP_STRING,'('' RL1,RL2,RL3:'',(3X,3D11.3))') RL1,RL2
     $          ,RL3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        DO i=1,3
          DO j=1,3
            TGP(i,j)=0.0d0
          ENDDO
        ENDDO
        CALL ENERGY(nr,CG,DW,RL1,RL2,RL3,0.0d0,0.0d0,0.0d0,YG, ERROR,
     $       *9999)
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
     $          ((TGP(i,j),J=1,3),I=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        CALL TRAN(3,TGP,TG,RM)

C-TG=fn of fibre and transv strains----------------------------------

      ELSE IF(KTYP55(nr).EQ.3) THEN !TG=func of fibre and transv strains
C       Compute undeformed covariant metric tensor
        DO i=1,3
          DO j=1,3
            AXL(i,j)=0.0d0
          ENDDO !j
          AXL(i,i)=1.0d0/AXU(i,i)
        ENDDO !i
C!!! NOTE can only use resid strains for pole zero law until init extn
C!!!      mat params are set up for other problem types
        IF(KTYP56(nr).EQ.3) THEN !pole-zero law
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
                   AZL(i,j)=AZL(i,j)+ Fgrowth(i,m)*AZL_tmp(m,n)
     $                  *Fgrowth(n,j)
                ENDDO !n
              ENDDO !m
            ENDDO !j
          ENDDO !i
C         Recompute AZU,AZ from transformed AZL
          CALL INVERT(NIT(nb),AZL,AZU,AZ)
        ENDIF !ktyp56=3 (pole-zero law)
C       Compute Green strain tensor wrt material fibre coords
        DO i=1,3
          DO j=1,3
            EG(i,j)=0.50d0*(AZL(i,j)-AXL(i,j))/DSQRT(AXL(i,i)*AXL(j,j))
          ENDDO !j
        ENDDO !i
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
           WRITE(OP_STRING,'('' EG: '',3D12.4,/(5X,3D12.4))') ((EG(i,j)
     $          ,j=1,3),i=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(KTYP56(nr).EQ.6) THEN !linear viscous relation
C         write(*,'('' linear viscous relation'')')
C         terms for the viscous stress.
C         DD1 is 2nd deriv(wrt nu1)*1st deriv.
C         DD2 is 2nd deriv(wrt nu2)*1st deriv.
          DD1(1,1)=ZG(1,3)*ZG(1,2)+ZG(2,3)*ZG(2,2)+ZG(3,3)*ZG(3,2)
          DD1(2,1)=ZG(1,6)*ZG(1,2)+ZG(2,6)*ZG(2,2)+ZG(3,6)*ZG(3,2)
          DD1(3,1)=ZG(1,9)*ZG(1,2)+ZG(2,9)*ZG(2,2)+ZG(3,9)*ZG(3,2)
          DD1(1,2)=ZG(1,3)*ZG(1,4)+ZG(2,3)*ZG(2,4)+ZG(3,3)*ZG(3,4)
          DD1(2,2)=ZG(1,6)*ZG(1,4)+ZG(2,6)*ZG(2,4)+ZG(3,6)*ZG(3,4)
          DD1(3,2)=ZG(1,9)*ZG(1,4)+ZG(2,9)*ZG(2,4)+ZG(3,9)*ZG(3,4)
          DD1(1,3)=ZG(1,3)*ZG(1,7)+ZG(2,3)*ZG(2,7)+ZG(3,3)*ZG(3,7)
          DD1(2,3)=ZG(1,6)*ZG(1,7)+ZG(2,6)*ZG(2,7)+ZG(3,6)*ZG(3,7)
          DD1(3,3)=ZG(1,9)*ZG(1,7)+ZG(2,9)*ZG(2,7)+ZG(3,9)*ZG(3,7)
C
          DD2(1,1)=ZG(1, 6)*ZG(1,2)+ZG(2, 6)*ZG(2,2)+ZG(3, 6)*ZG(3,2)
          DD2(2,1)=ZG(1, 5)*ZG(1,2)+ZG(2, 5)*ZG(2,2)+ZG(3, 5)*ZG(3,2)
          DD2(3,1)=ZG(1,10)*ZG(1,2)+ZG(2,10)*ZG(2,2)+ZG(3,10)*ZG(3,2)
          DD2(1,2)=ZG(1, 6)*ZG(1,4)+ZG(2, 6)*ZG(2,4)+ZG(3, 6)*ZG(3,4)
          DD2(2,2)=ZG(1, 5)*ZG(1,4)+ZG(2, 5)*ZG(2,4)+ZG(3, 5)*ZG(3,4)
          DD2(3,2)=ZG(1,10)*ZG(1,4)+ZG(2,10)*ZG(2,4)+ZG(3,10)*ZG(3,4)
          DD2(1,3)=ZG(1, 6)*ZG(1,7)+ZG(2, 6)*ZG(2,7)+ZG(3, 6)*ZG(3,7)
          DD2(2,3)=ZG(1, 5)*ZG(1,7)+ZG(2, 5)*ZG(2,7)+ZG(3, 5)*ZG(3,7)
          DD2(3,3)=ZG(1,10)*ZG(1,7)+ZG(2,10)*ZG(2,7)+ZG(3,10)*ZG(3,7)
          fibre_angle=XG(4,1)
!     DO M=1,3 !strain rate term DO N=1,3 D(M,N)=
!     0.5d0*(DD1(M,N)+DD1(N,M))*CG(4)*cos(fibre_angle) ' +
!     0.5d0*(DD2(M,N)+DD2(N,M))*CG(4)*sin(fibre_angle) ENDDO ENDDO
          dNudXi(1,1)=DCOS(fibre_angle)
          dNudXi(2,1)=-DSIN(fibre_angle)
          dNudXi(3,1)=0.d0
          DO m=1,3 !strain rate term - derivs taken wrt Xi1
            DO n=1,3
               D(m,n)=CG(4)*((DD1(m,n)+DD1(n,m))*dNudXi(1,1) +(DD2(m,n)
     $              +DD2(n,m))*dNudXi(2,1))
            ENDDO !n
          ENDDO !m
          IF(KTYP52(nr).EQ.5) THEN !incomp+inext
            hyd_press=ZG(NH_LOC(0,nx)-1,1)
          ELSE
            hyd_press=ZG(NH_LOC(0,nx),1)
          ENDIF
          TG(1,1)=2.d0*(2.d0*CG(2)-CG(3))*D(1,1)+2.0d0*hyd_press
          TG(1,2)=2.d0*CG(2)*D(1,2)
          TG(2,2)=2.d0*CG(3)*D(2,2)+2.0d0*hyd_press
          TG(2,1)=TG(1,2)
          TG(1,3)=2.d0*CG(2)*D(1,3)
          TG(2,3)=2.d0*CG(3)*D(2,3)
          TG(3,1)=2.d0*CG(2)*D(1,3)
          TG(3,2)=2.d0*CG(3)*D(2,3)
          TG(3,3)=2.d0*CG(3)*D(3,3)+2.0d0*hyd_press
            
          IF(KTYP52(nr).EQ.5) THEN !incomp+inext
            W4=ZG(NH_LOC(0,nx),1)  !is fibre stress term
            TG(1,1)=TG(1,1)+W4*AZU(1,1) !is inext fibre stress term
            IF(DOP) THEN
               WRITE(OP_STRING,'('' AZU(1,1)='',D11.3
     $              ,' //''' TG(1,1)='',D11.3)') AZU(1,1),TG(1,1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF !KTYP52(nr)=5
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
             WRITE(OP_STRING,'('' TG: '',3D12.4,/(5X,3D12.4))') ((TG(i,j
     $            ),j=1,3),i=1,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF !DOP

        ELSE IF(KTYP56(nr).NE.6) THEN !all others
           CALL ENERGY(nr,CG,DW,EG(1,1),EG(2,2),EG(3,3),EG(1,2),EG(1,3),
     $          EG(2,3),YG,ERROR,*9999)
          IF(KTYP52(nr).EQ.1) THEN !compressible
            W3=0.0d0 !Check this
          ELSE                     !incompressible
            W3=ZG(NH_LOC(0,nx),1)
          ENDIF
          IF(STRESSTYPE(1:14).EQ.'Total_no_hydro') THEN
            W3=0.0d0
          ENDIF
C news VJ 17Dec2003 Adding an if block for Gauss point stresses (grid coupling) 
          IF(KTYP54(nr).NE.3) THEN !not Gauss point stresses (grid coupling)
            DO i=1,3
              TG(i,i)=DW(i)/AXL(i,i)+2.0d0*W3*AZU(i,i)
            ENDDO
            TG(1,2)=DW(4)/DSQRT(AXL(1,1)*AXL(2,2))+2.0d0*W3*AZU(1,2)
            TG(1,3)=DW(5)/DSQRT(AXL(1,1)*AXL(3,3))+2.0d0*W3*AZU(1,3)
            TG(2,3)=DW(6)/DSQRT(AXL(2,2)*AXL(3,3))+2.0d0*W3*AZU(2,3)
            TG(2,1)=TG(1,2)
            TG(3,1)=TG(1,3)
            TG(3,2)=TG(2,3)
C following is the added code. problem is the DW array stores the deviatoric components of
C the strain energy function in alpha numeric order. But in CMISS it is 11,22,33 then 12,13,23
C delete this else if block and change the code to be more general
          ELSE IF(KTYP54(nr).EQ.3) THEN !Gauss point stresses (grid coupling)
            TG(1,1)=DW(1)/AXL(1,1)+2.0d0*W3*AZU(1,1)
            TG(2,2)=DW(4)/AXL(2,2)+2.0d0*W3*AZU(2,2)
            TG(3,3)=DW(6)/AXL(3,3)+2.0d0*W3*AZU(3,3)
C Refer to Holger Schmid.et.al technical note on definition of stress in terms of 
C derivatives of the strain energy function. "Ambiguities in Hyperelastic Constitutive Law Formulation
C The formula for the stress components of a material depend on the formulation of the strain energy
C density function. If the strain energy function may be expressed in terms of 6 or 9 different strain
C components - W(E11,E22,E33,E12,E13,E21,E23) or W(E11,E22,E33,E12,E13,E23,E21,E32,E31). If expressed in terms
C of 6 components then 
C TGij deviatoric = dWdEij for i=j i=1...3, j=1...3
C TGij deviatoric = 0.5dWdEij for i<>j i=1...3, j=1...3
C else if expressed in terms of 9 components
C TGij deviatoric = dWdEij for i,j..i=1...3, j=1...3
C Hence it has been decided, in order for people to have the freedom to express W in their own way, DW array
C will contain the deviatoric component of stress calculated in the CellML file. 
            TG(1,2)=DW(2)/DSQRT(AXL(1,1)*AXL(2,2))+ 2.0d0*W3*AZU(1,2)
            TG(1,3)=DW(3)/DSQRT(AXL(1,1)*AXL(3,3))+ 2.0d0*W3*AZU(1,3)
            TG(2,3)=DW(5)/DSQRT(AXL(2,2)*AXL(3,3))+ 2.0d0*W3*AZU(2,3)
            TG(2,1)=TG(1,2)
            TG(3,1)=TG(1,3)
            TG(3,2)=TG(2,3)
          ENDIF !ktyp54
C newe VJ
        ENDIF !ktyp56

      ENDIF !ktyp55

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
         WRITE(OP_STRING,'('' TG:'',12X,3D12.4,/(16X,3D12.4))') ((TG(i,j
     $        ),j=1,3),i=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZGTG53')
      RETURN
 9999 CALL ERRORS('ZGTG53',ERROR)
      CALL EXITS('ZGTG53')
      RETURN 1
      END


