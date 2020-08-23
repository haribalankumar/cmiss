      SUBROUTINE MG_RESIDUAL(na,niq,NITB,NAQ,NLQ,NWQ,NXQ,
     '  ADAPTIVE,ISOTROPIC,TYPE,DT,GCHQ,GUQ,PROPQ,YQ,ERROR,*)

C#### Subroutine: MG_RESIDUAL
C###  Description:
C###    MG_RESIDUAL is the multigrid residual operator. It returns
C###    the residual of YQ(nq,1,na) in YQ(nq,niq,na) for given level na.
C###    na=1 for fine grid solution, =2,3.. for coarse grid solution.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER na,niq,NITB,
     '  NAQ(NQM),NLQ(NQM),NWQ(8,0:NQM),NXQ(-NIM:NIM,0:4,0:NQM)
      REAL*8 DT,GCHQ(3,NQM),GUQ(3,3,NQM),PROPQ(3,3,4,2,NQM),
     '  YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL ADAPTIVE,ISOTROPIC
!     Local Variables
      INTEGER i,ib,IJ,IK,j,jb,k,mq,ni,nii,nij,nik,nj,nq
      REAL*8 AA,D,DD,dUdX(3),d2UdX2(3,3),
     '  GCHQ1,GCHQ2,GUQ11,GUQ22,GUQ12,IPROP,SUM1,SUM2,SUM3,
     '  U,UQ(-1:1,-1:1,-1:1)

      CALL ENTERS('MG_RESIDUAL',*9999)

      DO i=-1,1
        DO j=-1,1
          DO k=-1,1
            UQ(i,j,k)=0.0d0
          ENDDO
        ENDDO
      ENDDO

C     DD=(0.5*DBLE(2**na)/DBLE(2**nat(1)))**2 !grid spacing
      D=2.d0/DBLE(2**na)
      DD=D*D
      DO nq=1,NQT !Loop over grid points
        IF(NAQ(nq).EQ.0         !nq is active in grid na
     '    .AND.(.NOT.ADAPTIVE
     '          .OR.(ADAPTIVE.AND.NLQ(nq).EQ.na))) THEN

          IF(NWQ(1,nq).EQ.0) THEN !interior g.p.
            IK=MAX(0,NITB-2) !zero for 1,2D, one for 3D
            IJ=MIN(NITB-1,1) !zero for 1D, one for 2,3D
            DO nik=-IK,IK
              DO nij=-IJ,IJ
                DO nii=-1,1
                  mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq))) !adj g.p.
                  IF(mq.GT.0) THEN
                    UQ(nii,nij,nik)=YQ(mq,1,na)       !u at mq
                  ELSE
                    write(*,'('' Error: nq='',I6,'
     '                //''' nii='',I2,'' nij='',I2,'' nik='',I2,'
     '                //''' outside grid in MG_RESIDUAL'')')
     '                nq,nii,nij,nik
                  ENDIF
               ENDDO !nii
             ENDDO !nij
           ENDDO !nik

C         Compute residual
            IF(ISOTROPIC) THEN !2D isotropy
              GUQ11=GUQ(1,1,nq)*DD
              GUQ22=GUQ(2,2,nq)*DD
              GUQ12=GUQ(1,2,nq)*DD
              GCHQ1=GCHQ( 1,nq)*D
              GCHQ2=GCHQ( 2,nq)*D
              AA=4.0D0*GUQ11*(UQ(-1, 0, 0)+UQ( 1, 0, 0))
     '          +4.0D0*GUQ22*(UQ( 0,-1, 0)+UQ( 0, 1, 0))
     '          +2.0D0*GUQ12*(UQ(-1,-1, 0)+UQ( 1, 1, 0)
     '                       -UQ( 1,-1, 0)-UQ(-1, 1, 0))
     '                +GCHQ1*(UQ(-1, 0, 0)-UQ( 1, 0, 0))
     '                +GCHQ2*(UQ( 0,-1, 0)-UQ( 0, 1, 0))
              YQ(nq,niq,na)=8.0D0*(GUQ11+GUQ22)*UQ(0,0,0)-
     '                             (YQ(nq,5,na)+AA)
c rect. grid  YQ(nq,niq,na)=(4.0D0*UQ(0,0,0)-UQ(-1,0,0)-UQ(1,0,0)
c    '                             -UQ(0,-1,0)-UQ(0,1,0))/DD

            ELSE IF(.NOT.ISOTROPIC) THEN !general anisotropic
C           Compute U,k by 1st f.d.s about nq
              dUdX(1)=UQ(1,0,0)-UQ(-1,0,0)
              dUdX(2)=UQ(0,1,0)-UQ(0,-1,0)
              dUdX(3)=UQ(0,0,1)-UQ(0,0,-1)
              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RESIDUAL_1)
                WRITE(OP_STRING,'('' du/dxi = '',3F10.5)')
     '            (dUdX(ni),ni=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_RESIDUAL_1)
              ENDIF
C           Compute U,ij=d2U/dXi(i)dXi(j) by taking 2nd f.d.s about nq
              DO ib=1,NITB
                DO jb=1,NITB
                  IF(ib.EQ.jb) THEN
                    IF(ib.EQ.1) THEN !U,11
                      d2UdX2(ib,jb)=(UQ( 1,0,0)-2.0d0*UQ(0,0,0)+
     '                               UQ(-1,0,0))*4.0d0
                    ELSE IF(ib.EQ.2) THEN !U,22
                      d2UdX2(ib,jb)=(UQ(0, 1,0)-2.0d0*UQ(0,0,0)+
     '                               UQ(0,-1,0))*4.0d0
                    ELSE IF(ib.EQ.3) THEN !U,33
                      d2UdX2(ib,jb)=(UQ(0,0, 1)-2.0d0*UQ(0,0,0)+
     '                               UQ(0,0,-1))*4.0d0
                    ENDIF
                  ELSE
                    IF(ib+jb.EQ.3) THEN !U,12 and U,21
                      d2UdX2(ib,jb)= UQ( 1,1,0)-UQ( 1,-1,0)
     '                              -UQ(-1,1,0)+UQ(-1,-1,0)
                    ELSE IF(ib+jb.EQ.4) THEN !U,13 and U,31
                      d2UdX2(ib,jb)= UQ( 1,0,1)-UQ( 1,0,-1)
     '                              -UQ(-1,0,1)+UQ(-1,0,-1)
                    ELSE IF(ib+jb.EQ.5) THEN !U,23 and U,32
                      d2UdX2(ib,jb)= UQ(0,1, 1)-UQ(0,-1, 1)
     '                              -UQ(0,1,-1)+UQ(0,-1,-1)
                    ENDIF
                  ENDIF !ib=jb
                ENDDO !jb
              ENDDO !ib

              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RESIDUAL_2)
                WRITE(OP_STRING,'('' d2U/dXi2'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                DO nj=1,NITB
                  WRITE(OP_STRING,'(3F10.5)') (d2UdX2(ni,nj),ni=1,nitb)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO
CC$OMP END CRITICAL(MG_RESIDUAL_2)
              ENDIF !dop

C***  Compute residual YQ(nq,niq,na)=IPROP
              IPROP=0.0d0
              DO i=1,NITB
                DO j=1,NITB
                  SUM1=0.0d0
                  SUM2=0.0d0
                  DO k=1,NITB
                    SUM1=SUM1+PROPQ(k,i,j+1,1,nq)*dUdX(k)
                    SUM2=SUM2+PROPQ(k,i,1,1,nq)*d2UdX2(j,k)
                  ENDDO
                  IPROP=IPROP-(SUM1*D+SUM2*DD)*GUQ(i,j,nq)
                ENDDO !j
              ENDDO !i
              DO j=1,NITB
                SUM3=0.0d0
                DO i=1,NITB
                  SUM3=SUM3+PROPQ(i,j,1,1,nq)*dUdX(i)
                ENDDO !i
                IPROP=IPROP+SUM3*D*GCHQ(j,nq)
              ENDDO !j

              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RESIDUAL_3)
                WRITE(OP_STRING,'('' Iprop = '',F15.5,'' uA/mm^3'')')
     '            IPROP*1.d3
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_RESIDUAL_3)
              ENDIF !dop

              IF(TYPE(1:6).EQ.'STATIC') THEN
                YQ(nq,niq,na)=IPROP-YQ(nq,5,na)
              ELSE IF(TYPE(1:7).EQ.'DYNAMIC') THEN
C               (YQ(6) holds previous timestep solution)
                U=UQ(0,0,0)
                YQ(nq,niq,na)=IPROP-YQ(nq,5,na)+(U-YQ(nq,6,na))/DT
              ELSE IF(TYPE(1:6).EQ.'ACTIVE') THEN
                U=UQ(0,0,0)
                YQ(nq,niq,na)=IPROP-YQ(nq,5,na)
     '            + (U-YQ(nq,6,na))/DT + U*(1.d0-U)*(0.15d0-U)
              ENDIF !type

            ENDIF !isotropic/anisotropic

          ELSE !Boundary point
            IF(NWQ(5,nq).EQ.1) THEN      !Dirichlet b.c.
              YQ(nq,niq,na)=0.0D0
            ELSE IF(NWQ(5,nq).EQ.2) THEN !Neumann b.c.
              IF(NWQ(2,nq).GT.0) THEN !2 interior pts defined
                YQ(nq,niq,na)=3.0D0*YQ(nq,1,na)+YQ(NWQ(2,nq),1,na)
     '                        -4.0D0*YQ(NWQ(1,nq),1,na)
              ELSE                       !only 1 interior pt defined
                YQ(nq,niq,na)=YQ(nq,1,na)-YQ(NWQ(1,nq),1,na)
              ENDIF
            ENDIF
          ENDIF !NWQ

          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RESIDUAL_4)
            WRITE(OP_STRING,
     '       '('' Residual YQ('',I6,'','',I1,'','',I1,'')='',D12.4)')
     '        nq,niq,na,YQ(nq,niq,na)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_RESIDUAL_4)
          ENDIF

        ENDIF !NAQ=0
      ENDDO !Global grid point nq

      CALL EXITS('MG_RESIDUAL')
      RETURN
 9999 CALL ERRORS('MG_RESIDUAL',ERROR)
      CALL EXITS('MG_RESIDUAL')
      RETURN 1
      END


